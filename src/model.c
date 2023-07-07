#include <stdio.h>
#include <stdlib.h>

#include "grafika.h"
#include "tex.h"
#include "obj.h"
#include "timer.h"
#include "matematika.h"

typedef vec3 triangle[3];

typedef struct rasterstate
{
    mat4 proj;
    mat4 view;
    mat4 model;

    obj_t obj;
    tex_t tex;
} rasterstate_t;

static rasterstate_t state = {0};

static float zbuffer[GRAFIKA_SCREEN_HEIGHT * GRAFIKA_SCREEN_WIDTH];

static inline float edgefunc(const vec3 a, const vec3 b, const vec3 c)
{
    return (c[0] - a[0]) * (b[1] - a[1]) - (c[1] - a[1]) * (b[0] - a[0]);
}

void drawtriangle(const triangle t, const vec3 texcoords[3], const mat4 MVP)
{
    vec3 screenspace[3] = {0};
    vec4 clipspace[3]   = {0};
    vec3 ndc[3]         = {0};

    // convert to clip space
    for (int i = 0; i < 3; ++i)
    {
        //  clip space position
        m4mulv4(MVP, (vec4){t[i][0], t[i][1], t[i][2], 1.0f}, clipspace[i]);

        // clipping (is this correct?)
        const float x = fabsf(clipspace[i][0]);
        const float y = fabsf(clipspace[i][1]);
        const float w = fabsf(clipspace[i][3]);

        if ((-w <= x && x <= w) || (-w <= y && y <= w))
            continue;
        else
            return;
    }

    // perspective division (clip to ndc)
    for (size_t i = 0; i < 3; i++)
    {
        ndc[i][0] = clipspace[i][0] / clipspace[i][3];
        ndc[i][1] = clipspace[i][1] / clipspace[i][3];
        ndc[i][2] = clipspace[i][2] / clipspace[i][3];
    }

    // back face culling (surface normal, can this be done with area?)
    // vec3 sub10, sub20, normal;
    // v3sub(ndc[1], ndc[0], sub10);
    // v3sub(ndc[2], ndc[0], sub20);
    // cross(sub10, sub20, normal);
    // if (normal[2] > 0.0f)
    //    return;

    for (int i = 0; i < 3; ++i)
    {
        screenspace[i][0] = (ndc[i][0] + 1.0f) * 0.5f * GRAFIKA_SCREEN_WIDTH;
        screenspace[i][1] = (ndc[i][1] + 1.0f) * 0.5f * GRAFIKA_SCREEN_HEIGHT;
        screenspace[i][2] = (ndc[i][2] + 1.0f) * 0.5f;
    }

    // calculate bounding rectangle
    float fminX = fminf(screenspace[0][0], fminf(screenspace[1][0], screenspace[2][0]));
    float fminY = fminf(screenspace[0][1], fminf(screenspace[1][1], screenspace[2][1]));
    float fmaxX = fmaxf(screenspace[0][0], fmaxf(screenspace[1][0], screenspace[2][0]));
    float fmaxY = fmaxf(screenspace[0][1], fmaxf(screenspace[1][1], screenspace[2][1]));

    // clip to screen space
    int minX = max(0, min((int)floorf(fminX), GRAFIKA_SCREEN_WIDTH - 1));
    int minY = max(0, min((int)floorf(fminY), GRAFIKA_SCREEN_HEIGHT - 1));
    int maxX = max(0, min((int)floorf(fmaxX), GRAFIKA_SCREEN_WIDTH - 1));
    int maxY = max(0, min((int)floorf(fmaxY), GRAFIKA_SCREEN_HEIGHT - 1));

    vec3 col1, col2, col3;
    v3div((vec3){1.0f, 0.0f, 0.0f}, clipspace[0][3], col1);
    v3div((vec3){0.0f, 1.0f, 0.0f}, clipspace[1][3], col2);
    v3div((vec3){0.0f, 0.0f, 1.0f}, clipspace[2][3], col3);

    vec3 uv[3];
    v3div(texcoords[0], clipspace[0][3], uv[0]);
    v3div(texcoords[1], clipspace[1][3], uv[1]);
    v3div(texcoords[2], clipspace[2][3], uv[2]);

    // Rasterize
    triangle tri = {
        {screenspace[0][0], screenspace[0][1], 1.0f / clipspace[0][3]},
        {screenspace[1][0], screenspace[1][1], 1.0f / clipspace[1][3]},
        {screenspace[2][0], screenspace[2][1], 1.0f / clipspace[2][3]},
    };

    float area     = edgefunc(tri[0], tri[1], tri[2]);
    float inv_area = 1.0f / area;

    for (int i = minY; i <= maxY; ++i)
    {
        for (int j = minX; j <= maxX; ++j)
        {
            vec3 point = {0.5f + (float)j, 0.5f + (float)i, 0.0f};

            float w0 = edgefunc(tri[1], tri[2], point);
            float w1 = edgefunc(tri[2], tri[0], point);
            float w2 = edgefunc(tri[0], tri[1], point);

            if (w0 < 0.0f || w1 < 0.0f || w2 < 0.0f)
                continue;

            w0 *= inv_area;
            w1 *= inv_area;
            w2 *= inv_area;

            // correction factor
            float cf = 1.0f / (w0 * tri[0][2] + w1 * tri[1][2] + w2 * tri[2][2]);

            const int   index = (i * GRAFIKA_SCREEN_WIDTH) + j;
            const float newZ  = w0 * screenspace[0][2] + w1 * screenspace[1][2] + w2 * screenspace[2][2];
            const float oldZ  = zbuffer[index];
            if (newZ < oldZ)
                continue;

            zbuffer[index] = newZ;

            float u = (uv[0][0] * w0 + uv[1][0] * w1 + uv[2][0] * w2) * cf;
            float v = (uv[0][1] * w0 + uv[1][1] * w1 + uv[2][1] * w2) * cf;

            u = fabsf(u);
            v = fabsf(v);

            u *= (float)state.tex.w - 1;
            v *= (float)state.tex.h - 1;

            unsigned char *texcolour   = state.tex.data + (((int)u + state.tex.w * (int)v) * state.tex.bpp);
            uint32_t       pixelcolour = (0xFF << 24) + (texcolour[2] << 16) + (texcolour[1] << 8) + (texcolour[0] << 0);

            // Draw pixel
            // 0xFF 00 00 00
            // uint32_t pixelcolour = (0xFF << 24) + (red << 16) + (blu << 8) + (gre << 0);
            grafika_setpixel(j, i, pixelcolour);
        }
    }
}

static inline void update(void)
{
    mat4 MVP;
    m4mulm4(state.view, state.model, MVP);
    m4mulm4(state.proj, MVP, MVP);

    obj_t obj = state.obj;

#pragma omp parallel for shared(obj)
    for (int i = 0; i < obj.num_f_rows; i++)
    {
        triangle t     = {0};
        vec3     tc[3] = {0};
        for (int j = 0; j < 3; j++)
        {
            vertindices_t indices = obj.indices[i * 3 + j];

            const int vertexIndex = indices.v_idx;
            const int texIndex    = indices.vt_idx;

            t[j][0] = obj.pos[3 * vertexIndex + 0];
            t[j][1] = obj.pos[3 * vertexIndex + 1];
            t[j][2] = obj.pos[3 * vertexIndex + 2];

            tc[j][0] = obj.texs[2 * texIndex + 0];
            tc[j][1] = obj.texs[2 * texIndex + 1];
        }

        drawtriangle(t, tc, MVP);
    }
}

int main(int argc, char *argv[])
{
    grafika_startup();

    m4perspective(DEG2RAD(60.0f), (float)GRAFIKA_SCREEN_WIDTH / (float)GRAFIKA_SCREEN_HEIGHT, 0.1f, 100.0f, state.proj);

    // state.obj = obj_load("res/cube.obj");
    //  state.obj = obj_load("plane.obj");
    //  state.obj = obj_load("bunny.obj");
    // state.obj = obj_load("res/Dog House/Doghouse.obj");
    // state.obj = obj_load("res/Lorry/lorry.obj");
    state.obj = obj_load("res/Camera/Camera.obj");

    // state.tex = tex_load("res/Dog House/Doghouse_PBR_BaseColor.png", true);
    // state.tex = tex_load("res/Lorry/texture_atlas.png", true);
    state.tex = tex_load("res/Camera/texture_atlas.png", true);
    // state.tex = tex_load("res/wood.png", false);
    //    state.tex = tex_load("metal.png", false);

    // Calculate model height
    BoundingBox_t bbox    = state.obj.bbox;
    float         centerx = -(bbox.min[0] + bbox.max[0]) * 0.5f;
    float         centery = -(bbox.min[1] + bbox.max[1]) * 0.5f;
    float         centerz = -(bbox.min[2] + bbox.max[2]) * 0.5f;

    mat4 trans;
    m4transmake(centerx, centery, centerz, trans);

    // Scale model height to 1
    float model_scale = 1.0f / (bbox.max[1] - bbox.min[1]);

    mat4 scale;
    m4scalemake(model_scale, model_scale, model_scale, scale);

    float rotationAngleX = 0.0f, rotationAngleY = 0.0f;
    float scrollAmount = -2.0f;

    timer_t frame_timer = {0};
    TIMER_START(frame_timer);

    int    frame_counter          = 0;
    double frame_time_accumulated = 0.0;

    while (!rend.quit)
    {
        for (int i = 0; i < GRAFIKA_SCREEN_WIDTH * GRAFIKA_SCREEN_HEIGHT; ++i)
            zbuffer[i] = 1.0f;

        SDL_Event event;
        while (SDL_PollEvent(&event))
        {
            if (SDL_QUIT == event.type)
            {
                rend.quit = true;
            }

            if (event.type == SDL_MOUSEMOTION && (event.motion.state & SDL_BUTTON(SDL_BUTTON_LEFT)))
            {
                rotationAngleY += event.motion.xrel * 0.4f;
                rotationAngleX += event.motion.yrel * 0.4f;
            }

            if (event.type == SDL_MOUSEWHEEL)
            {
                scrollAmount += (float)event.wheel.y;
                if (scrollAmount > -1.5f)
                    scrollAmount = -1.5f;
            }
        }

        vec3 eye    = {0.0f, 0.0f, scrollAmount};
        vec3 target = {0.0f, 0.0f, -1.0f};
        vec3 up     = {0.0f, -1.0f, 0.0f};
        m4lookat(eye, target, up, state.view);

        // Create rotation matrices
        mat4 rotx, roty;
        rotateMatrix(rotx, DEG2RAD(rotationAngleX), 'x');
        rotateMatrix(roty, DEG2RAD(rotationAngleY), 'y');

        // Combine rotations by multiplying matrices (SRT)
        mat4 rot;
        m4mulm4(rotx, roty, rot);
        m4mulm4(rot, trans, state.model);
        m4mulm4(scale, state.model, state.model);

        grafika_clear();
        update();
        grafika_present();

        TIMER_UPDATE(frame_timer);

        if (frame_counter == 64)
        {
            const double avg_time = frame_time_accumulated / frame_counter;

            char buff[16] = {0};
            sprintf_s(buff, sizeof(buff), "%0.2fms", avg_time);
            SDL_SetWindowTitle(rend.window, buff);

            frame_counter          = 0;
            frame_time_accumulated = 0.0;
        }
        else
        {
            frame_time_accumulated += TIMER_ELAPSED_MS(frame_timer);
            frame_counter++;
        }
    }

    obj_destroy(&state.obj);
    grafika_destroy();

    return 0;
}