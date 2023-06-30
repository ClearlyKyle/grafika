#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "grafika.h"
#include "tex.h"
#include "obj.h"

#define max(a, b) (((a) > (b)) ? (a) : (b))
#define min(a, b) (((a) < (b)) ? (a) : (b))

typedef float vec3[3];
typedef float vec4[4];
typedef vec4  mat4[4];
typedef vec3  triangle[3];

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

static inline float v3dot(const vec3 v0, const vec3 v1)
{
    return (v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2]);
}

static inline void v3sub(const vec3 v0, const vec3 v1, vec3 res)
{
    res[0] = v0[0] - v1[0];
    res[1] = v0[1] - v1[1];
    res[2] = v0[2] - v1[2];
}

static inline void v3div(const vec3 v, const float val, vec3 res)
{
    res[0] = v[0] / val;
    res[1] = v[1] / val;
    res[2] = v[2] / val;
}

static inline float edgefunc(const vec3 a, const vec3 b, const vec3 c)
{
    return (c[0] - a[0]) * (b[1] - a[1]) - (c[1] - a[1]) * (b[0] - a[0]);
}

static inline void m4mulm4(const mat4 m1, const mat4 m2, mat4 dest)
{
    float a00 = m1[0][0], a01 = m1[0][1], a02 = m1[0][2], a03 = m1[0][3],
          a10 = m1[1][0], a11 = m1[1][1], a12 = m1[1][2], a13 = m1[1][3],
          a20 = m1[2][0], a21 = m1[2][1], a22 = m1[2][2], a23 = m1[2][3],
          a30 = m1[3][0], a31 = m1[3][1], a32 = m1[3][2], a33 = m1[3][3],

          b00 = m2[0][0], b01 = m2[0][1], b02 = m2[0][2], b03 = m2[0][3],
          b10 = m2[1][0], b11 = m2[1][1], b12 = m2[1][2], b13 = m2[1][3],
          b20 = m2[2][0], b21 = m2[2][1], b22 = m2[2][2], b23 = m2[2][3],
          b30 = m2[3][0], b31 = m2[3][1], b32 = m2[3][2], b33 = m2[3][3];

    dest[0][0] = a00 * b00 + a10 * b01 + a20 * b02 + a30 * b03;
    dest[0][1] = a01 * b00 + a11 * b01 + a21 * b02 + a31 * b03;
    dest[0][2] = a02 * b00 + a12 * b01 + a22 * b02 + a32 * b03;
    dest[0][3] = a03 * b00 + a13 * b01 + a23 * b02 + a33 * b03;
    dest[1][0] = a00 * b10 + a10 * b11 + a20 * b12 + a30 * b13;
    dest[1][1] = a01 * b10 + a11 * b11 + a21 * b12 + a31 * b13;
    dest[1][2] = a02 * b10 + a12 * b11 + a22 * b12 + a32 * b13;
    dest[1][3] = a03 * b10 + a13 * b11 + a23 * b12 + a33 * b13;
    dest[2][0] = a00 * b20 + a10 * b21 + a20 * b22 + a30 * b23;
    dest[2][1] = a01 * b20 + a11 * b21 + a21 * b22 + a31 * b23;
    dest[2][2] = a02 * b20 + a12 * b21 + a22 * b22 + a32 * b23;
    dest[2][3] = a03 * b20 + a13 * b21 + a23 * b22 + a33 * b23;
    dest[3][0] = a00 * b30 + a10 * b31 + a20 * b32 + a30 * b33;
    dest[3][1] = a01 * b30 + a11 * b31 + a21 * b32 + a31 * b33;
    dest[3][2] = a02 * b30 + a12 * b31 + a22 * b32 + a32 * b33;
    dest[3][3] = a03 * b30 + a13 * b31 + a23 * b32 + a33 * b33;
}

static inline void m4mulv4(const mat4 m, const vec4 v, vec4 res)
{
    res[0] = m[0][0] * v[0] + m[1][0] * v[1] + m[2][0] * v[2] + m[3][0] * v[3];
    res[1] = m[0][1] * v[0] + m[1][1] * v[1] + m[2][1] * v[2] + m[3][1] * v[3];
    res[2] = m[0][2] * v[0] + m[1][2] * v[1] + m[2][2] * v[2] + m[3][2] * v[3];
    res[3] = m[0][3] * v[0] + m[1][3] * v[1] + m[2][3] * v[2] + m[3][3] * v[3];
}

void normalize(vec3 v)
{
    float length = sqrtf(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);

    if (length == 0.0f)
    {
        v[0] = v[1] = v[2] = 0.0f;
        return;
    }

    v[0] /= length;
    v[1] /= length;
    v[2] /= length;
}

void cross(const vec3 v1, const vec3 v2, vec3 res)
{
    res[0] = v1[1] * v2[2] - v1[2] * v2[1];
    res[1] = v1[2] * v2[0] - v1[0] * v2[2];
    res[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

void subtract(const vec3 v1, const vec3 v2, vec3 res)
{
    res[0] = v1[0] - v2[0];
    res[1] = v1[1] - v2[1];
    res[2] = v1[2] - v2[2];
}

void v3add(const vec3 v1, const vec3 v2, vec3 res)
{
    res[0] = v1[0] + v2[0];
    res[1] = v1[1] + v2[1];
    res[2] = v1[2] + v2[2];
}

static inline void m4identity(mat4 m)
{
    memset(m, 0, sizeof(mat4));
    m[0][0] = m[1][1] = m[2][2] = m[3][3] = 1.0f;
}

void rotateMatrix(mat4 matrix, float angle, char axis)
{
    m4identity(matrix);

    const float c = cosf(angle);
    const float s = sinf(angle);

    switch (axis)
    {
    case 'x':
    case 'X':
        matrix[1][1] = c;
        matrix[1][2] = -s;
        matrix[2][1] = s;
        matrix[2][2] = c;
        break;

    case 'y':
    case 'Y':
        matrix[0][0] = c;
        matrix[0][2] = s;
        matrix[2][0] = -s;
        matrix[2][2] = c;
        break;

    case 'z':
    case 'Z':
        matrix[0][0] = c;
        matrix[0][1] = -s;
        matrix[1][0] = s;
        matrix[1][1] = c;
        break;

    default:
        // Invalid axis specified
        break;
    }
}

void m4lookat(const vec3 eye, const vec3 dir, const vec3 up, mat4 res)
{
    vec3 target;
    v3add(eye, dir, target);

    vec3 forward, right, new_up;

    subtract(target, eye, forward);
    normalize(forward);

    cross(forward, up, right);
    normalize(right);

    cross(right, forward, new_up);
    // normalize(new_up);

    res[0][0] = right[0];
    res[0][1] = new_up[0];
    res[0][2] = -forward[0];
    res[0][3] = 0.0f;

    res[1][0] = right[1];
    res[1][1] = new_up[1];
    res[1][2] = -forward[1];
    res[1][3] = 0.0f;

    res[2][0] = right[2];
    res[2][1] = new_up[2];
    res[2][2] = -forward[2];
    res[2][3] = 0.0f;

    res[3][0] = -v3dot(right, eye);
    res[3][1] = -v3dot(new_up, eye);
    res[3][2] = v3dot(forward, eye);
    res[3][3] = 1.0f;
}

static inline void m4transmake(float x, float y, float z, mat4 res)
{
    m4identity(res);
    res[3][0] = x;
    res[3][1] = y;
    res[3][2] = z;
}

static inline void m4scalemake(float x, float y, float z, mat4 res)
{
    m4identity(res);
    res[0][0] = x;
    res[1][1] = y;
    res[2][2] = z;
}

static inline void m4perspective(float fovy, float aspect, float nearZ, float farZ, mat4 res)
{
    memset(res, 0, sizeof(mat4));

    const float f  = 1.0f / tanf(fovy * 0.5f);
    const float fn = 1.0f / (nearZ - farZ);

#if 1 // [-1, 1]
    res[0][0] = f / aspect;
    res[1][1] = f;
    res[2][2] = (nearZ + farZ) * fn;
    res[2][3] = -1.0f;
    res[3][2] = 2.0f * nearZ * farZ * fn;
#else // [0, 1]
    res[0][0] = f / aspect;
    res[1][1] = f;
    res[2][2] = farZ * fn;
    res[2][3] = -1.0f;
    res[3][2] = nearZ * farZ * fn;
#endif
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

#define _PI          3.14159265358979323846264338327950288 /* pi */
#define _PIf         ((float)_PI)
#define DEG2RAD(DEG) ((DEG)*_PIf / 180.0f)

static inline void update(void)
{
    mat4 MVP;
    m4mulm4(state.view, state.model, MVP);
    m4mulm4(state.proj, MVP, MVP);

    obj_t obj = state.obj;
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

    state.obj = obj_load("res/cube.obj");
    // state.obj = obj_load("plane.obj");
    // state.obj = obj_load("bunny.obj");
    // state.obj = obj_load("../res/Dog House/Doghouse.obj");

    // state.tex = tex_load("../res/Dog House/Doghouse_PBR_BaseColor.png", true);
    state.tex = tex_load("res/wood.png", false);
    //  state.tex = tex_load("metal.png", false);

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
                // if (scrollAmount > -2.0f)
                //     scrollAmount = -2.0f;
            }
        }
        // printf("tranz : %f\n", transz);

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
    }

    obj_destroy(&state.obj);
    grafika_destroy();

    return 0;
}