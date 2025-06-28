#define GRAFIKA_TITLE ("grafika - edging")
#include "common.h"

struct triangle
{
    vec3 pos[3];
    vec2 tex[3];
};

struct edging
{
    mat4             model;
    mat4             MVP;
    vec3             cam_pos;
    struct tex       tex;
    struct obj       obj;
    struct triangle *triangles;
};

struct edging raster_state = {0};

static void draw_triangle(const struct triangle t)
{
    // convert to clip space
    vec4 clip_space[3];
    for (int i = 0; i < 3; ++i)
    {
        vec4 pos = {t.pos[i][0], t.pos[i][1], t.pos[i][2], 1.0f};
        m4_mul_v4(raster_state.MVP, pos, clip_space[i]);

        // clipping
        if (clip_space[i][0] < -clip_space[i][3] || clip_space[i][0] > clip_space[i][3] ||
            clip_space[i][1] < -clip_space[i][3] || clip_space[i][1] > clip_space[i][3] ||
            clip_space[i][2] < -clip_space[i][3] || clip_space[i][2] > clip_space[i][3])
        {
            return; // triangle is outside
        }
    }

    // perspective division (clip to ndc)
    vec3 ndc[3] = {0}, w_vals = {0};
    for (int i = 0; i < 3; ++i)
    {
        w_vals[i] = 1.0f / clip_space[i][3]; // 1.0f / w
        ndc[i][0] = clip_space[i][0] * w_vals[i];
        ndc[i][1] = clip_space[i][1] * w_vals[i];
        ndc[i][2] = clip_space[i][2] * w_vals[i];
    }

    // NOTE : something happens here, using (ndc[i][1] + 1.0f) doesn't cuase a triangle winding switch
    //  but using (1.0f - ndc[i][1]) will flip CCW triangles into CW, and thus invert our object
    // NOTE : we are using a pixel buffer given to us by SDL, and SDL has its screen space origin at
    //  the "Top Left", and since we are doing it this way, then in our raster loop, we are dealing with CW
    //  triangles :O, SDL pixel buffer origin is top-left with Y going down
    vec3 screen_space[3];
    for (int i = 0; i < 3; ++i)
    {
        screen_space[i][0] = (ndc[i][0] + 1.0f) * 0.5f * GRAFIKA_SCREEN_WIDTH;
        // screen_space[i][1] = (ndc[i][1] + 1.0f) * 0.5f * GRAFIKA_SCREEN_HEIGHT; // origin Bottom Left
        screen_space[i][1] = (1.0f - ndc[i][1]) * 0.5f * GRAFIKA_SCREEN_HEIGHT; // origin Top Left
        screen_space[i][2] = (ndc[i][2] + 1.0f) * 0.5f;
    }

    // back face culling with signed area
    // positive area - triangle is CCW
    // negative area - triangle is CW
    float area = v3_edgefunc(screen_space[0], screen_space[1], screen_space[2]);
    if (area >= 0.0f) return;

    // calculate bounding rectangle
    int AABB[4] = {0};
    AABB_make(screen_space, AABB);

    float inv_area = 1.0f / area;

    // pre fetch tex data
    const int      tex_w    = raster_state.tex.w;
    const int      tex_h    = raster_state.tex.h;
    const int      tex_bpp  = raster_state.tex.bpp;
    unsigned char *tex_data = raster_state.tex.data;

    float *depth_buffer = rend.depth_buffer;

    // rasterize
    for (int y = AABB[1]; y <= AABB[3]; ++y)
    {
        for (int x = AABB[0]; x <= AABB[2]; ++x)
        {
            vec3 point = {0.5f + (float)x, 0.5f + (float)y, 0.0f};

            float w0 = v3_edgefunc(screen_space[1], screen_space[2], point);
            float w1 = v3_edgefunc(screen_space[2], screen_space[0], point);
            float w2 = v3_edgefunc(screen_space[0], screen_space[1], point);

            if (w0 > 0.0f || w1 > 0.0f || w2 > 0.0f) continue;

            w0 *= inv_area;
            w1 *= inv_area;
            w2 *= inv_area;

            const int index = (y * GRAFIKA_SCREEN_WIDTH) + x;

            const float depth = w0 * screen_space[0][2] + w1 * screen_space[1][2] + w2 * screen_space[2][2];

            float      *old_z = depth_buffer + index;
            const float new_z = depth;

            if (new_z > *old_z) continue;
            *old_z = new_z;

            // Weights, see what i did there ;)
            const float wait0 = w0 * w_vals[0];
            const float wait1 = w1 * w_vals[1];
            const float wait2 = w2 * w_vals[2];

            // correction factor
            const float cf = 1.0f / (wait0 + wait1 + wait2);

            float u = (t.tex[0][0] * wait0 + t.tex[1][0] * wait1 + t.tex[2][0] * wait2) * cf;
            float v = (t.tex[0][1] * wait0 + t.tex[1][1] * wait1 + t.tex[2][1] * wait2) * cf;

            u = SDL_clamp(u, 0.0f, 1.0f);
            v = SDL_clamp(v, 0.0f, 1.0f);

            u *= (float)tex_w - 1;
            v *= (float)tex_h - 1;

            unsigned char *tex_colour = tex_data + (((int)u + tex_w * (int)v) * tex_bpp);

            uint32_t pixel_colour = ((uint32_t)tex_colour[0] << 16) |
                                    ((uint32_t)tex_colour[1] << 8) |
                                    (uint32_t)tex_colour[2];

            grafika_setpixel((uint32_t)x, (uint32_t)y, pixel_colour);
        }
    }
}

void draw_onstart(struct arena *arena)
{
    struct obj obj = raster_state.obj;

    ASSERT(obj.mats, "Object must have atleast a diffuse texture\n");
    raster_state.tex = tex_load(obj.mats[0].map_Kd);

    raster_state.triangles = arena_alloc_aligned(arena, sizeof(struct triangle) * obj.num_f_rows, 16);
}

void draw_object(struct arena *arena)
{
    UNUSED(arena);

    struct triangle    *triangles = raster_state.triangles;
    struct obj          obj       = raster_state.obj;
    float              *pPos      = obj.pos;
    float              *pTex      = obj.texs;
    struct vertindices *pIndices  = obj.indices;

#if 1
#pragma omp parallel
    {
#pragma omp for
        for (size_t i = 0; i < obj.num_f_rows; ++i)
        {
            struct vertindices indices0   = pIndices[i * 3 + 0];
            const int          pos_index0 = indices0.v_idx;

            triangles[i].pos[0][0] = pPos[3 * pos_index0 + 0];
            triangles[i].pos[0][1] = pPos[3 * pos_index0 + 1];
            triangles[i].pos[0][2] = pPos[3 * pos_index0 + 2];

            const int tex_index0   = indices0.vt_idx;
            triangles[i].tex[0][0] = pTex[2 * tex_index0 + 0];
            triangles[i].tex[0][1] = pTex[2 * tex_index0 + 1];

            struct vertindices indices1   = pIndices[i * 3 + 1];
            const int          pos_index1 = indices1.v_idx;

            triangles[i].pos[1][0] = pPos[3 * pos_index1 + 0];
            triangles[i].pos[1][1] = pPos[3 * pos_index1 + 1];
            triangles[i].pos[1][2] = pPos[3 * pos_index1 + 2];

            const int tex_index1   = indices1.vt_idx;
            triangles[i].tex[1][0] = pTex[2 * tex_index1 + 0];
            triangles[i].tex[1][1] = pTex[2 * tex_index1 + 1];

            struct vertindices indices2   = pIndices[i * 3 + 2];
            const int          pos_index2 = indices2.v_idx;

            triangles[i].pos[2][0] = pPos[3 * pos_index2 + 0];
            triangles[i].pos[2][1] = pPos[3 * pos_index2 + 1];
            triangles[i].pos[2][2] = pPos[3 * pos_index2 + 2];

            const int tex_index2   = indices2.vt_idx;
            triangles[i].tex[2][0] = pTex[2 * tex_index2 + 0];
            triangles[i].tex[2][1] = pTex[2 * tex_index2 + 1];
        }

#pragma omp for
        for (size_t i = 0; i < obj.num_f_rows; ++i)
        {
            draw_triangle(triangles[i]);
        }
    }
#else
#pragma omp parallel
    {
#pragma omp for
        for (size_t i = 0; i < obj.num_f_rows; ++i)
        {
            struct triangle t;

            for (size_t j = 0; j < 3; ++j)
            {
                struct vertindices indices = pIndices[i * 3 + j];

                const int posIndex = 3 * indices.v_idx;
                const int texIndex = 2 * indices.vt_idx;

#ifdef CCW_TRIANGLES
                const size_t storeIndex = j; // CCW triangles
#else
                const size_t storeIndex = 2 - j; // CW triangles (which blender uses)
#endif

                t.pos[storeIndex][0] = pPos[posIndex + 0];
                t.pos[storeIndex][1] = pPos[posIndex + 1];
                t.pos[storeIndex][2] = pPos[posIndex + 2];

                t.tex[storeIndex][0] = pTex[texIndex + 0];
                t.tex[storeIndex][1] = pTex[texIndex + 1];
            }
            draw_triangle(t);
        }
    }
#endif
}

void draw_onexit(void)
{
    tex_destroy(&raster_state.tex);
    // raster_state.obj       // cleaned up by arena
    // raster_state.triangles // cleaned up by arena
}