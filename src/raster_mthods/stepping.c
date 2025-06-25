#include "common.h"

struct triangle
{
    vec3 pos[3];
    vec2 tex[3];
};

struct stepping
{
    mat4       model;
    mat4       MVP;
    vec3       cam_pos;
    tex_t      tex;
    struct obj obj;
};

struct stepping raster_state = {0};

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
    vec3 ndc[3], w_vals;
    for (size_t i = 0; i < 3; i++)
    {
        w_vals[i] = 1.0f / clip_space[i][3]; // 1.0f / w
        ndc[i][0] = clip_space[i][0] * w_vals[i];
        ndc[i][1] = clip_space[i][1] * w_vals[i];
        ndc[i][2] = clip_space[i][2] * w_vals[i];
    }

    vec3 screen_space[3];
    for (int i = 0; i < 3; ++i)
    {
        screen_space[i][0] = (ndc[i][0] + 1.0f) * 0.5f * GRAFIKA_SCREEN_WIDTH;
        screen_space[i][1] = (1.0f - ndc[i][1]) * 0.5f * GRAFIKA_SCREEN_HEIGHT;
        screen_space[i][2] = (ndc[i][2] + 1.0f) * 0.5f;
    }

    const float dY0 = screen_space[2][1] - screen_space[1][1], dX0 = screen_space[1][0] - screen_space[2][0];
    const float dY1 = screen_space[0][1] - screen_space[2][1], dX1 = screen_space[2][0] - screen_space[0][0];
    const float dY2 = screen_space[1][1] - screen_space[0][1], dX2 = screen_space[0][0] - screen_space[1][0];

    // back face culling
    float signed_area = (dX1 * dY2 - dY1 * dX2);

    // If cross_z > 0, the triangle is counter-clockwise (CCW).
    // If cross_z < 0, the triangle is clockwise (CW).
    // If cross_z == 0, the triangle is degenerate (its points are collinear).
    if (signed_area > 0.0f) return;

    const float C0 = (screen_space[2][0] * screen_space[1][1]) - (screen_space[2][1] * screen_space[1][0]);
    const float C1 = (screen_space[0][0] * screen_space[2][1]) - (screen_space[0][1] * screen_space[2][0]);
    const float C2 = (screen_space[1][0] * screen_space[0][1]) - (screen_space[1][1] * screen_space[0][0]);

    // calculate bounding rectangle
    int AABB[4] = {0};
    AABB_make(screen_space, AABB);

    // const vec3 P     = {minX + 0.5f, minY + 0.5f, 0.0f};
    float alpha = (dY0 * (float)AABB[0]) + (dX0 * (float)AABB[1]) + C0;
    float betaa = (dY1 * (float)AABB[0]) + (dX1 * (float)AABB[1]) + C1;
    float gamma = (dY2 * (float)AABB[0]) + (dX2 * (float)AABB[1]) + C2;

    const float inv_area = 1.0f / signed_area;

    screen_space[1][2] = (screen_space[1][2] - screen_space[0][2]) * inv_area;
    screen_space[2][2] = (screen_space[2][2] - screen_space[0][2]) * inv_area;

    const float zstep = dY1 * screen_space[1][2] + dY2 * screen_space[2][2];

    // pre fetch tex data
    const int      tex_w    = raster_state.tex.w;
    const int      tex_h    = raster_state.tex.h;
    const int      tex_bpp  = raster_state.tex.bpp;
    unsigned char *tex_data = raster_state.tex.data;

    float *depth_buffer = rend.depth_buffer;

    for (int y = AABB[1]; y <= AABB[3]; ++y)
    {
        // barycentric coordinates at start of row
        float w0 = alpha;
        float w1 = betaa;
        float w2 = gamma;

        float depth = screen_space[0][2] + (screen_space[1][2] * betaa) + (screen_space[2][2] * gamma);

        for (int x = AABB[0]; x <= AABB[2]; ++x,            //
                                            w0 += dY0,      //
                                            w1 += dY1,      //
                                            w2 += dY2,      // one step to the right
                                            depth += zstep) // step the depth
        {
            // if (w0 < 0.0f || w1 < 0.0f || w2 < 0.0f) continue;
            if (w0 > 0.0f || w1 > 0.0f || w2 > 0.0f) continue;

            const int index = (y * GRAFIKA_SCREEN_WIDTH) + x;

            float      *oldZ = depth_buffer + index;
            const float invZ = depth;

            // depth test
            if (invZ > *oldZ) continue;
            *oldZ = invZ;

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
        // step one row
        alpha += dX0;
        betaa += dX1;
        gamma += dX2;
    }
}

void draw_onstart(struct arena *arena)
{
    UNUSED(arena);

    struct obj obj = raster_state.obj;

    ASSERT(obj.mats, "Object must have atleast a diffuse texture\n");
    raster_state.tex = tex_load(obj.mats[0].map_Kd);
}

void draw_object(struct arena *arena)
{
    UNUSED(arena);

#pragma omp parallel
    {
        const struct obj obj = raster_state.obj;

        float              *obj_pos     = obj.pos;
        float              *obj_tex     = obj.texs;
        struct vertindices *obj_indices = obj.indices;

#pragma omp for
        for (size_t i = 0; i < obj.num_f_rows; ++i)
        {
            struct triangle t = {0};

            for (size_t j = 0; j < 3; ++j)
            {
                struct vertindices indices = obj_indices[i * 3 + j];

                const int pos_index = indices.v_idx;
                const int tex_index = indices.vt_idx;

#ifdef CCW_TRIANGLES
                const size_t store_index = j; // CCW triangles
#else
                const size_t store_index = 2 - j; // CW triangles (which blender uses)
#endif

                t.pos[store_index][0] = obj_pos[3 * pos_index + 0];
                t.pos[store_index][1] = obj_pos[3 * pos_index + 1];
                t.pos[store_index][2] = obj_pos[3 * pos_index + 2];

                t.tex[store_index][0] = obj_tex[2 * tex_index + 0];
                t.tex[store_index][1] = obj_tex[2 * tex_index + 1];
            }
            draw_triangle(t);
        }
    }
}

void draw_onexit(void)
{
    tex_destroy(&raster_state.tex);
    // raster_state.obj       // cleaned up by arena
}