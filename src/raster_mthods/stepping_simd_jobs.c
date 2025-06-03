#include "common.h"

#define RABOTY_IMPLEMENTATION
#include "../raboty.h"

struct triangle
{
    vec3 pos[3];
    vec2 tex[3];
};

struct stepping_simd_jobs
{
    mat4       model;
    mat4       MVP;
    vec3       cam_pos;
    tex_t      tex;
    struct obj obj;
};

struct stepping_simd_jobs raster_state = {0};

static void draw_triangle(const struct triangle t, int region[4])
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

    // back face culling (surface normal)
    vec3 sub10, sub20, normal;
    v3_sub(ndc[1], ndc[0], sub10);
    v3_sub(ndc[2], ndc[0], sub20);
    v3_cross(sub10, sub20, normal);
    if (normal[2] > 0.0f)
        return;

    vec3 screen_space[3];
    for (int i = 0; i < 3; ++i)
    {
        screen_space[i][0] = (ndc[i][0] + 1.0f) * 0.5f * GRAFIKA_SCREEN_WIDTH;
        screen_space[i][1] = (ndc[i][1] + 1.0f) * 0.5f * GRAFIKA_SCREEN_HEIGHT;
        screen_space[i][2] = (ndc[i][2] + 1.0f) * 0.5f;
    }

    // compute the triangles bounding box
    int min_x = INT_MAX, min_y = INT_MAX;
    int max_x = INT_MIN, max_y = INT_MIN;
    for (int i = 0; i < 3; ++i)
    {
        if (screen_space[i][0] < min_x) min_x = (int)screen_space[i][0];
        if (screen_space[i][1] < min_y) min_y = (int)screen_space[i][1];
        if (screen_space[i][0] > max_x) max_x = (int)screen_space[i][0];
        if (screen_space[i][1] > max_y) max_y = (int)screen_space[i][1];
    }

    // clamp to region
    int AABB[4] = {0};
    AABB[0]     = (min_x < region[0]) ? region[0] : min_x;
    AABB[1]     = (min_y < region[1]) ? region[1] : min_y;
    AABB[2]     = (max_x > region[2]) ? region[2] : max_x;
    AABB[3]     = (max_y > region[3]) ? region[3] : max_y;

    AABB[0] = AABB[0] & ~3;       // round down AABB[0] to the nearest multiple of 4
    AABB[2] = (AABB[2] + 3) & ~3; // round up AABB[2] to the next multiple of 4 - 1

    const float dY0 = screen_space[2][1] - screen_space[1][1], dX0 = screen_space[1][0] - screen_space[2][0];
    const float dY1 = screen_space[0][1] - screen_space[2][1], dX1 = screen_space[2][0] - screen_space[0][0];
    const float dY2 = screen_space[1][1] - screen_space[0][1], dX2 = screen_space[0][0] - screen_space[1][0];

    const float C0 = (screen_space[2][0] * screen_space[1][1]) - (screen_space[2][1] * screen_space[1][0]);
    const float C1 = (screen_space[0][0] * screen_space[2][1]) - (screen_space[0][1] * screen_space[2][0]);
    const float C2 = (screen_space[1][0] * screen_space[0][1]) - (screen_space[1][1] * screen_space[0][0]);

    // step vectors (for 4-pixel stride)
    const __m128 A0 = _mm_set1_ps(dY0);
    const __m128 A1 = _mm_set1_ps(dY1);
    const __m128 A2 = _mm_set1_ps(dY2);

    const __m128 B0 = _mm_set1_ps(dX0);
    const __m128 B1 = _mm_set1_ps(dX1);
    const __m128 B2 = _mm_set1_ps(dX2);

    __m128 x_step_amount    = _mm_setr_ps(0.5f, 1.5f, 2.5f, 3.5f);
    __m128 x_starting_pixel = _mm_set1_ps((float)AABB[0] + 0.5f);
    __m128 mm_x             = _mm_add_ps(x_step_amount, x_starting_pixel);
    __m128 A0_start         = _mm_mul_ps(A0, mm_x);
    __m128 A1_start         = _mm_mul_ps(A1, mm_x);
    __m128 A2_start         = _mm_mul_ps(A2, mm_x);

    __m128 y_starting_pixel = _mm_set1_ps((float)AABB[1] + 0.5f);
    __m128 B0_start         = _mm_mul_ps(B0, y_starting_pixel);
    __m128 B1_start         = _mm_mul_ps(B1, y_starting_pixel);
    __m128 B2_start         = _mm_mul_ps(B2, y_starting_pixel);

    __m128 E0 = _mm_add_ps(_mm_add_ps(A0_start, B0_start), _mm_set1_ps(C0));
    __m128 E1 = _mm_add_ps(_mm_add_ps(A1_start, B1_start), _mm_set1_ps(C1));
    __m128 E2 = _mm_add_ps(_mm_add_ps(A2_start, B2_start), _mm_set1_ps(C2));

    __m128 B0_inc = B0;
    __m128 B1_inc = B1;
    __m128 B2_inc = B2;

    __m128 mm_4   = _mm_set1_ps(4.0f);
    __m128 A0_inc = _mm_mul_ps(A0, mm_4);
    __m128 A1_inc = _mm_mul_ps(A1, mm_4);
    __m128 A2_inc = _mm_mul_ps(A2, mm_4);

    const __m128 mm_w_vals0 = _mm_set1_ps(w_vals[0]);
    const __m128 mm_w_vals1 = _mm_set1_ps(w_vals[1]);
    const __m128 mm_w_vals2 = _mm_set1_ps(w_vals[2]);

    const __m128 U0 = _mm_set1_ps(t.tex[0][0]);
    const __m128 U1 = _mm_set1_ps(t.tex[1][0]);
    const __m128 U2 = _mm_set1_ps(t.tex[2][0]);

    const __m128 V0 = _mm_set1_ps(t.tex[0][1]);
    const __m128 V1 = _mm_set1_ps(t.tex[1][1]);
    const __m128 V2 = _mm_set1_ps(t.tex[2][1]);

    const float inv_area = 1.0f / (dX1 * dY2 - dY1 * dX2);
    screen_space[1][2]   = (screen_space[1][2] - screen_space[0][2]) * inv_area;
    screen_space[2][2]   = (screen_space[2][2] - screen_space[0][2]) * inv_area;

    const float depth_step = dY1 * screen_space[1][2] + dY2 * screen_space[2][2];

    // pre fetch tex data
    const int      tex_w    = raster_state.tex.w;
    const int      tex_h    = raster_state.tex.h;
    const int      tex_bpp  = raster_state.tex.bpp;
    unsigned char *tex_data = raster_state.tex.data;

    float *depth_buffer = rend.depth_buffer;

    __m128 zstep = _mm_set1_ps(depth_step * 4);
    __m128 ss0   = _mm_set1_ps(screen_space[0][2]);
    __m128 ss1   = _mm_set1_ps(screen_space[1][2]);
    __m128 ss2   = _mm_set1_ps(screen_space[2][2]);

    __m128 mm_1 = _mm_set1_ps(1.0f);
    __m128 mm_0 = _mm_setzero_ps();
    //__m128 sign_mask = _mm_set1_ps(-0.0f); // sets only the sign bit (0x80000000)

    __m128 mm_tex_w = _mm_set1_ps((float)tex_w - 1);
    __m128 mm_tex_h = _mm_set1_ps((float)tex_h - 1);

    __m128i mmi_tex_w   = _mm_set1_epi32(tex_w);
    __m128i mmi_tex_bpp = _mm_set1_epi32(tex_bpp);

    const __m128i shuffle_mask = _mm_setr_epi8(2, 1, 0, -1,   // texel0 (B, G, R, A)
                                               6, 5, 4, -1,   // texel1
                                               10, 9, 8, -1,  // texel2
                                               14, 13, 12, -1 // texel3
    );

    // TIMED_BLOCK_BEGIN(raster_pixels);
    for (int y  = AABB[1]; y <= AABB[3]; ++y,
             E0 = _mm_add_ps(E0, B0_inc),
             E1 = _mm_add_ps(E1, B1_inc),
             E2 = _mm_add_ps(E2, B2_inc))
    {
        __m128 alpha = E0;
        __m128 betaa = E1;
        __m128 gamma = E2;

        __m128 z = _mm_add_ps(_mm_add_ps(ss0, _mm_mul_ps(ss1, betaa)), _mm_mul_ps(ss2, gamma));

        for (int x     = AABB[0]; x <= AABB[2]; x += 4,
                 alpha = _mm_add_ps(alpha, A0_inc),
                 betaa = _mm_add_ps(betaa, A1_inc),
                 gamma = _mm_add_ps(gamma, A2_inc),
                 z     = _mm_add_ps(z, zstep))
        {
            __m128 EdgeFuncTestResult = _mm_and_ps(_mm_and_ps(
                                                       _mm_cmpgt_ps(alpha, mm_0),
                                                       _mm_cmpgt_ps(betaa, mm_0)),
                                                   _mm_cmpgt_ps(gamma, mm_0));

            if (_mm_movemask_ps(EdgeFuncTestResult) == 0) continue;

            __m128 wait0 = _mm_mul_ps(alpha, mm_w_vals0); // Bary A
            __m128 wait1 = _mm_mul_ps(betaa, mm_w_vals1); // B
            __m128 wait2 = _mm_mul_ps(gamma, mm_w_vals2); // C

            const __m128 cf = _mm_rcp_ps(_mm_add_ps(_mm_add_ps(wait0, wait1), wait2));

            wait0 = _mm_mul_ps(wait0, cf);
            wait1 = _mm_mul_ps(wait1, cf);
            wait2 = _mm_mul_ps(wait2, cf);

            const size_t pixel_index = (const size_t)(y * GRAFIKA_SCREEN_WIDTH + x);

            float       *depth_location = &depth_buffer[pixel_index];
            const __m128 current_depth  = _mm_load_ps(&rend.depth_buffer[pixel_index]);

            const __m128 depth_test = _mm_cmple_ps(z, current_depth);

            if (_mm_movemask_ps(depth_test) == 0) continue;

            const __m128 sseWriteMask = _mm_and_ps(depth_test, EdgeFuncTestResult);
            __m128       blended      = _mm_blendv_ps(current_depth, z, sseWriteMask);

            _mm_store_ps(depth_location, blended);

            __m128 pixU = _mm_add_ps(_mm_mul_ps(U2, wait2),
                                     _mm_add_ps(
                                         _mm_mul_ps(U1, wait1),
                                         _mm_mul_ps(U0, wait0)));
            //__m128 pixU = _mm_fmadd_ps(U2, wait2, _mm_fmadd_ps(U1, wait1, _mm_mul_ps(U0, wait0)));
            pixU = _mm_max_ps(_mm_min_ps(pixU, mm_1), mm_0);
            pixU = _mm_mul_ps(pixU, mm_tex_w);

            __m128 pixV = _mm_add_ps(_mm_mul_ps(V2, wait2),
                                     _mm_add_ps(
                                         _mm_mul_ps(V1, wait1),
                                         _mm_mul_ps(V0, wait0)));
            //__m128 pixV = _mm_fmadd_ps(V2, wait2, _mm_fmadd_ps(V1, wait1, _mm_mul_ps(V0, wait0)));
            pixV = _mm_max_ps(_mm_min_ps(pixV, mm_1), mm_0);
            pixV = _mm_mul_ps(pixV, mm_tex_h);

            // (U + tex.width * V) * tex.bpp
            __m128i tex_offset = _mm_mullo_epi32(_mm_add_epi32(
                                                     _mm_cvtps_epi32(pixU),
                                                     _mm_mullo_epi32(_mm_cvtps_epi32(pixV), mmi_tex_w)),
                                                 mmi_tex_bpp);

#if 0
            // const int32_t *offsets = (const int32_t *)&texOffset;
            // for (size_t i = 0; i < 4; i++)
            //{
            //     if (sseWriteMask[i])
            //     {
            //         unsigned char *texcolour0 = texdata + offsets[i];
            //         unsigned char *dst        = (unsigned char *)&rend.pixels[pixel_index + i];
            //         dst[0]                    = texcolour0[2];
            //         dst[1]                    = texcolour0[1];
            //         dst[2]                    = texcolour0[0];
            //     }
            // }

            const int32_t *offsets = (const int32_t *)&texOffset;

            __m128i texel0 = _mm_load_si128((__m128i *)(texdata + offsets[0]));
            __m128i texel1 = _mm_load_si128((__m128i *)(texdata + offsets[1]));
            __m128i texel2 = _mm_load_si128((__m128i *)(texdata + offsets[2]));
            __m128i texel3 = _mm_load_si128((__m128i *)(texdata + offsets[3]));

            //__m128i final_colour = _mm_setr_epi32((int)texel0[0],
            //                                      (int)texel1[0],
            //                                      (int)texel2[0],
            //                                      (int)texel3[0]);

            __m128i final_colour = _mm_shuffle_epi8(
                _mm_setr_epi32((int)texel0[0],
                               (int)texel1[0],
                               (int)texel2[0],
                               (int)texel3[0]),
                shuffle_mask);

                // unsigned char *texcolour0 = texdata + _mm_extract_epi32(texOffset, 0);
                // unsigned char *texcolour1 = texdata + _mm_extract_epi32(texOffset, 1);
                // unsigned char *texcolour2 = texdata + _mm_extract_epi32(texOffset, 2);
                // unsigned char *texcolour3 = texdata + _mm_extract_epi32(texOffset, 3);

                //__m128i final_colour = _mm_setr_epi8((char)texcolour0[2], (char)texcolour0[1], (char)texcolour0[0], -1,
                //                                     (char)texcolour1[2], (char)texcolour1[1], (char)texcolour1[0], -1,
                //                                     (char)texcolour2[2], (char)texcolour2[1], (char)texcolour2[0], -1,
                //                                     (char)texcolour3[2], (char)texcolour3[1], (char)texcolour3[0], -1);

                // uint32_t *pixel_location = &rend.pixels[pixel_index];
#else
            // load 4 offsets from texOffset (assuming texOffset is __m128i)
            const int32_t *offsets = (const int32_t *)&tex_offset;

            __m128i texels = _mm_setr_epi32(*(const int32_t *)(tex_data + offsets[0]),
                                            *(const int32_t *)(tex_data + offsets[1]),
                                            *(const int32_t *)(tex_data + offsets[2]),
                                            *(const int32_t *)(tex_data + offsets[3]));

            __m128i final_colour = _mm_shuffle_epi8(texels, shuffle_mask);
#endif

#if 0 /* Fabian method */
            uint32_t *pixel_location = &rend.pixels[pixel_index];
            // const __m128i original_pixel_data = _mm_load_si128((__m128i *)pixel_location);

            const __m128i write_mask    = _mm_castps_si128(sseWriteMask);
            const __m128i masked_output = _mm_or_si128(_mm_and_si128(write_mask, final_colour),
                                                       _mm_andnot_si128(write_mask, *(__m128i *)pixel_location));

            _mm_store_si128((__m128i *)pixel_location, masked_output);
#elif 1
            uint32_t *pixel_location = &rend.pixels[pixel_index];
            // const __m128i original_pixel_data = _mm_load_si128((__m128i *)pixel_location);
            const __m128i write_mask    = _mm_castps_si128(sseWriteMask);
            const __m128i masked_output = _mm_blendv_epi8(*(__m128i *)pixel_location, final_colour, write_mask);
            _mm_store_si128((__m128i *)pixel_location, masked_output);
#else
            const __m128i write_mask = _mm_castps_si128(sseWriteMask);
            _mm_maskstore_epi32((int *)&rend.pixels[pixel_index], write_mask, final_colour);
#endif
        }
    }
    // TIMED_BLOCK_END(raster_pixels);
}

struct raster_job_data
{
    struct triangle *triangle_list;
    uint32_t         triangle_count;
    int              region[4];
};

static void raster_job(void *data)
{
    struct raster_job_data *d = (struct raster_job_data *)data;
    ASSERT(d, "Raster job data is inavlid\n");

    for (uint32_t i = 0; i < d->triangle_count; i++)
    {
        draw_triangle(d->triangle_list[i], d->region);
    }
}

void draw_onstart(struct arena *arena)
{
    UNUSED(arena);

    jobs_init();

    struct obj obj = raster_state.obj;

    ASSERT(obj.mats, "Object must have atleast a diffuse texture\n");
    raster_state.tex = tex_load(obj.mats[0].map_Kd);
}

void draw_object(struct arena *arena)
{
    UNUSED(arena);

    struct arena tmp_arena = arena_create_tmp(arena);

    const struct obj obj = raster_state.obj;

    float              *obj_pos     = obj.pos;
    float              *obj_tex     = obj.texs;
    struct vertindices *obj_indices = obj.indices;

    for (size_t i = 0; i < obj.num_f_rows; ++i)
    {
        struct triangle *t = ma_push_struct(&tmp_arena, struct triangle);

        for (size_t j = 0; j < 3; ++j)
        {
            struct vertindices indices = obj_indices[i * 3 + j];

            const int pos_index = indices.v_idx;
            const int tex_index = indices.vt_idx;
#if 0
                const size_t store_index = j; // CCW triangles
#else
            const size_t store_index = 2 - j; // CW triangles (which blender uses)
#endif
            t->pos[store_index][0] = obj_pos[3 * pos_index + 0];
            t->pos[store_index][1] = obj_pos[3 * pos_index + 1];
            t->pos[store_index][2] = obj_pos[3 * pos_index + 2];

            t->tex[store_index][0] = obj_tex[2 * tex_index + 0];
            t->tex[store_index][1] = obj_tex[2 * tex_index + 1];
        }
    }

    // 1 - have a list of all triangles to be rendered
    // 2 - pass all triangles to all threads
    // 3 - give each thread a "region" to work on

    int tile_count_x = 2;
    int tile_count_y = 2;

    int tile_w = GRAFIKA_SCREEN_WIDTH / tile_count_x;
    int tile_h = GRAFIKA_SCREEN_HEIGHT / tile_count_y;

    tile_w = ((tile_w + 3) / 4) * 4; // rouding up to nearest 4

    for (int32_t tile_y = 0; tile_y < tile_count_y; tile_y++)
    {
        for (int32_t tile_x = 0; tile_x < tile_count_x; tile_x++)
        {
            int render_region[4];
            render_region[0] = (tile_x * tile_w);
            render_region[2] = (render_region[0] + tile_w);

            if (tile_x == (tile_count_x - 1))
            {
                render_region[2] = GRAFIKA_SCREEN_WIDTH; // clip to max screen width
            }

            render_region[1] = (tile_y * tile_h);
            render_region[3] = (render_region[1] + tile_h);

            if (tile_y == (tile_count_y - 1))
            {
                render_region[3] = GRAFIKA_SCREEN_HEIGHT; // clip to max screen height
            }

            struct raster_job_data *rd = ma_push_struct(&tmp_arena, struct raster_job_data);

            // LOG("Region : %d, %d, %d, %d\n", render_region[0], render_region[1],
            //     render_region[2], render_region[3]);

            rd->region[0] = render_region[0];
            rd->region[1] = render_region[1];
            rd->region[2] = render_region[2];
            rd->region[3] = render_region[3];

            rd->triangle_count = obj.num_f_rows;

            rd->triangle_list = (struct triangle *)tmp_arena.base;

            job_submit((struct job){.arguments = rd, .function = raster_job});

            // raster_job(rd);
        }
    }

    jobs_complete_all_work();
}

void draw_onexit(void)
{
    tex_destroy(&raster_state.tex);
    // raster_state.obj       // cleaned up by arena

    jobs_shutdown();
}