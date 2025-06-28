#define GRAFIKA_TITLE ("grafika - matrix_simd")
#include "common.h"

struct matrix_simd
{
    mat4       model;
    mat4       MVP;
    vec3       cam_pos;
    struct tex tex;
    struct obj obj;
    float     *transformed_vertices;
};

struct matrix_simd raster_state = {0};

// TODO : where does this come from?
#ifdef LH_COORDINATE_SYSTEM
#define VP_MATRIX                                                                                  \
    (mat4)                                                                                         \
    {                                                                                              \
        {0.5f * (float)GRAFIKA_SCREEN_WIDTH, 0.0f, 0.0f, 0.0f},                                    \
            {0.0f, -0.5f * (float)GRAFIKA_SCREEN_HEIGHT, 0.0f, 0.0f},                              \
            {0.0f, 0.0f, 1.0f, 0.0f},                                                              \
            {0.5f * (float)GRAFIKA_SCREEN_WIDTH, 0.5f * (float)GRAFIKA_SCREEN_HEIGHT, 0.0f, 1.0f}, \
    }
#else
#define VP_MATRIX                                                                                  \
    (mat4)                                                                                         \
    {                                                                                              \
        {0.5f * (float)GRAFIKA_SCREEN_WIDTH, 0.0f, 0.0f, 0.0f},                                    \
            {0.0f, 0.5f * (float)GRAFIKA_SCREEN_HEIGHT, 0.0f, 0.0f},                               \
            {0.0f, 0.0f, 1.0f, 0.0f},                                                              \
            {0.5f * (float)GRAFIKA_SCREEN_WIDTH, 0.5f * (float)GRAFIKA_SCREEN_HEIGHT, 0.0f, 1.0f}, \
    }
#endif

static inline __m128 interpolate_values(const __m128 F[3], const __m128 V[3])
{
    return _mm_add_ps(
        _mm_mul_ps(V[2], F[2]),
        _mm_add_ps(
            _mm_mul_ps(V[1], F[1]),
            _mm_mul_ps(V[0], F[0])));
}

static void draw_triangle(vec4 verts[3], vec2 tex_coords[3])
{
    // from the paper, calculate out "A" matrix
    const float a0 = (verts[1][1] * verts[2][3]) - (verts[2][1] * verts[1][3]);
    const float a1 = (verts[2][1] * verts[0][3]) - (verts[0][1] * verts[2][3]);
    const float a2 = (verts[0][1] * verts[1][3]) - (verts[1][1] * verts[0][3]);

    const float b0 = (verts[2][0] * verts[1][3]) - (verts[1][0] * verts[2][3]);
    const float b1 = (verts[0][0] * verts[2][3]) - (verts[2][0] * verts[0][3]);
    const float b2 = (verts[1][0] * verts[0][3]) - (verts[0][0] * verts[1][3]);

    const float c0 = (verts[1][0] * verts[2][1]) - (verts[2][0] * verts[1][1]);
    const float c1 = (verts[2][0] * verts[0][1]) - (verts[0][0] * verts[2][1]);
    const float c2 = (verts[0][0] * verts[1][1]) - (verts[1][0] * verts[0][1]);

    const float detM = (c0 * verts[0][3]) + (c1 * verts[1][3]) + (c2 * verts[2][3]);

    // if (detM >= 0.0f) return; // reject CCW triangles
    if (detM < 0.0f) return; // reject CW triangles

    // set up edge functions (this is A = adj(M))
    vec3 E0 = {a0, b0, c0};
    vec3 E1 = {a1, b1, c1};
    vec3 E2 = {a2, b2, c2};

    // after our cum matrix multiplication, our verts are in some
    // screen space, but we need to divide by w to get correct bounding boxes?
    vec3 bounds[3];
    for (int i = 0; i < 3; ++i)
    {
        bounds[i][0] = verts[i][0] / verts[i][3];
        bounds[i][1] = verts[i][1] / verts[i][3];
    }

    // calculate bounding rectangle
    int AABB[4] = {0};
    AABB_make(bounds, AABB);

    // align to multiples of 4 for better simd'ing
    AABB[0] = AABB[0] - (AABB[0] & (4 - 1)); // p & (a - 1) = p % a
    AABB[2] = AABB[2] + (AABB[2] & (4 - 1));

    __m128 U[3], V[3], Z[3];
    U[0] = _mm_set1_ps(tex_coords[0][0]);
    U[1] = _mm_set1_ps(tex_coords[1][0]);
    U[2] = _mm_set1_ps(tex_coords[2][0]);

    V[0] = _mm_set1_ps(tex_coords[0][1]);
    V[1] = _mm_set1_ps(tex_coords[1][1]);
    V[2] = _mm_set1_ps(tex_coords[2][1]);

    Z[0] = _mm_set1_ps(verts[0][2]);
    Z[1] = _mm_set1_ps(verts[1][2]);
    Z[2] = _mm_set1_ps(verts[2][2]);

    __m128 edge0_A = _mm_set_ps1(E0[0]);
    __m128 edge1_A = _mm_set_ps1(E1[0]);
    __m128 edge2_A = _mm_set_ps1(E2[0]);

    __m128 edge0_B = _mm_set_ps1(E0[1]);
    __m128 edge1_B = _mm_set_ps1(E1[1]);
    __m128 edge2_B = _mm_set_ps1(E2[1]);

    // compute E(x, y) = (x * a) + (y * b) + c, at block origin once
    __m128 E0_origin = _mm_set1_ps((E0[0] * (float)AABB[0] + 0.5f) + (E0[1] * (float)AABB[1] + 0.5f) + E0[2]);
    __m128 E1_origin = _mm_set1_ps((E1[0] * (float)AABB[0] + 0.5f) + (E1[1] * (float)AABB[1] + 0.5f) + E1[2]);
    __m128 E2_origin = _mm_set1_ps((E2[0] * (float)AABB[0] + 0.5f) + (E2[1] * (float)AABB[1] + 0.5f) + E2[2]);

    // generate masks used for tie-breaking rules (not to double-shade along shared edges)
    // We are only checking for true conditions
    // return (E[0] < 0.0f) || (E[0] == 0.0f && E[1] <= 0.0f);
    __m128 edge0_AB_tie_breaker_result = _mm_and_ps(
        _mm_cmplt_ps(edge0_A, _mm_setzero_ps()), /* (E[0] > 0.0f) -> true */
        _mm_or_ps(
            _mm_cmpneq_ps(edge0_A, _mm_setzero_ps()), /* (E[0] == 0.0f && E[1] > 0.0f) -> true */
            _mm_cmpgt_ps(edge0_B, _mm_setzero_ps())));

    __m128 edge1_AB_tie_breaker_result = _mm_and_ps(_mm_cmplt_ps(edge1_A, _mm_setzero_ps()),
                                                    _mm_or_ps(_mm_cmpneq_ps(edge1_A, _mm_setzero_ps()),
                                                              _mm_cmpgt_ps(edge1_B, _mm_setzero_ps())));

    __m128 edge2_AB_tie_breaker_result = _mm_and_ps(_mm_cmplt_ps(edge2_A, _mm_setzero_ps()),
                                                    _mm_or_ps(_mm_cmpneq_ps(edge2_A, _mm_setzero_ps()),
                                                              _mm_cmpgt_ps(edge2_B, _mm_setzero_ps())));

    __m128 starting_x = _mm_set_ps(3.5f, 2.5f, 1.5f, 0.5f);
    __m128 step_y     = _mm_set1_ps(0.5f);

    // pre fetch tex data
    const int      tex_w    = raster_state.tex.w;
    const int      tex_h    = raster_state.tex.h;
    const int      tex_bpp  = raster_state.tex.bpp;
    unsigned char *tex_data = raster_state.tex.data;

    float *depth_buffer = rend.depth_buffer;

    for (int pix_y = AABB[1]; pix_y <= AABB[3];
         pix_y++,
             step_y = _mm_add_ps(step_y, _mm_set1_ps(1.0f)))
    {
        __m128 step_x = starting_x;

        for (int pix_x = AABB[0]; pix_x <= AABB[2];
             pix_x += 4,
                 step_x = _mm_add_ps(step_x, _mm_set1_ps(4.0f)))
        {
            // we are stepping 4 pixels at a time in the x direction
            // a * s
            __m128 edge0_as = _mm_mul_ps(edge0_A, step_x); // TODO : move this out the loop
            __m128 edge1_as = _mm_mul_ps(edge1_A, step_x); // "add" instead
            __m128 edge2_as = _mm_mul_ps(edge2_A, step_x);

            // b * t
            __m128 edge0_bt = _mm_mul_ps(edge0_B, step_y);
            __m128 edge1_bt = _mm_mul_ps(edge1_B, step_y);
            __m128 edge2_bt = _mm_mul_ps(edge2_B, step_y);

            // E(x + s, y + t) = E(x, y) + (a * s) + (t * b)
            __m128 edge_func_0 = _mm_add_ps(E0_origin, _mm_add_ps(edge0_as, edge0_bt));
            __m128 edge_func_1 = _mm_add_ps(E1_origin, _mm_add_ps(edge1_as, edge1_bt));
            __m128 edge_func_2 = _mm_add_ps(E2_origin, _mm_add_ps(edge2_as, edge2_bt));

            // tie Breaker edge testing
#if 1
            // NOTE : since we are using the edge test result as a mask, we want the result to be the
            //          pixels we actually want to draw, not which have "failed" the test, so our simd
            //          test will be the inverse of the single threaded version
            // return (edge_value < 0.0f) || (edge_value == 0.0f && ab_test_result);

            __m128 edge0_positive = _mm_cmpgt_ps(edge_func_0, _mm_setzero_ps());
            __m128 edge1_positive = _mm_cmpgt_ps(edge_func_1, _mm_setzero_ps());
            __m128 edge2_positive = _mm_cmpgt_ps(edge_func_2, _mm_setzero_ps());

            __m128 edge0_zero = _mm_cmpneq_ps(edge_func_0, _mm_setzero_ps());
            __m128 edge1_zero = _mm_cmpneq_ps(edge_func_1, _mm_setzero_ps());
            __m128 edge2_zero = _mm_cmpneq_ps(edge_func_2, _mm_setzero_ps());

            __m128 edge0_func_mask = _mm_and_ps(edge0_positive,
                                                _mm_or_ps(edge0_AB_tie_breaker_result, edge0_zero));

            __m128 edge1_func_mask = _mm_and_ps(edge1_positive,
                                                _mm_or_ps(edge1_AB_tie_breaker_result, edge1_zero));

            __m128 edge2_func_mask = _mm_and_ps(edge2_positive,
                                                _mm_or_ps(edge2_AB_tie_breaker_result, edge2_zero));
#else
            __m128 edge0_func_mask = _mm_cmpgt_ps(edge_func_0, _mm_setzero_ps());
            __m128 edge1_func_mask = _mm_cmpgt_ps(edge_func_1, _mm_setzero_ps());
            __m128 edge2_func_mask = _mm_cmpgt_ps(edge_func_2, _mm_setzero_ps());
#endif
            __m128 edge_test_result = _mm_and_ps(edge0_func_mask, _mm_and_ps(edge1_func_mask, edge2_func_mask));

            int maskInt = _mm_movemask_ps(edge_test_result);
            if (maskInt == 0) continue;

            // calculate interpolation coefficients
            __m128 F[3] = {0};
            F[0]        = edge_func_0;
            F[1]        = edge_func_1;
            F[2]        = edge_func_2;

            // R(x, y) = F0(x, y) + F1(x, y) + F2(x, y)
            // r = 1 / (F0(x, y) + F1(x, y) + F2(x, y))
            __m128 r = _mm_add_ps(F[0], _mm_add_ps(F[1], F[2]));
            r        = _mm_rcp_ps(r);

            F[0] = _mm_mul_ps(F[0], r);
            F[1] = _mm_mul_ps(F[1], r);
            F[2] = _mm_mul_ps(F[2], r);

            const size_t pixel_index = (const size_t)(pix_y * GRAFIKA_SCREEN_WIDTH + pix_x);

            // interpolate depth value
            const __m128 interpolated_depth = interpolate_values(F, Z);

            float       *depth_location = &depth_buffer[pixel_index];
            const __m128 current_depth  = _mm_load_ps(depth_location);

            const __m128 depth_test_result = _mm_cmple_ps(interpolated_depth, current_depth);

            if ((uint16_t)_mm_movemask_ps(depth_test_result) == 0x0) continue;

            // AND depth mask & coverage mask for quads of fragments
            const __m128 write_mask = _mm_and_ps(depth_test_result, edge_test_result);

            // write interpolated Z values
            _mm_maskmoveu_si128(
                _mm_castps_si128(interpolated_depth),
                _mm_castps_si128(write_mask),
                (char *)depth_location);

            // interpolate texture coordinates
            __m128 pixU = interpolate_values(F, U);
            __m128 pixV = interpolate_values(F, V);

            // clamp the vector to the range [0.0f, 1.0f]
            pixU = _mm_max_ps(_mm_min_ps(pixU, _mm_set1_ps(1.0f)), _mm_setzero_ps());
            pixV = _mm_max_ps(_mm_min_ps(pixV, _mm_set1_ps(1.0f)), _mm_setzero_ps());

            pixU = _mm_mul_ps(pixU, _mm_set1_ps((float)(tex_w - 1)));
            pixV = _mm_mul_ps(pixV, _mm_set1_ps((float)(tex_h - 1)));

            // (U + tex.width * V) * tex.bpp
            __m128i tex_offset = _mm_mullo_epi32(
                _mm_set1_epi32(tex_bpp),
                _mm_add_epi32(
                    _mm_cvtps_epi32(pixU),
                    _mm_mullo_epi32(
                        _mm_set1_epi32(tex_w),
                        _mm_cvtps_epi32(pixV))));

#if 1
            uint32_t offset[4] = {0};
            offset[0]          = (uint32_t)_mm_extract_epi32(tex_offset, 0);
            offset[1]          = (uint32_t)_mm_extract_epi32(tex_offset, 1);
            offset[2]          = (uint32_t)_mm_extract_epi32(tex_offset, 2);
            offset[3]          = (uint32_t)_mm_extract_epi32(tex_offset, 3);
#else
            uint32_t *offset = (uint32_t *)&tex_offset; // I think this is unsafe
#endif
            char *sample[4] = {0};
            for (int i = 0; i < 4; ++i)
                sample[i] = (char *)(tex_data + offset[i]);

            __m128i final_colour = _mm_set_epi8(-1, sample[3][0], sample[3][1], sample[3][2],
                                                -1, sample[2][0], sample[2][1], sample[2][2],
                                                -1, sample[1][0], sample[1][1], sample[1][2],
                                                -1, sample[0][0], sample[0][1], sample[0][2]);

            uint32_t *pixel_location = &rend.pixels[pixel_index];

#if 0 // Fabian method
            const __m128i original_pixel_data = _mm_loadu_si128((__m128i *)pixel_location);

            const __m128i write_mask    = _mm_castps_si128(sseWriteMask);
            const __m128i masked_output = _mm_or_si128(_mm_and_si128(write_mask, final_colour),
                                                       _mm_andnot_si128(write_mask, original_pixel_data));

            _mm_storeu_si128((__m128i *)pixel_location, masked_output);
#else
            _mm_maskstore_epi32(
                (int *)pixel_location,
                _mm_castps_si128(write_mask),
                final_colour);

#endif
        }
    }
}

void draw_onstart(struct arena *arena)
{
    UNUSED(arena);

    struct obj obj = raster_state.obj;

    ASSERT(obj.mats, "Object must have atleast a diffuse texture\n");
    raster_state.tex = tex_load(obj.mats[0].map_Kd);

    raster_state.transformed_vertices = arena_alloc_aligned(arena, sizeof(vec4) * obj.num_pos, 16);
}

void draw_object(struct arena *arena)
{
    UNUSED(arena);

    mat4 cum_matrix = {0};
    m4_mul_m4(VP_MATRIX, raster_state.MVP, cum_matrix);

    float *trans_verts = raster_state.transformed_vertices;

    const struct obj obj = raster_state.obj;
    float           *pos = obj.pos;

    // #pragma omp parallel
    {
        //  transform all verts to screen space
#pragma omp parallel for
        for (size_t i = 0; i < obj.num_pos; i++)
        {
            vec4 vertex = {pos[i * 3 + 0],
                           pos[i * 3 + 1],
                           pos[i * 3 + 2],
                           1.0f};
            m4_mul_v4(cum_matrix, vertex, &trans_verts[i * 4]);
        }

        float              *obj_pos     = trans_verts;
        float              *obj_tex     = obj.texs;
        struct vertindices *obj_indices = obj.indices;

#pragma omp parallel for
        for (size_t i = 0; i < obj.num_f_rows; ++i)
        {
            vec4 vertices[3]   = {0};
            vec2 tex_coords[3] = {0};

            for (size_t j = 0; j < 3; ++j)
            {
                struct vertindices indices = obj_indices[i * 3 + j];

                const int pos_index = indices.v_idx;
                const int tex_index = indices.vt_idx;

#ifdef CCW_TRIANGLES
                const size_t storeIndex = j; // CCW triangles
#else
                const size_t storeIndex = 2 - j; // CW triangles (which blender uses)
#endif

                float *p                = obj_pos + 4 * pos_index;
                vertices[storeIndex][0] = p[0];
                vertices[storeIndex][1] = p[1];
                vertices[storeIndex][2] = p[2];
                vertices[storeIndex][3] = p[3];

                tex_coords[storeIndex][0] = obj_tex[2 * tex_index + 0];
                tex_coords[storeIndex][1] = obj_tex[2 * tex_index + 1];
            }
            draw_triangle(vertices, tex_coords);
        }
    }
}

void draw_onexit(void)
{
    tex_destroy(&raster_state.tex);
    // raster_state.obj                     // cleaned up by arena
    // raster_state.transformed_vertices    // cleaned up by arena
}