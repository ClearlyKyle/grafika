#include "common.h"

struct rasterstate raster_state         = {0};
static float      *transformed_vertices = NULL;

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

    // the sign of the determinant gives the orientation of the triangle: a counterclockwise order gives a positive determinant
    if (detM >= 0.0f)
        return;

    // set up edge functions (this is A = adj(M))
    vec3 E0 = {a0, b0, c0};
    vec3 E1 = {a1, b1, c1};
    vec3 E2 = {a2, b2, c2};

    // projection : projected = vert / vert.w
    vec3 proj[3] = {0};
    for (int i = 0; i < 3; i++)
    {
        proj[i][0] = verts[i][0] / verts[i][3];
        proj[i][1] = verts[i][1] / verts[i][3];
        proj[i][2] = verts[i][2] / verts[i][3];
    }

    // get the bounding box of the triangle
    float bb_minX = proj[0][0];
    float bb_minY = proj[0][1];
    float bb_maxX = proj[0][0];
    float bb_maxY = proj[0][1];

    for (int i = 1; i < 3; ++i)
    {
        if (proj[i][0] < bb_minX) bb_minX = proj[i][0];
        if (proj[i][1] < bb_minY) bb_minY = proj[i][1];
        if (proj[i][0] > bb_maxX) bb_maxX = proj[i][0];
        if (proj[i][1] > bb_maxY) bb_maxY = proj[i][1];
    }

    const int screen_width_minus_one  = GRAFIKA_SCREEN_WIDTH - 1;
    const int screen_height_minus_one = GRAFIKA_SCREEN_HEIGHT - 1;

    // clamp values to screen space
    int minX = max(0, min((int)bb_minX, screen_width_minus_one));
    int minY = max(0, min((int)bb_minY, screen_height_minus_one));
    int maxX = max(0, min((int)bb_maxX, screen_width_minus_one));
    int maxY = max(0, min((int)bb_maxY, screen_height_minus_one));

    // align to multiples of 4 for better simd'ing
    minX = minX - (minX & (4 - 1)); // p & (a - 1) = p % a
    maxX = maxX + (minX & (4 - 1));

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
    __m128 E0_origin = _mm_set1_ps((E0[0] * (float)minX) + (E0[1] * (float)minY) + E0[2]);
    __m128 E1_origin = _mm_set1_ps((E1[0] * (float)minX) + (E1[1] * (float)minY) + E1[2]);
    __m128 E2_origin = _mm_set1_ps((E2[0] * (float)minX) + (E2[1] * (float)minY) + E2[2]);

    // generate masks used for tie-breaking rules (not to double-shade along shared edges)
    // We are only checking for true conditions
    __m128 edge0_AB_tie_breaker_result = _mm_or_ps(
        _mm_cmpgt_ps(edge0_A, _mm_setzero_ps()), /* (E[0] > 0.0f) -> true */
        _mm_and_ps(
            _mm_cmpeq_ps(edge0_A, _mm_setzero_ps()), /* (E[0] == 0.0f && E[1] > 0.0f) -> true */
            _mm_cmpge_ps(edge0_B, _mm_setzero_ps())));

    __m128 edge1_AB_tie_breaker_result = _mm_or_ps(_mm_cmpgt_ps(edge1_A, _mm_setzero_ps()),
                                                   _mm_and_ps(_mm_cmpge_ps(edge1_B, _mm_setzero_ps()),
                                                              _mm_cmpeq_ps(edge1_A, _mm_setzero_ps())));

    __m128 edge2_AB_tie_breaker_result = _mm_or_ps(_mm_cmpgt_ps(edge2_A, _mm_setzero_ps()),
                                                   _mm_and_ps(_mm_cmpge_ps(edge2_B, _mm_setzero_ps()),
                                                              _mm_cmpeq_ps(edge2_A, _mm_setzero_ps())));

    __m128 starting_x = _mm_set_ps(3.5f, 2.5f, 1.5f, 0.5f);
    __m128 step_y     = _mm_set1_ps(0.5f);

    // pre fetch tex data
    const int      tex_w    = raster_state.tex.w;
    const int      tex_h    = raster_state.tex.h;
    const int      tex_bpp  = raster_state.tex.bpp;
    unsigned char *tex_data = raster_state.tex.data;

    float *depth_buffer = rend.depth_buffer;

    for (int pix_y = minY; pix_y <= maxY;
         pix_y++,
             step_y = _mm_add_ps(step_y, _mm_set1_ps(1.0f)))
    {
        __m128 step_x = starting_x;

        for (int pix_x = minX; pix_x <= maxX;
             pix_x += 4,
                 step_x = _mm_add_ps(step_x, _mm_set1_ps(4.0f)))
        {
            // we are stepping 4 pixels at a time in the x direction
            // a * s
            __m128 edge0_as = _mm_mul_ps(edge0_A, step_x);
            __m128 edge1_as = _mm_mul_ps(edge1_A, step_x);
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
            __m128 edge0_positive  = _mm_cmplt_ps(edge_func_0, _mm_setzero_ps());
            __m128 edge0_negative  = _mm_cmpgt_ps(edge_func_0, _mm_setzero_ps());
            __m128 edge0_func_mask = _mm_or_ps(edge0_positive,
                                               _mm_andnot_ps(edge0_negative,
                                                             edge0_AB_tie_breaker_result));
            __m128 edge1_positive  = _mm_cmplt_ps(edge_func_1, _mm_setzero_ps());
            __m128 edge1_negative  = _mm_cmpgt_ps(edge_func_1, _mm_setzero_ps());
            __m128 edge1_func_mask = _mm_or_ps(edge1_positive,
                                               _mm_andnot_ps(edge1_negative,
                                                             edge1_AB_tie_breaker_result));

            __m128 edge2_positive  = _mm_cmplt_ps(edge_func_2, _mm_setzero_ps());
            __m128 edge2_negative  = _mm_cmpgt_ps(edge_func_2, _mm_setzero_ps());
            __m128 edge2_func_mask = _mm_or_ps(edge2_positive,
                                               _mm_andnot_ps(edge2_negative,
                                                             edge2_AB_tie_breaker_result));
#else
            __m128 edge0_func_mask = _mm_cmplt_ps(edge_func_0, _mm_setzero_ps());
            __m128 edge1_func_mask = _mm_cmplt_ps(edge_func_1, _mm_setzero_ps());
            __m128 edge2_func_mask = _mm_cmplt_ps(edge_func_2, _mm_setzero_ps());
#endif
            // Combine resulting masks of all three edges
            __m128 edgeFuncTestResult = _mm_and_ps(edge0_func_mask, _mm_and_ps(edge1_func_mask, edge2_func_mask));

            uint16_t maskInt = (uint16_t)_mm_movemask_ps(edgeFuncTestResult);

            if (maskInt == 0x0) continue;

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
            const __m128 sseZInterpolated = interpolate_values(F, Z);

            float       *pDepth          = &depth_buffer[pixel_index];
            const __m128 sseDepthCurrent = _mm_loadu_ps(pDepth);

            const __m128 sseDepthRes = _mm_cmple_ps(sseZInterpolated, sseDepthCurrent);

            if ((uint16_t)_mm_movemask_ps(sseDepthRes) == 0x0) continue;

            // AND depth mask & coverage mask for quads of fragments
            const __m128 sseWriteMask = _mm_and_ps(sseDepthRes, edgeFuncTestResult);

            // write interpolated Z values
            _mm_maskmoveu_si128(
                _mm_castps_si128(sseZInterpolated),
                _mm_castps_si128(sseWriteMask),
                (char *)pDepth);

            // interpolate texture coordinates
            __m128 pixU = interpolate_values(F, U);
            __m128 pixV = interpolate_values(F, V);

            // clamp the vector to the range [0.0f, 1.0f]
            pixU = _mm_max_ps(_mm_min_ps(pixU, _mm_set1_ps(1.0f)), _mm_setzero_ps());
            pixV = _mm_max_ps(_mm_min_ps(pixV, _mm_set1_ps(1.0f)), _mm_setzero_ps());

            pixU = _mm_mul_ps(pixU, _mm_set1_ps((float)tex_w));
            pixV = _mm_mul_ps(pixV, _mm_set1_ps((float)tex_h));

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

            __m128i final_colour = _mm_setr_epi8(sample[0][2], sample[0][0], sample[0][1], -1,
                                                 sample[1][2], sample[1][0], sample[1][1], -1,
                                                 sample[2][2], sample[2][0], sample[2][1], -1,
                                                 sample[3][2], sample[3][0], sample[3][1], -1);

            uint32_t *pixel_location = &rend.pixels[pixel_index];

#if 1 // Fabian method
            const __m128i original_pixel_data = _mm_loadu_si128((__m128i *)pixel_location);

            const __m128i write_mask    = _mm_castps_si128(sseWriteMask);
            const __m128i masked_output = _mm_or_si128(_mm_and_si128(write_mask, final_colour),
                                                       _mm_andnot_si128(write_mask, original_pixel_data));

            _mm_storeu_si128((__m128i *)pixel_location, masked_output);
#else
            // Mask-store 4-sample fragment values
            _mm_maskstore_epi32(
                (int *)pixel_location,
                _mm_castps_si128(sseWriteMask),
                final_colour);
#endif
        }
    }
}

void draw_object(struct arena *arena)
{
    mat4 cum_matrix = {0};
    m4_mul_m4(VP_MATRIX, raster_state.MVP, cum_matrix);

    const struct obj obj = raster_state.obj;

    if (!transformed_vertices)
        transformed_vertices = arena_alloc_aligned(arena, sizeof(vec4) * obj.num_pos, 16);

    float *pos = obj.pos;

#pragma omp parallel
    {
        //  transform all verts to screen space
#pragma omp parallel for
        for (size_t i = 0; i < obj.num_pos; i++)
        {
            vec4 vertex = {pos[i * 3 + 0],
                           pos[i * 3 + 1],
                           pos[i * 3 + 2],
                           1.0f};
            m4_mul_v4(cum_matrix, vertex, &transformed_vertices[i * 4]);
        }

        float         *obj_pos     = transformed_vertices;
        float         *obj_tex     = obj.texs;
        struct vertindices *obj_indices = obj.indices;

#pragma omp for
        for (size_t i = 0; i < obj.num_f_rows; ++i)
        {
            vec4 vertices[3]   = {0};
            vec2 tex_coords[3] = {0};

            for (size_t j = 0; j < 3; ++j)
            {
                struct vertindices indices = obj_indices[i * 3 + j];

                const int pos_index = indices.v_idx;
                const int tex_index = indices.vt_idx;
#if 0
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