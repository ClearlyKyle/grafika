#ifndef __AVX_H__
#define __AVX_H__

#include <xmmintrin.h>
#include <immintrin.h>

#include "common.h"

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

static float *transformed_vertices = NULL;

static inline bool Tie_Breaker_AB_Test(const vec3 E)
{
    if (E[0] > 0.0f)
        return true;
    if (E[0] < 0.0f)
        return false;
    if (E[0] == 0.0f && E[1] < 0.0f)
        return false;
    return true;
}

static inline bool Edge_Tie_Breaker(const float edge_value, const bool AB_test_result)
{
    // Apply tie-breaking rules on shared vertices in order to avoid double-shading fragments
    if (edge_value > 0.0f)
        return true;
    if (edge_value < 0.0f)
        return false;

    return AB_test_result;
}

static inline __m128 Interpolate_Vertex_Value(const __m128 F[3], const __m128 V[3])
{
    return _mm_add_ps(
        _mm_mul_ps(V[2], F[2]),
        _mm_add_ps(
            _mm_mul_ps(V[1], F[1]),
            _mm_mul_ps(V[0], F[0])));
}

static void draw_triangle(vec4 verts[3], vec2 texcoords[3])
{
    vec4 ss_v[3] = {0}; // screen space vertices
    for (int i = 0; i < 3; ++i)
    {
        ss_v[i][0] = verts[i][0];
        ss_v[i][1] = verts[i][1];
        ss_v[i][2] = verts[i][2];
        ss_v[i][3] = verts[i][3];
    }

    /* From the paper, calculate out "A" matrix */
    const float a0 = (ss_v[1][1] * ss_v[2][3]) - (ss_v[2][1] * ss_v[1][3]);
    const float a1 = (ss_v[2][1] * ss_v[0][3]) - (ss_v[0][1] * ss_v[2][3]);
    const float a2 = (ss_v[0][1] * ss_v[1][3]) - (ss_v[1][1] * ss_v[0][3]);

    const float b0 = (ss_v[2][0] * ss_v[1][3]) - (ss_v[1][0] * ss_v[2][3]);
    const float b1 = (ss_v[0][0] * ss_v[2][3]) - (ss_v[2][0] * ss_v[0][3]);
    const float b2 = (ss_v[1][0] * ss_v[0][3]) - (ss_v[0][0] * ss_v[1][3]);

    const float c0 = (ss_v[1][0] * ss_v[2][1]) - (ss_v[2][0] * ss_v[1][1]);
    const float c1 = (ss_v[2][0] * ss_v[0][1]) - (ss_v[0][0] * ss_v[2][1]);
    const float c2 = (ss_v[0][0] * ss_v[1][1]) - (ss_v[1][0] * ss_v[0][1]);

    const float detM = (c0 * ss_v[0][3]) + (c1 * ss_v[1][3]) + (c2 * ss_v[2][3]);
    /*
    The sign of the determinant gives the orientation of the triangle: a counterclockwise order gives a positive determinant
    */
    if (detM >= 0.0f)
        return;

    // Set up edge functions (this is A = adj(M))
    vec3 E0 = {a0, b0, c0};
    vec3 E1 = {a1, b1, c1};
    vec3 E2 = {a2, b2, c2};

    /* Projection */
    // projected = vert / vert.w;
    vec3 proj[3] = {0};
    for (int i = 0; i < 3; i++)
    {
        proj[i][0] = ss_v[i][0] / ss_v[i][3];
        proj[i][1] = ss_v[i][1] / ss_v[i][3];
        proj[i][2] = ss_v[i][2] / ss_v[i][3];
    }

    /* Get the bounding box of the triangle */
    float fminX = proj[0][0];
    float fminY = proj[0][1];
    float fmaxX = proj[0][0];
    float fmaxY = proj[0][1];

    for (int i = 1; i < 3; ++i)
    {
        /* Update minimum and maximum values */
        if (proj[i][0] < fminX)
            fminX = proj[i][0];
        if (proj[i][1] < fminY)
            fminY = proj[i][1];
        if (proj[i][0] > fmaxX)
            fmaxX = proj[i][0];
        if (proj[i][1] > fmaxY)
            fmaxY = proj[i][1];
    }

    /* Precompute constants */
    const int screen_width_minus_one  = GRAFIKA_SCREEN_WIDTH - 1;
    const int screen_height_minus_one = GRAFIKA_SCREEN_HEIGHT - 1;

    /* Clamp values to valid range */
    int minX = max(0, min((int)fminX, screen_width_minus_one));
    int minY = max(0, min((int)fminY, screen_height_minus_one));
    int maxX = max(0, min((int)fmaxX, screen_width_minus_one));
    int maxY = max(0, min((int)fmaxY, screen_height_minus_one));

    // align to multiples of 4 for better simd'ing
    minX = minX - (minX & (4 - 1)); // p & (a - 1) = p % a
    maxX = maxX + (minX & (4 - 1));

    // Setup texture coordinates
    __m128 U[3], V[3], Z[3];
    U[0] = _mm_set1_ps(texcoords[0][0]);
    U[1] = _mm_set1_ps(texcoords[1][0]);
    U[2] = _mm_set1_ps(texcoords[2][0]);

    V[0] = _mm_set1_ps(texcoords[0][1]);
    V[1] = _mm_set1_ps(texcoords[1][1]);
    V[2] = _mm_set1_ps(texcoords[2][1]);

    Z[0] = _mm_set1_ps(ss_v[0][2]);
    Z[1] = _mm_set1_ps(ss_v[1][2]);
    Z[2] = _mm_set1_ps(ss_v[2][2]);

    __m128 Edge0_A = _mm_set_ps1(E0[0]);
    __m128 Edge1_A = _mm_set_ps1(E1[0]);
    __m128 Edge2_A = _mm_set_ps1(E2[0]);

    __m128 Edge0_B = _mm_set_ps1(E0[1]);
    __m128 Edge1_B = _mm_set_ps1(E1[1]);
    __m128 Edge2_B = _mm_set_ps1(E2[1]);

    // Compute E(x, y) = (x * a) + (y * b) + c, at block origin once
    __m128 E0_origin = _mm_set1_ps((E0[0] * (float)minX) + (E0[1] * (float)minY) + E0[2]);
    __m128 E1_origin = _mm_set1_ps((E1[0] * (float)minX) + (E1[1] * (float)minY) + E1[2]);
    __m128 E2_origin = _mm_set1_ps((E2[0] * (float)minX) + (E2[1] * (float)minY) + E2[2]);

    // Generate masks used for tie-breaking rules (not to double-shade along shared edges)
    // We are only checking for true conditions
    __m128 Edge0_AB_Tie_Breaker_Result = _mm_or_ps(
        _mm_cmpgt_ps(Edge0_A, _mm_setzero_ps()), /* (E[0] > 0.0f) -> true */
        _mm_and_ps(
            _mm_cmpeq_ps(Edge0_A, _mm_setzero_ps()), /* (E[0] == 0.0f && E[1] > 0.0f) -> true */
            _mm_cmpge_ps(Edge0_B, _mm_setzero_ps())));

    __m128 Edge1_AB_Tie_Breaker_Result = _mm_or_ps(_mm_cmpgt_ps(Edge1_A, _mm_setzero_ps()),
                                                   _mm_and_ps(_mm_cmpge_ps(Edge1_B, _mm_setzero_ps()), _mm_cmpeq_ps(Edge1_A, _mm_setzero_ps())));

    __m128 Edge2_AB_Tie_Breaker_Result = _mm_or_ps(_mm_cmpgt_ps(Edge2_A, _mm_setzero_ps()),
                                                   _mm_and_ps(_mm_cmpge_ps(Edge2_B, _mm_setzero_ps()), _mm_cmpeq_ps(Edge2_A, _mm_setzero_ps())));

    __m128 starting_x = _mm_set_ps(3.5f, 2.5f, 1.5f, 0.5f);
    __m128 step_y     = _mm_set1_ps(0.5f);

    // pre fetch tex data
    const int      texw    = state.tex.w;
    const int      texh    = state.tex.h;
    const int      texbpp  = state.tex.bpp;
    unsigned char *texdata = state.tex.data;

    float *pDepthBuffer = rend.depth_buffer;

    // Start rasterizing by looping over pixels to output a per-pixel color
    for (int pix_y = minY; pix_y <= maxY;
         pix_y++,
             step_y = _mm_add_ps(step_y, _mm_set1_ps(1.0f))) /* from min to max pixels in the y */
    {
        __m128 step_x = starting_x;

        for (int pix_x = minX; pix_x <= maxX;
             pix_x += 4,
                 step_x = _mm_add_ps(step_x, _mm_set1_ps(4.0f)))
        {
            // We are stepping 4 pixels at a time in the x direction
            // a * s
            __m128 Edge0_as = _mm_mul_ps(Edge0_A, step_x);
            __m128 Edge1_as = _mm_mul_ps(Edge1_A, step_x);
            __m128 Edge2_as = _mm_mul_ps(Edge2_A, step_x);

            // b * t
            __m128 Edge0_bt = _mm_mul_ps(Edge0_B, step_y);
            __m128 Edge1_bt = _mm_mul_ps(Edge1_B, step_y);
            __m128 Edge2_bt = _mm_mul_ps(Edge2_B, step_y);

            // E(x + s, y + t) = E(x, y) + (a * s) + (t * b)
            __m128 Edge_Func0 = _mm_add_ps(E0_origin, _mm_add_ps(Edge0_as, Edge0_bt));
            __m128 Edge_Func1 = _mm_add_ps(E1_origin, _mm_add_ps(Edge1_as, Edge1_bt));
            __m128 Edge_Func2 = _mm_add_ps(E2_origin, _mm_add_ps(Edge2_as, Edge2_bt));

            /*
                Apply tie-breaking rules on shared vertices in order to avoid double-shading fragments
                if (edge_value > 0.0f)
                    return true;
                if (edge_value < 0.0f)
                    return false;
                if (E[0] > 0.0f)
                    return true;
                if (E[0] < 0.0f)
                    return false;
                if (E[0] == 0.0f && E[1] < 0.0f)
                    return false;
                return true;
            */

#if 1 /* Tie Breaker Edge testing */
            // Edge 0 test
            __m128 Edge0Positive = _mm_cmplt_ps(Edge_Func0, _mm_setzero_ps());
            __m128 Edge0Negative = _mm_cmpgt_ps(Edge_Func0, _mm_setzero_ps());
            __m128 Edge0FuncMask = _mm_or_ps(Edge0Positive,
                                             _mm_andnot_ps(Edge0Negative,
                                                           Edge0_AB_Tie_Breaker_Result));
            // Edge 1 test
            __m128 Edge1Positive = _mm_cmplt_ps(Edge_Func1, _mm_setzero_ps());
            __m128 Edge1Negative = _mm_cmpgt_ps(Edge_Func1, _mm_setzero_ps());
            __m128 Edge1FuncMask = _mm_or_ps(Edge1Positive,
                                             _mm_andnot_ps(Edge1Negative,
                                                           Edge1_AB_Tie_Breaker_Result));
            // Edge 2 test
            __m128 Edge2Positive = _mm_cmplt_ps(Edge_Func2, _mm_setzero_ps());
            __m128 Edge2Negative = _mm_cmpgt_ps(Edge_Func2, _mm_setzero_ps());
            __m128 Edge2FuncMask = _mm_or_ps(Edge2Positive,
                                             _mm_andnot_ps(Edge2Negative,
                                                           Edge2_AB_Tie_Breaker_Result));
#else
            __m128 Edge0FuncMask = _mm_cmplt_ps(Edge_Func0, _mm_setzero_ps());
            __m128 Edge1FuncMask = _mm_cmplt_ps(Edge_Func1, _mm_setzero_ps());
            __m128 Edge2FuncMask = _mm_cmplt_ps(Edge_Func2, _mm_setzero_ps());
#endif
            // Combine resulting masks of all three edges
            __m128 EdgeFuncTestResult = _mm_and_ps(Edge0FuncMask, _mm_and_ps(Edge1FuncMask, Edge2FuncMask));

            uint16_t maskInt = (uint16_t)_mm_movemask_ps(EdgeFuncTestResult);

            if (maskInt == 0x0)
                continue;

            /* Calculate Interpolation Coefficients */
            __m128 F[3] = {0};
            F[0]        = Edge_Func0;
            F[1]        = Edge_Func1;
            F[2]        = Edge_Func2;

            // R(x, y) = F0(x, y) + F1(x, y) + F2(x, y)
            // r = 1 / (F0(x, y) + F1(x, y) + F2(x, y))
            __m128 r = _mm_add_ps(F[0], _mm_add_ps(F[1], F[2]));
            r        = _mm_rcp_ps(r);

            F[0] = _mm_mul_ps(F[0], r);
            F[1] = _mm_mul_ps(F[1], r);
            F[2] = _mm_mul_ps(F[2], r);

            const size_t pixel_index = (const size_t)(pix_y * GRAFIKA_SCREEN_WIDTH + pix_x);

            // Interpolate depth values prior to depth test
            const __m128 sseZInterpolated = Interpolate_Vertex_Value(F, Z);

            // Load current depth buffer contents
            // float       *pDepth          = pDepthBuffer + pixel_index;
            float       *pDepth          = &pDepthBuffer[pixel_index];
            const __m128 sseDepthCurrent = _mm_loadu_ps(pDepth);

            // Perform LESS_THAN_EQUAL depth test
            const __m128 sseDepthRes = _mm_cmple_ps(sseZInterpolated, sseDepthCurrent);

            if ((uint16_t)_mm_movemask_ps(sseDepthRes) == 0x0)
                continue;

            // AND depth mask & coverage mask for quads of fragments
            const __m128 sseWriteMask = _mm_and_ps(sseDepthRes, EdgeFuncTestResult);

            // Write interpolated Z values
            _mm_maskmoveu_si128(
                _mm_castps_si128(sseZInterpolated),
                _mm_castps_si128(sseWriteMask),
                (char *)pDepth);

            /* Interpolate Texture Coordinates */
            __m128 pixU = Interpolate_Vertex_Value(F, U);
            __m128 pixV = Interpolate_Vertex_Value(F, V);

            // clamp the vector to the range [0.0f, 1.0f]
            pixU = _mm_max_ps(_mm_min_ps(pixU, _mm_set1_ps(1.0f)), _mm_setzero_ps());
            pixV = _mm_max_ps(_mm_min_ps(pixV, _mm_set1_ps(1.0f)), _mm_setzero_ps());

            pixU = _mm_mul_ps(pixU, _mm_set1_ps((float)texw));
            pixV = _mm_mul_ps(pixV, _mm_set1_ps((float)texh));

            // (U + tex.width * V) * tex.bpp
            __m128i texOffset = _mm_mullo_epi32(
                _mm_set1_epi32(texbpp),
                _mm_add_epi32(
                    _mm_cvtps_epi32(pixU),
                    _mm_mullo_epi32(
                        _mm_set1_epi32(texw),
                        _mm_cvtps_epi32(pixV))));

#if 0
            uint32_t offset[4] = {0};
            offset[0]          = (uint32_t)_mm_extract_epi32(texOffset, 0);
            offset[1]          = (uint32_t)_mm_extract_epi32(texOffset, 1);
            offset[2]          = (uint32_t)_mm_extract_epi32(texOffset, 2);
            offset[3]          = (uint32_t)_mm_extract_epi32(texOffset, 3);
#else
            uint32_t *offset = (uint32_t *)&texOffset; // I think this is unsafe
#endif
            char *sample[4] = {0};
            for (int i = 0; i < 4; ++i)
                sample[i] = (char *)(texdata + offset[i]);

            // A B G R
            __m128i final_colour = _mm_setr_epi8(sample[0][0], sample[0][1], sample[0][2], -1,
                                                 sample[1][0], sample[1][1], sample[1][2], -1,
                                                 sample[2][0], sample[2][1], sample[2][2], -1,
                                                 sample[3][0], sample[3][1], sample[3][2], -1);

            uint32_t *pixel_location = &rend.pixels[pixel_index];

#if 1 /* Fabian method */
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

static void draw_onexit(void)
{
    SAFE_FREE(transformed_vertices);
}

static void draw_object(void)
{
    mat4 cumMatrix = {0};
    m4_mul_m4(VP_MATRIX, state.MVP, cumMatrix);

    const obj_t object = state.obj;

    if (!transformed_vertices)
        transformed_vertices = malloc(sizeof(vec4) * object.num_pos);

    float *pos = object.pos;
#pragma omp parallel for
    // transform all verts to screen space
    for (size_t i = 0; i < object.num_pos; i++)
    {
        vec4 vertex = {pos[i * 3 + 0],
                       pos[i * 3 + 1],
                       pos[i * 3 + 2],
                       1.0f};
        m4_mul_v4(cumMatrix, vertex, &transformed_vertices[i * 4]);
    }

#pragma omp parallel
    {
        float         *pPos     = transformed_vertices;
        float         *pTex     = object.texs;
        vertindices_t *pIndices = object.indices;

#pragma omp for
        for (size_t i = 0; i < object.num_f_rows; ++i)
        {
            vec4 vertices[3]  = {0};
            vec2 texcoords[3] = {0};

            for (size_t j = 0; j < 3; ++j)
            {
                vertindices_t indices = pIndices[i * 3 + j];

                const int posIndex = indices.v_idx;
                const int texIndex = indices.vt_idx;
#if 0
                const size_t storeIndex = j; // CCW triangles
#else
                const size_t storeIndex = 2 - j; // CW triangles (which blender uses)
#endif
                float *pP               = pPos + 4 * posIndex;
                vertices[storeIndex][0] = pP[0];
                vertices[storeIndex][1] = pP[1];
                vertices[storeIndex][2] = pP[2];
                vertices[storeIndex][3] = pP[3];

                texcoords[storeIndex][0] = pTex[2 * texIndex + 0];
                texcoords[storeIndex][1] = pTex[2 * texIndex + 1];
            }
            draw_triangle(vertices, texcoords);
        }
    }
}

#endif // __AVX_H__