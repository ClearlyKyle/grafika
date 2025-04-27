#ifndef __STEPPING_H__
#define __STEPPING_H__

#include "common.h"

// static void draw_triangle(const struct triangle t)
//{
//     // convert to clip space
//     vec4 clipspace[3] = {0};
//     for (int i = 0; i < 3; ++i)
//     {
//         vec4 pos = {t.pos[i][0], t.pos[i][1], t.pos[i][2], 1.0f};
//         m4_mul_v4(raster_state.MVP, pos, clipspace[i]);

//        // clipping
//        if (clipspace[i][0] < -clipspace[i][3] || clipspace[i][0] > clipspace[i][3] ||
//            clipspace[i][1] < -clipspace[i][3] || clipspace[i][1] > clipspace[i][3] ||
//            clipspace[i][2] < -clipspace[i][3] || clipspace[i][2] > clipspace[i][3])
//        {
//            return; // triangle is outside
//        }
//    }

//    // perspective division (clip to ndc)
//    vec3 ndc[3] = {0}, w_vals = {0};
//    for (size_t i = 0; i < 3; i++)
//    {
//        w_vals[i] = 1.0f / clipspace[i][3]; // 1.0f / w
//        ndc[i][0] = clipspace[i][0] * w_vals[i];
//        ndc[i][1] = clipspace[i][1] * w_vals[i];
//        ndc[i][2] = clipspace[i][2] * w_vals[i];
//    }

//    // back face culling (surface normal)
//    vec3 sub10 = {0}, sub20 = {0}, normal = {0};
//    v3_sub(ndc[1], ndc[0], sub10);
//    v3_sub(ndc[2], ndc[0], sub20);
//    v3_cross(sub10, sub20, normal);
//    if (normal[2] > 0.0f)
//        return;

//    vec3 screenspace[3] = {0};
//    for (int i = 0; i < 3; ++i)
//    {
//        screenspace[i][0] = (ndc[i][0] + 1.0f) * 0.5f * GRAFIKA_SCREEN_WIDTH;
//        screenspace[i][1] = (ndc[i][1] + 1.0f) * 0.5f * GRAFIKA_SCREEN_HEIGHT;
//        screenspace[i][2] = (ndc[i][2] + 1.0f) * 0.5f;
//    }

//    // calculate bounding rectangle
//    int AABB[4] = {0};
//    AABB_make(screenspace, AABB);

//    AABB[0] = AABB[0] & ~3;       // round down AABB[0] to the nearest multiple of 4
//    AABB[2] = (AABB[2] + 3) & ~3; // round up AABB[2] to the next multiple of 4 - 1

//    ASSERT(AABB[0] % 4 == 0, "starting x is not multiple of 4\n");
//    ASSERT(AABB[2] % 4 == 0, "ending x is not multiple of 4\n");

//    const float dY0 = screenspace[2][1] - screenspace[1][1], dX0 = screenspace[1][0] - screenspace[2][0];
//    const float dY1 = screenspace[0][1] - screenspace[2][1], dX1 = screenspace[2][0] - screenspace[0][0];
//    const float dY2 = screenspace[1][1] - screenspace[0][1], dX2 = screenspace[0][0] - screenspace[1][0];

//    const float C0 = (screenspace[2][0] * screenspace[1][1]) - (screenspace[2][1] * screenspace[1][0]);
//    const float C1 = (screenspace[0][0] * screenspace[2][1]) - (screenspace[0][1] * screenspace[2][0]);
//    const float C2 = (screenspace[1][0] * screenspace[0][1]) - (screenspace[1][1] * screenspace[0][0]);

//    float alpha0 = (dY0 * ((float)AABB[0] + 0.5f)) + (dX0 * ((float)AABB[1])) + C0;
//    float alpha1 = (dY0 * ((float)AABB[0] + 1.5f)) + (dX0 * ((float)AABB[1])) + C0;
//    float alpha2 = (dY0 * ((float)AABB[0] + 2.5f)) + (dX0 * ((float)AABB[1])) + C0;
//    float alpha3 = (dY0 * ((float)AABB[0] + 3.5f)) + (dX0 * ((float)AABB[1])) + C0;

//    float betaa0 = (dY1 * ((float)AABB[0] + 0.5f)) + (dX1 * ((float)AABB[1])) + C1;
//    float betaa1 = (dY1 * ((float)AABB[0] + 1.5f)) + (dX1 * ((float)AABB[1])) + C1;
//    float betaa2 = (dY1 * ((float)AABB[0] + 2.5f)) + (dX1 * ((float)AABB[1])) + C1;
//    float betaa3 = (dY1 * ((float)AABB[0] + 3.5f)) + (dX1 * ((float)AABB[1])) + C1;

//    float gamma0 = (dY2 * ((float)AABB[0] + 0.5f)) + (dX2 * ((float)AABB[1])) + C2;
//    float gamma1 = (dY2 * ((float)AABB[0] + 1.5f)) + (dX2 * ((float)AABB[1])) + C2;
//    float gamma2 = (dY2 * ((float)AABB[0] + 2.5f)) + (dX2 * ((float)AABB[1])) + C2;
//    float gamma3 = (dY2 * ((float)AABB[0] + 3.5f)) + (dX2 * ((float)AABB[1])) + C2;

//    // const float inv_area = 1.0f / (dX1 * dY2 - dY1 * dX2);
//    //  screenspace[1][2] = (screenspace[1][2] - screenspace[0][2]) * inv_area;
//    //  screenspace[2][2] = (screenspace[2][2] - screenspace[0][2]) * inv_area;

//    // const float zstep = dY1 * screenspace[1][2] + dY2 * screenspace[2][2];

//    // pre fetch tex data
//    const int      texw    = raster_state.tex.w;
//    const int      texh    = raster_state.tex.h;
//    const int      texbpp  = raster_state.tex.bpp;
//    unsigned char *texdata = raster_state.tex.data;

//    // float *depth_buffer = rend.depth_buffer;

//    // Store edge 0 equation coefficients
//    const __m128 A0 = _mm_set1_ps(dY0);
//    const __m128 A1 = _mm_set1_ps(dY1);
//    const __m128 A2 = _mm_set1_ps(dY2);

//    const __m128 B0 = _mm_set1_ps(dX0);
//    const __m128 B1 = _mm_set1_ps(dX1);
//    const __m128 B2 = _mm_set1_ps(dX2);

//    __m128 E0_origin = _mm_setr_ps(alpha0, alpha1, alpha2, alpha3);
//    __m128 E1_origin = _mm_setr_ps(betaa0, betaa1, betaa2, betaa3);
//    __m128 E2_origin = _mm_setr_ps(gamma0, gamma1, gamma2, gamma3);

//    const __m128 mm_w_vals0 = _mm_set1_ps(w_vals[0]);
//    const __m128 mm_w_vals1 = _mm_set1_ps(w_vals[1]);
//    const __m128 mm_w_vals2 = _mm_set1_ps(w_vals[2]);

//    // const __m128 triArea = _mm_mul_ps(_mm_mul_ps(B2, A1), _mm_mul_ps(B1, A2));
//    //// const __m128  oneOverTriArea = _mm_rcp_ps(_mm_cvtepi32_ps(triArea));
//    // const __m128 oneOverTriArea = _mm_div_ps(_mm_set1_ps(1.0f), triArea);

//    const float  triArea        = dX1 * dY2 - dY1 * dX2;
//    const __m128 oneOverTriArea = _mm_set1_ps(1.0f / triArea);

//    const __m128 x_step = _mm_set1_ps(4.0f);
//    __m128       y_step = _mm_set1_ps(1.0f);

//    // rasterize
//    for (int y         = AABB[1]; y < AABB[3]; ++y,
//             E0_origin = _mm_add_ps(E0_origin, B0),
//             E1_origin = _mm_add_ps(E1_origin, B1),
//             E2_origin = _mm_add_ps(E2_origin, B2))
//    {
//        // Barycentric coordinates at start of row
//        // float alpha_start = alpha;
//        // float beta_start  = betaa;
//        // float gamma_start = gamma;

//        //__m128 step_x = starting_x;

//        __m128 w0 = E0_origin;
//        __m128 w1 = E1_origin;
//        __m128 w2 = E2_origin;

//        // float depth_start = screenspace[0][2] + (screenspace[1][2] * betaa) + (screenspace[2][2] * gamma);

//        for (int x  = AABB[0]; x < AABB[2]; x += 4,
//                 w0 = _mm_add_ps(w0, _mm_mul_ps(A0, x_step)),
//                 w1 = _mm_add_ps(w1, _mm_mul_ps(A1, x_step)),
//                 w2 = _mm_add_ps(w2, _mm_mul_ps(A2, x_step)))
//        {
//            // Compute barycentric coordinates for 4 pixels
//            //__m128 fx = _mm_set_ps(x + 3.0f, x + 2.0f, x + 1.0f, x + 0.0f);
//            //__m128 fy = _mm_set1_ps((float)y);

//            ////__m128 w0 = _mm_add_ps(_mm_set1_ps(alpha_start), _mm_mul_ps(_mm_set1_ps(dY0), fx));
//            ////__m128 w1 = _mm_add_ps(_mm_set1_ps(beta_start), _mm_mul_ps(_mm_set1_ps(dY1), fx));
//            ////__m128 w2 = _mm_add_ps(_mm_set1_ps(gamma_start), _mm_mul_ps(_mm_set1_ps(dY2), fx));

//            //// Compute barycentric coordinates
//            //__m128 w0 = _mm_add_ps(_mm_mul_ps(_mm_set1_ps(dY0), fx), _mm_mul_ps(_mm_set1_ps(dX0), fy));
//            // w0        = _mm_add_ps(w0, _mm_set1_ps(C0));
//            //__m128 w1 = _mm_add_ps(_mm_mul_ps(_mm_set1_ps(dY1), fx), _mm_mul_ps(_mm_set1_ps(dX1), fy));
//            // w1        = _mm_add_ps(w1, _mm_set1_ps(C1));
//            //__m128 w2 = _mm_add_ps(_mm_mul_ps(_mm_set1_ps(dY2), fx), _mm_mul_ps(_mm_set1_ps(dX2), fy));
//            // w2        = _mm_add_ps(w2, _mm_set1_ps(C2));

//            //// check if inside triangle: (w0, w1, w2) >= 0
//            //__m128 mask0 = _mm_cmplt_ps(w0, _mm_setzero_ps());
//            //__m128 mask1 = _mm_cmplt_ps(w1, _mm_setzero_ps());
//            //__m128 mask2 = _mm_cmplt_ps(w2, _mm_setzero_ps());

//            //__m128 inside_mask = _mm_or_ps(_mm_or_ps(mask0, mask1), mask2);

//            // if (_mm_movemask_ps(inside_mask) != 0)
//            //     continue; // All outside

//            // We are stepping 4 pixels at a time in the x direction
//            // a * s
//            // const __m128 Edge0_as = _mm_mul_ps(A0, step_x);
//            // const __m128 Edge1_as = _mm_mul_ps(A1, step_x);
//            // const __m128 Edge2_as = _mm_mul_ps(A2, step_x);

//            //// b * t
//            // const __m128 Edge0_bt = _mm_mul_ps(B0, step_y);
//            // const __m128 Edge1_bt = _mm_mul_ps(B1, step_y);
//            // const __m128 Edge2_bt = _mm_mul_ps(B2, step_y);

//            // E(x + s, y + t) = E(x, y) + (a * s) + (t * b)
//            // const __m128 Edge_Func0 = _mm_add_ps(E0_origin, _mm_add_ps(Edge0_as, Edge0_bt));
//            // const __m128 Edge_Func1 = _mm_add_ps(E1_origin, _mm_add_ps(Edge1_as, Edge1_bt));
//            // const __m128 Edge_Func2 = _mm_add_ps(E2_origin, _mm_add_ps(Edge2_as, Edge2_bt));

//            const __m128 Edge0FuncMask = _mm_cmpgt_ps(w0, _mm_setzero_ps());
//            const __m128 Edge1FuncMask = _mm_cmpgt_ps(w1, _mm_setzero_ps());
//            const __m128 Edge2FuncMask = _mm_cmpgt_ps(w2, _mm_setzero_ps());

//            // Combine resulting masks of all three edges
//            const __m128 EdgeFuncTestResult = _mm_and_ps(Edge0FuncMask, _mm_and_ps(Edge1FuncMask, Edge2FuncMask));

//            const uint16_t maskInt = (uint16_t)_mm_movemask_ps(EdgeFuncTestResult);

//            // if (maskInt == 0x0)
//            //     continue;

//            if (maskInt && 0x000F)
//            {
//                grafika_setpixel(x + 0, y, 0xFF112233);
//            }
//            if (maskInt && 0x00F0)
//            {
//                grafika_setpixel(x + 1, y, 0xFF112233);
//            }
//            if (maskInt && 0x0F00)
//            {
//                grafika_setpixel(x + 2, y, 0xFF112233);
//            }
//            if (maskInt && 0xF000)
//            {
//                grafika_setpixel(x + 3, y, 0xFF112233);
//            }

//            //// depth test
//            //__m128 depth_new = _mm_add_ps(_mm_set1_ps(depth_start), _mm_mul_ps(_mm_set1_ps(zstep), fx));

//            // const int    index         = (x * GRAFIKA_SCREEN_WIDTH) + y;
//            // float *const depth_current = depth_buffer + index;

//            // const __m128 depth_old  = _mm_loadu_ps(depth_current);
//            // const __m128 depth_test = _mm_cmplt_ps(depth_new, depth_old);

//            // if (_mm_movemask_ps(depth_test) == 0)
//            //     continue;

//            // const __m128 depth_write_mask = _mm_and_ps(depth_new, inside_mask);

//            // const __m128 finaldepth = _mm_blendv_ps(depth_old, depth_new, depth_write_mask);
//            //_mm_store_ps(depth_current, finaldepth);

//            //// Weights, see what i did there ;)
//            // const float wait0 = w0 * w_vals[0];
//            // const float wait1 = w1 * w_vals[1];
//            // const float wait2 = w2 * w_vals[2];

//            //__m128 wait0 = _mm_mul_ps(w0, oneOverTriArea); // Bary A
//            //__m128 wait1 = _mm_mul_ps(w1, oneOverTriArea); // B
//            //__m128 wait2 = _mm_mul_ps(w2, oneOverTriArea); // C

//            // wait0 = _mm_mul_ps(wait0, mm_w_vals0);
//            // wait1 = _mm_mul_ps(wait1, mm_w_vals1);
//            // wait2 = _mm_mul_ps(wait2, mm_w_vals2);

//            //// calculate correction factor
//            //// R(x, y) = F0(x, y) + F1(x, y) + F2(x, y)
//            //// r = 1 / (F0(x, y) + F1(x, y) + F2(x, y))

//            //// const float cf = 1.0f / (wait0 + wait1 + wait2);
//            // const __m128 cf = _mm_rcp_ps(_mm_add_ps(_mm_add_ps(wait0, wait1), wait2));

//            // wait0 = _mm_mul_ps(wait0, cf);
//            // wait1 = _mm_mul_ps(wait1, cf);
//            // wait2 = _mm_mul_ps(wait2, cf);

//            // float u0 = (t.tex[0][0] * wait0[0] + t.tex[1][0] * wait1[0] + t.tex[2][0] * wait2[0]);
//            // float u1 = (t.tex[0][0] * wait0[1] + t.tex[1][0] * wait1[1] + t.tex[2][0] * wait2[1]);
//            // float u2 = (t.tex[0][0] * wait0[2] + t.tex[1][0] * wait1[2] + t.tex[2][0] * wait2[2]);
//            // float u3 = (t.tex[0][0] * wait0[3] + t.tex[1][0] * wait1[3] + t.tex[2][0] * wait2[3]);

//            // float v0 = (t.tex[0][1] * wait0[0] + t.tex[1][1] * wait1[0] + t.tex[2][1] * wait2[0]);
//            // float v1 = (t.tex[0][1] * wait0[1] + t.tex[1][1] * wait1[1] + t.tex[2][1] * wait2[1]);
//            // float v2 = (t.tex[0][1] * wait0[2] + t.tex[1][1] * wait1[2] + t.tex[2][1] * wait2[2]);
//            // float v3 = (t.tex[0][1] * wait0[3] + t.tex[1][1] * wait1[3] + t.tex[2][1] * wait2[3]);

//            // if (u0 > 0.0f && v0 > 0.0f)
//            //{
//            //     u0 = SDL_clamp(u0, 0.0f, 1.0f);
//            //     v0 = SDL_clamp(v0, 0.0f, 1.0f);

//            //    float new_u = u0 * (float)(texw - 1);
//            //    float new_v = v0 * (float)(texh - 1);

//            //    unsigned char *texcolour0 = texdata + (((int)new_u + texw * (int)new_v) * texbpp);

//            //    uint32_t pixelcolour0 = ((uint32_t)texcolour0[0] << 16) | ((uint32_t)texcolour0[1] << 8) | (uint32_t)texcolour0[2];
//            //    grafika_setpixel(x + 0, y, pixelcolour0);
//            //}
//            // if (u1 > 0.0f && v1 > 0.0f)
//            //{
//            //    u1 = SDL_clamp(u1, 0.0f, 1.0f);
//            //    v1 = SDL_clamp(v1, 0.0f, 1.0f);

//            //    float new_u = u1 * (float)(texw - 1);
//            //    float new_v = v1 * (float)(texh - 1);

//            //    unsigned char *texcolour0 = texdata + (((int)new_u + texw * (int)new_v) * texbpp);

//            //    uint32_t pixelcolour0 = ((uint32_t)texcolour0[0] << 16) | ((uint32_t)texcolour0[1] << 8) | (uint32_t)texcolour0[2];
//            //    grafika_setpixel(x + 1, y, pixelcolour0);
//            //}
//            // if (u2 > 0.0f && v2 > 0.0f)
//            //{
//            //    u2 = SDL_clamp(u2, 0.0f, 1.0f);
//            //    v2 = SDL_clamp(v2, 0.0f, 1.0f);

//            //    float new_u = u2 * (float)(texw - 1);
//            //    float new_v = v2 * (float)(texh - 1);

//            //    unsigned char *texcolour0 = texdata + (((int)new_u + texw * (int)new_v) * texbpp);

//            //    uint32_t pixelcolour0 = ((uint32_t)texcolour0[0] << 16) | ((uint32_t)texcolour0[1] << 8) | (uint32_t)texcolour0[2];
//            //    grafika_setpixel(x + 2, y, pixelcolour0);
//            //}
//            // if (u3 > 0.0f && v3 > 0.0f)
//            //{
//            //    u3 = SDL_clamp(u3, 0.0f, 1.0f);
//            //    v3 = SDL_clamp(v3, 0.0f, 1.0f);

//            //    float new_u = u3 * (float)(texw - 1);
//            //    float new_v = v3 * (float)(texh - 1);

//            //    unsigned char *texcolour0 = texdata + (((int)new_u + texw * (int)new_v) * texbpp);

//            //    uint32_t pixelcolour0 = ((uint32_t)texcolour0[0] << 16) | ((uint32_t)texcolour0[1] << 8) | (uint32_t)texcolour0[2];
//            //    grafika_setpixel(x + 3, y, pixelcolour0);
//            //}
//        }
//    }
//}
static void draw_triangle(const struct triangle t)
{
    // convert to clip space
    vec4 clipspace[3] = {0};
    for (int i = 0; i < 3; ++i)
    {
        vec4 pos = {t.pos[i][0], t.pos[i][1], t.pos[i][2], 1.0f};
        m4_mul_v4(raster_state.MVP, pos, clipspace[i]);

        // clipping
        if (clipspace[i][0] < -clipspace[i][3] || clipspace[i][0] > clipspace[i][3] ||
            clipspace[i][1] < -clipspace[i][3] || clipspace[i][1] > clipspace[i][3] ||
            clipspace[i][2] < -clipspace[i][3] || clipspace[i][2] > clipspace[i][3])
        {
            return; // triangle is outside
        }
    }

    // perspective division (clip to ndc)
    vec3 ndc[3] = {0}, w_vals = {0};
    for (size_t i = 0; i < 3; i++)
    {
        w_vals[i] = 1.0f / clipspace[i][3]; // 1.0f / w
        ndc[i][0] = clipspace[i][0] * w_vals[i];
        ndc[i][1] = clipspace[i][1] * w_vals[i];
        ndc[i][2] = clipspace[i][2] * w_vals[i];
    }

    // back face culling (surface normal)
    vec3 sub10 = {0}, sub20 = {0}, normal = {0};
    v3_sub(ndc[1], ndc[0], sub10);
    v3_sub(ndc[2], ndc[0], sub20);
    v3_cross(sub10, sub20, normal);
    if (normal[2] > 0.0f)
        return;

    vec3 screenspace[3] = {0};
    for (int i = 0; i < 3; ++i)
    {
        screenspace[i][0] = (ndc[i][0] + 1.0f) * 0.5f * GRAFIKA_SCREEN_WIDTH;
        screenspace[i][1] = (ndc[i][1] + 1.0f) * 0.5f * GRAFIKA_SCREEN_HEIGHT;
        screenspace[i][2] = (ndc[i][2] + 1.0f) * 0.5f;
    }

    // calculate bounding rectangle
    int AABB[4] = {0};
    AABB_make(screenspace, AABB);

    AABB[0] = AABB[0] & ~3;       // round down AABB[0] to the nearest multiple of 4
    AABB[2] = (AABB[2] + 3) & ~3; // round up AABB[2] to the next multiple of 4 - 1

    ASSERT(AABB[0] % 4 == 0, "starting x is not multiple of 4\n");
    ASSERT(AABB[2] % 4 == 0, "ending x is not multiple of 4\n");

    const float dY0 = screenspace[2][1] - screenspace[1][1], dX0 = screenspace[1][0] - screenspace[2][0];
    const float dY1 = screenspace[0][1] - screenspace[2][1], dX1 = screenspace[2][0] - screenspace[0][0];
    const float dY2 = screenspace[1][1] - screenspace[0][1], dX2 = screenspace[0][0] - screenspace[1][0];

    const float C0 = (screenspace[2][0] * screenspace[1][1]) - (screenspace[2][1] * screenspace[1][0]);
    const float C1 = (screenspace[0][0] * screenspace[2][1]) - (screenspace[0][1] * screenspace[2][0]);
    const float C2 = (screenspace[1][0] * screenspace[0][1]) - (screenspace[1][1] * screenspace[0][0]);

    // Step vectors (for 4-pixel stride)
    const __m128 A0 = _mm_set1_ps(dY0);
    const __m128 A1 = _mm_set1_ps(dY1);
    const __m128 A2 = _mm_set1_ps(dY2);

    __m128 A0_start = _mm_mul_ps(A0, _mm_add_ps(_mm_setr_ps(0.5f, 1.5f, 2.5f, 3.5f), _mm_set1_ps((float)AABB[0] + 0.5f)));
    __m128 A1_start = _mm_mul_ps(A1, _mm_add_ps(_mm_setr_ps(0.5f, 1.5f, 2.5f, 3.5f), _mm_set1_ps((float)AABB[0] + 0.5f)));
    __m128 A2_start = _mm_mul_ps(A2, _mm_add_ps(_mm_setr_ps(0.5f, 1.5f, 2.5f, 3.5f), _mm_set1_ps((float)AABB[0] + 0.5f)));

    const __m128 B0 = _mm_set1_ps(dX0);
    const __m128 B1 = _mm_set1_ps(dX1);
    const __m128 B2 = _mm_set1_ps(dX2);

    __m128 B0_start = _mm_mul_ps(B0, _mm_set1_ps((float)AABB[1] + 0.5f));
    __m128 B1_start = _mm_mul_ps(B1, _mm_set1_ps((float)AABB[1] + 0.5f));
    __m128 B2_start = _mm_mul_ps(B2, _mm_set1_ps((float)AABB[1] + 0.5f));

    __m128 E0 = _mm_add_ps(_mm_add_ps(A0_start, B0_start), _mm_set1_ps(C0));
    __m128 E1 = _mm_add_ps(_mm_add_ps(A1_start, B1_start), _mm_set1_ps(C1));
    __m128 E2 = _mm_add_ps(_mm_add_ps(A2_start, B2_start), _mm_set1_ps(C2));

    __m128 B0_inc = B0;
    __m128 B1_inc = B1;
    __m128 B2_inc = B2;

    __m128 A0_inc = _mm_mul_ps(A0, _mm_set1_ps(4.0f));
    __m128 A1_inc = _mm_mul_ps(A1, _mm_set1_ps(4.0f));
    __m128 A2_inc = _mm_mul_ps(A2, _mm_set1_ps(4.0f));

    const __m128 mm_w_vals0 = _mm_set1_ps(w_vals[0]);
    const __m128 mm_w_vals1 = _mm_set1_ps(w_vals[1]);
    const __m128 mm_w_vals2 = _mm_set1_ps(w_vals[2]);

    // pre fetch tex data
    const int      texw    = raster_state.tex.w;
    const int      texh    = raster_state.tex.h;
    const int      texbpp  = raster_state.tex.bpp;
    unsigned char *texdata = raster_state.tex.data;

    const float inv_area = 1.0f / (dX1 * dY2 - dY1 * dX2);
    screenspace[1][2]    = (screenspace[1][2] - screenspace[0][2]) * inv_area;
    screenspace[2][2]    = (screenspace[2][2] - screenspace[0][2]) * inv_area;

    const float depth_step = dY1 * screenspace[1][2] + dY2 * screenspace[2][2];

    // rasterize
    for (int y  = AABB[1]; y <= AABB[3]; ++y,
             E0 = _mm_add_ps(E0, B0_inc),
             E1 = _mm_add_ps(E1, B1_inc),
             E2 = _mm_add_ps(E2, B2_inc))
    {
        __m128 alpha = E0;
        __m128 betaa = E1;
        __m128 gamma = E2;

        float z0 = (screenspace[0][2] + screenspace[1][2] * betaa[0] + screenspace[2][2] * gamma[0]);
        float z1 = (screenspace[0][2] + screenspace[1][2] * betaa[1] + screenspace[2][2] * gamma[1]);
        float z2 = (screenspace[0][2] + screenspace[1][2] * betaa[2] + screenspace[2][2] * gamma[2]);
        float z3 = (screenspace[0][2] + screenspace[1][2] * betaa[3] + screenspace[2][2] * gamma[3]);

        for (int x     = AABB[0]; x <= AABB[2]; x += 4,
                 alpha = _mm_add_ps(alpha, A0_inc),
                 betaa = _mm_add_ps(betaa, A1_inc),
                 gamma = _mm_add_ps(gamma, A2_inc),
                 z0 += 4 * depth_step,
                 z1 += 4 * depth_step,
                 z2 += 4 * depth_step,
                 z3 += 4 * depth_step)
        {
            __m128 Edge0FuncMask = _mm_cmpgt_ps(alpha, _mm_setzero_ps());
            __m128 Edge1FuncMask = _mm_cmpgt_ps(betaa, _mm_setzero_ps());
            __m128 Edge2FuncMask = _mm_cmpgt_ps(gamma, _mm_setzero_ps());

            __m128 EdgeFuncTestResult = _mm_and_ps(Edge0FuncMask, _mm_and_ps(Edge1FuncMask, Edge2FuncMask));

            if (_mm_movemask_ps(EdgeFuncTestResult) == 0) continue;

            __m128 wait0 = _mm_mul_ps(alpha, mm_w_vals0); // Bary A
            __m128 wait1 = _mm_mul_ps(betaa, mm_w_vals1); // B
            __m128 wait2 = _mm_mul_ps(gamma, mm_w_vals2); // C

            // wait0 = _mm_mul_ps(wait0, mm_w_vals0);
            // wait1 = _mm_mul_ps(wait1, mm_w_vals1);
            // wait2 = _mm_mul_ps(wait2, mm_w_vals2);

            const __m128 cf = _mm_rcp_ps(_mm_add_ps(_mm_add_ps(wait0, wait1), wait2));

            wait0 = _mm_mul_ps(wait0, cf);
            wait1 = _mm_mul_ps(wait1, cf);
            wait2 = _mm_mul_ps(wait2, cf);

            const size_t pixel_index = (const size_t)(y * GRAFIKA_SCREEN_WIDTH + x);

            // const float z0 = (clipspace[0][3] * wait0[0] + clipspace[1][3] * wait1[0] + clipspace[2][3] * wait2[0]);
            // const float z1 = (clipspace[0][3] * wait0[1] + clipspace[1][3] * wait1[1] + clipspace[2][3] * wait2[1]);
            // const float z2 = (clipspace[0][3] * wait0[2] + clipspace[1][3] * wait1[2] + clipspace[2][3] * wait2[2]);
            // const float z3 = (clipspace[0][3] * wait0[3] + clipspace[1][3] * wait1[3] + clipspace[2][3] * wait2[3]);

            const __m128 interpolated_z = _mm_setr_ps(z0, z1, z2, z3);

            float       *pDepth          = &rend.depth_buffer[pixel_index];
            const __m128 sseDepthCurrent = _mm_load_ps(pDepth);
            // Perform LESS_THAN_EQUAL depth test
            const __m128 sseDepthRes = _mm_cmple_ps(interpolated_z, sseDepthCurrent);

            if (_mm_movemask_ps(sseDepthRes) == 0) continue;

            const __m128 sseWriteMask = _mm_and_ps(sseDepthRes, EdgeFuncTestResult);

            // Write interpolated Z values
            _mm_maskmoveu_si128(
                _mm_castps_si128(interpolated_z),
                _mm_castps_si128(sseWriteMask),
                (char *)pDepth);

            float u0 = (t.tex[0][0] * wait0[0] + t.tex[1][0] * wait1[0] + t.tex[2][0] * wait2[0]);
            float u1 = (t.tex[0][0] * wait0[1] + t.tex[1][0] * wait1[1] + t.tex[2][0] * wait2[1]);
            float u2 = (t.tex[0][0] * wait0[2] + t.tex[1][0] * wait1[2] + t.tex[2][0] * wait2[2]);
            float u3 = (t.tex[0][0] * wait0[3] + t.tex[1][0] * wait1[3] + t.tex[2][0] * wait2[3]);

            float v0 = (t.tex[0][1] * wait0[0] + t.tex[1][1] * wait1[0] + t.tex[2][1] * wait2[0]);
            float v1 = (t.tex[0][1] * wait0[1] + t.tex[1][1] * wait1[1] + t.tex[2][1] * wait2[1]);
            float v2 = (t.tex[0][1] * wait0[2] + t.tex[1][1] * wait1[2] + t.tex[2][1] * wait2[2]);
            float v3 = (t.tex[0][1] * wait0[3] + t.tex[1][1] * wait1[3] + t.tex[2][1] * wait2[3]);

            u0 = SDL_clamp(u0, 0.0f, 1.0f);
            u1 = SDL_clamp(u1, 0.0f, 1.0f);
            u2 = SDL_clamp(u2, 0.0f, 1.0f);
            u3 = SDL_clamp(u3, 0.0f, 1.0f);

            v0 = SDL_clamp(v0, 0.0f, 1.0f);
            v1 = SDL_clamp(v1, 0.0f, 1.0f);
            v2 = SDL_clamp(v2, 0.0f, 1.0f);
            v3 = SDL_clamp(v3, 0.0f, 1.0f);

            u0 *= (float)texw - 1;
            u1 *= (float)texw - 1;
            u2 *= (float)texw - 1;
            u3 *= (float)texw - 1;

            v0 *= (float)texh - 1;
            v1 *= (float)texh - 1;
            v2 *= (float)texh - 1;
            v3 *= (float)texh - 1;

            unsigned char *texcolour0 = texdata + (((int)u0 + texw * (int)v0) * texbpp);
            unsigned char *texcolour1 = texdata + (((int)u1 + texw * (int)v1) * texbpp);
            unsigned char *texcolour2 = texdata + (((int)u2 + texw * (int)v2) * texbpp);
            unsigned char *texcolour3 = texdata + (((int)u3 + texw * (int)v3) * texbpp);

            __m128i final_colour = _mm_setr_epi8(texcolour0[2], texcolour0[1], texcolour0[0], -1,
                                                 texcolour1[2], texcolour1[1], texcolour1[0], -1,
                                                 texcolour2[2], texcolour2[1], texcolour2[0], -1,
                                                 texcolour3[2], texcolour3[1], texcolour3[0], -1);

            uint32_t *pixel_location = &rend.pixels[pixel_index];

#if 1 /* Fabian method */
            const __m128i original_pixel_data = _mm_load_si128((__m128i *)pixel_location);

            const __m128i write_mask    = _mm_castps_si128(sseWriteMask);
            const __m128i masked_output = _mm_or_si128(_mm_and_si128(write_mask, final_colour),
                                                       _mm_andnot_si128(write_mask, original_pixel_data));

            _mm_store_si128((__m128i *)pixel_location, masked_output);
#else
            // Mask-store 4-sample fragment values
            _mm_maskstore_epi32(
                (int *)pixel_location,
                _mm_castps_si128(EdgeFuncTestResult),
                final_colour);
#endif
        }
    }
}

static void draw_onexit(void) {} /* Do nothing */

static void draw_object(struct arena *arena)
{
    // #pragma omp parallel
    {
        obj_t          obj      = raster_state.obj;
        float         *pPos     = obj.pos;
        float         *pTex     = obj.texs;
        vertindices_t *pIndices = obj.indices;

        // #pragma omp for
        for (size_t i = 0; i < obj.num_f_rows; ++i)
        {
            struct triangle t = {0};

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
                t.pos[storeIndex][0] = pPos[3 * posIndex + 0];
                t.pos[storeIndex][1] = pPos[3 * posIndex + 1];
                t.pos[storeIndex][2] = pPos[3 * posIndex + 2];

                t.tex[storeIndex][0] = pTex[2 * texIndex + 0];
                t.tex[storeIndex][1] = pTex[2 * texIndex + 1];
            }
            draw_triangle(t);
        }
    }
}

#endif // __STEPPING_H__