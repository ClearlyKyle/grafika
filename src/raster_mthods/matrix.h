#ifndef __MATRIX_H__
#define __MATRIX_H__

#include "common.h"

/*
Use this to precompute the tie braker edge conditions, use the E value before
looping through the pixels
*/
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

static inline float Interpolate_Values(vec3 littlef_values, vec3 attribute)
{
    return v3dot(littlef_values, attribute);
}

static void drawtriangle(triangle t, vec3 texcoords[3], mat4 MVP)
{
    // convert to clip space
    vec4 clipspace[3] = {0};
    for (int i = 0; i < 3; ++i)
    {
        //  clip space position
        m4mulv4((const float(*)[4])MVP, (vec4){t[i][0], t[i][1], t[i][2], 1.0f}, clipspace[i]);

        //// clipping (is this correct?)
        // const float x = fabsf(clipspace[i][0]);
        // const float y = fabsf(clipspace[i][1]);
        // const float w = fabsf(clipspace[i][3]);

        // if ((-w <= x && x <= w) || (-w <= y && y <= w))
        //     continue;
        // else
        //     return;
    }

    vec4        screenspace[3] = {0};
    const float half_width     = 0.5f * GRAFIKA_SCREEN_WIDTH;
    const float half_heigh     = 0.5f * GRAFIKA_SCREEN_HEIGHT;
    for (int i = 0; i < 3; ++i)
    {
        screenspace[i][0] = clipspace[i][0] * half_width + clipspace[i][3] * half_width;
        screenspace[i][1] = clipspace[i][1] * half_heigh + clipspace[i][3] * half_heigh;
        screenspace[i][2] = clipspace[i][2];
        screenspace[i][3] = clipspace[i][3];
    }

    /* From the paper, calculate out "A" matrix */
    const float a0 = (screenspace[1][1] * screenspace[2][3]) - (screenspace[2][1] * screenspace[1][3]);
    const float a1 = (screenspace[2][1] * screenspace[0][3]) - (screenspace[0][1] * screenspace[2][3]);
    const float a2 = (screenspace[0][1] * screenspace[1][3]) - (screenspace[1][1] * screenspace[0][3]);

    const float b0 = (screenspace[2][0] * screenspace[1][3]) - (screenspace[1][0] * screenspace[2][3]);
    const float b1 = (screenspace[0][0] * screenspace[2][3]) - (screenspace[2][0] * screenspace[0][3]);
    const float b2 = (screenspace[1][0] * screenspace[0][3]) - (screenspace[0][0] * screenspace[1][3]);

    const float c0 = (screenspace[1][0] * screenspace[2][1]) - (screenspace[2][0] * screenspace[1][1]);
    const float c1 = (screenspace[2][0] * screenspace[0][1]) - (screenspace[0][0] * screenspace[2][1]);
    const float c2 = (screenspace[0][0] * screenspace[1][1]) - (screenspace[1][0] * screenspace[0][1]);

    const float detM = (c0 * screenspace[0][3]) + (c1 * screenspace[1][3]) + (c2 * screenspace[2][3]);
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
        proj[i][0] = screenspace[i][0] / screenspace[i][3];
        proj[i][1] = screenspace[i][1] / screenspace[i][3];
        proj[i][2] = screenspace[i][2] / screenspace[i][3];
    }

    /* get the bounding box of the triangle */
    float fminX = fminf(proj[0][0], fminf(proj[1][0], proj[2][0]));
    float fminY = fminf(proj[0][1], fminf(proj[1][1], proj[2][1]));
    float fmaxX = fmaxf(proj[0][0], fmaxf(proj[1][0], proj[2][0]));
    float fmaxY = fmaxf(proj[0][1], fmaxf(proj[1][1], proj[2][1]));

    int minX = max(0, min((int)floorf(fminX), GRAFIKA_SCREEN_WIDTH - 1));
    int minY = max(0, min((int)floorf(fminY), GRAFIKA_SCREEN_HEIGHT - 1));
    int maxX = max(0, min((int)floorf(fmaxX), GRAFIKA_SCREEN_WIDTH - 1));
    int maxY = max(0, min((int)floorf(fmaxY), GRAFIKA_SCREEN_HEIGHT - 1));

    // Evaluaate edge equation at first tile origin
    const float edgeFunc0 = (E0[0] * (float)minX) + (E0[1] * (float)minY) + E0[2];
    const float edgeFunc1 = (E1[0] * (float)minX) + (E1[1] * (float)minY) + E1[2];
    const float edgeFunc2 = (E2[0] * (float)minX) + (E2[1] * (float)minY) + E2[2];

    /* Pre compute some tie breaker test results */
    const bool pre_comp_tie_E0 = Tie_Breaker_AB_Test(E0);
    const bool pre_comp_tie_E1 = Tie_Breaker_AB_Test(E1);
    const bool pre_comp_tie_E2 = Tie_Breaker_AB_Test(E2);

    // Start rasterizing by looping over pixels to output a per-pixel color
    for (int y = minY, step_y = 0; y <= maxY; y++, step_y++)
    {
        for (int x = minX, step_x = 0; x <= maxX; x++, step_x++)
        {
            // Step from edge function by multiples of the edge values
            // edgevalue * step_multiple

            // Evaluate edge functions at current fragment
            // E(x + s, y + t) = E(x, y) + sa + tb,
            const float edgeFuncTR0 = edgeFunc0 + ((E0[0] * step_x) + (E0[1] * step_y));
            const float edgeFuncTR1 = edgeFunc1 + ((E1[0] * step_x) + (E1[1] * step_y));
            const float edgeFuncTR2 = edgeFunc2 + ((E2[0] * step_x) + (E2[1] * step_y));

            // Check if the current point is inside a traingle using Tie breaker rules
            const bool TRForEdge0 = Edge_Tie_Breaker(edgeFuncTR0, pre_comp_tie_E0);
            if (TRForEdge0)
                continue;

            const bool TRForEdge1 = Edge_Tie_Breaker(edgeFuncTR1, pre_comp_tie_E1);
            if (TRForEdge1)
                continue;

            const bool TRForEdge2 = Edge_Tie_Breaker(edgeFuncTR2, pre_comp_tie_E2);
            if (TRForEdge2)
                continue;

            /* Calculate Interpolation Coefficients */
            const float F0 = edgeFuncTR0;
            const float F1 = edgeFuncTR1;
            const float F2 = edgeFuncTR2;

            const float r = 1.0f / (F0 + F1 + F2);

            vec3 littlef_values = {edgeFuncTR0 * r, edgeFuncTR1 * r, edgeFuncTR2 * r};

            /* Interpolate Depth value */
            const int index = (x * GRAFIKA_SCREEN_WIDTH) + y;

            const float depth = Interpolate_Values(littlef_values, (vec3){screenspace[0][2], screenspace[1][2], screenspace[2][2]});

            // const float invZ = 1.0f / depth;
            const float invZ = depth;

            float *oldZ = &zbuffer[index];

            // Perform depth test
            if (invZ > *oldZ)
                continue;

            // Depth test passed, update depth buffer value
            *oldZ = invZ;

            // Interpolate texture coordinates
            float u = Interpolate_Values(littlef_values, (vec4){texcoords[0][0], texcoords[1][0], texcoords[2][0]});
            float v = Interpolate_Values(littlef_values, (vec4){texcoords[0][1], texcoords[1][1], texcoords[2][1]});

            u *= (float)state.tex.w - 1;
            v *= (float)state.tex.h - 1;

            unsigned char *texcolour   = state.tex.data + (((int)u + state.tex.w * (int)v) * state.tex.bpp);
            uint32_t       pixelcolour = (0xFF << 24) + (texcolour[2] << 16) + (texcolour[1] << 8) + (texcolour[0] << 0);

            grafika_setpixel(x, y, pixelcolour);
        }
    }
}

#endif // __MATRIX_H__