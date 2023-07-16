#ifndef __STEPPING_H__
#define __STEPPING_H__

#include "common.h"

static void draw_triangle(const triangle_t t)
{
    // convert to clip space
    vec4 clipspace[3] = {0};
    for (int i = 0; i < 3; ++i)
    {
        vec4 pos = {t.pos[i][0], t.pos[i][1], t.pos[i][2], 1.0f};
        m4_mul_v4(state.MVP, pos, clipspace[i]);

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
    vec3 ndc[3], w_vals;
    for (size_t i = 0; i < 3; i++)
    {
        w_vals[i] = 1.0f / clipspace[i][3]; // 1.0f / w
        ndc[i][0] = clipspace[i][0] * w_vals[i];
        ndc[i][1] = clipspace[i][1] * w_vals[i];
        ndc[i][2] = clipspace[i][2] * w_vals[i];
    }

    // back face culling (surface normal)
    vec3 sub10, sub20, normal;
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
    int AABB[4];
    AABB_make(screenspace, AABB);

    const float dY0 = screenspace[2][1] - screenspace[1][1], dX0 = screenspace[1][0] - screenspace[2][0];
    const float dY1 = screenspace[0][1] - screenspace[2][1], dX1 = screenspace[2][0] - screenspace[0][0];
    const float dY2 = screenspace[1][1] - screenspace[0][1], dX2 = screenspace[0][0] - screenspace[1][0];

    const float C0 = (screenspace[2][0] * screenspace[1][1]) - (screenspace[2][1] * screenspace[1][0]);
    const float C1 = (screenspace[0][0] * screenspace[2][1]) - (screenspace[0][1] * screenspace[2][0]);
    const float C2 = (screenspace[1][0] * screenspace[0][1]) - (screenspace[1][1] * screenspace[0][0]);

    const float inv_area = 1.0f / (dX1 * dY2 - dY1 * dX2);

    /* Should be ints? */
    // const vec3 P     = {minX + 0.5f, minY + 0.5f, 0.0f};
    float alpha = (dY0 * ((float)AABB[0] + 0.5f)) + (dX0 * ((float)AABB[1] + 0.5f)) + C0;
    float betaa = (dY1 * ((float)AABB[0] + 0.5f)) + (dX1 * ((float)AABB[1] + 0.5f)) + C1;
    float gamma = (dY2 * ((float)AABB[0] + 0.5f)) + (dX2 * ((float)AABB[1] + 0.5f)) + C2;

    screenspace[1][2] = (screenspace[1][2] - screenspace[0][2]) * inv_area;
    screenspace[2][2] = (screenspace[2][2] - screenspace[0][2]) * inv_area;

    const float zstep = dY1 * screenspace[1][2] + dY2 * screenspace[2][2];

    // pre fetch tex data
    const int      texw    = state.tex.w;
    const int      texh    = state.tex.h;
    const int      texbpp  = state.tex.bpp;
    unsigned char *texdata = state.tex.data;

    float *pDepthBuffer = rend.depth_buffer;

    // rasterize
    for (int y = AABB[1]; y <= AABB[3]; ++y)
    {
        // Barycentric coordinates at start of row
        float w0 = alpha;
        float w1 = betaa;
        float w2 = gamma;

        float depth = screenspace[0][2] + (screenspace[1][2] * betaa) + (screenspace[2][2] * gamma);

        for (int x = AABB[0]; x <= AABB[2]; ++x,            //
                                            w0 += dY0,      //
                                            w1 += dY1,      //
                                            w2 += dY2,      // One step to the right
                                            depth += zstep) // Step the depth
        {
            if (w0 < 0.0f || w1 < 0.0f || w2 < 0.0f)
                continue;

            const int index = (x * GRAFIKA_SCREEN_WIDTH) + y;
            // const float invZ  = 1.0f / depth;
            const float invZ = depth;

            float *oldZ = pDepthBuffer + index;

            if (invZ > *oldZ)
                continue;

            *oldZ = invZ;

            // Weights, see what i did there ;)
            const float wait0 = w0 * w_vals[0];
            const float wait1 = w1 * w_vals[1];
            const float wait2 = w2 * w_vals[2];

            // correction factor
            const float cf = 1.0f / (wait0 + wait1 + wait2);

            float u = (t.tex[0][0] * wait0 + t.tex[1][0] * wait1 + t.tex[2][0] * wait2) * cf;
            float v = (t.tex[0][1] * wait0 + t.tex[1][1] * wait1 + t.tex[2][1] * wait2) * cf;

            // u = fabsf(u);
            // v = fabsf(v);

            u *= (float)texw - 1;
            v *= (float)texh - 1;

            unsigned char *texcolour   = texdata + (((int)u + texw * (int)v) * texbpp);
            uint32_t       pixelcolour = (0xFF << 24) + (texcolour[2] << 16) + (texcolour[1] << 8) + (texcolour[0] << 0);

            grafika_setpixel(x, y, pixelcolour);
        }
        // One row step
        alpha += dX0;
        betaa += dX1;
        gamma += dX2;
    }
}

static void draw_onexit(void) {} /* Do nothing */

static void draw_object(void)
{
#pragma omp parallel
    {
        obj_t          obj      = state.obj;
        float         *pPos     = obj.pos;
        float         *pTex     = obj.texs;
        vertindices_t *pIndices = obj.indices;

#pragma omp for
        for (size_t i = 0; i < obj.num_f_rows; ++i)
        {
            triangle_t t = {0};

            for (size_t j = 0; j < 3; ++j)
            {
                vertindices_t indices = pIndices[i * 3 + j];

                const int posIndex = indices.v_idx;
                const int texIndex = indices.vt_idx;

/* send triangles to 'draw_triangle' in CCW order, set the winding mode of the loaded model */
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