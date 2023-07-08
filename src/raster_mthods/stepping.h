#ifndef __STEPPING_H__
#define __STEPPING_H__

#include "../grafika.h"
#include "../tex.h"
#include "../obj.h"
#include "../timer.h"
#include "../matematika.h"

typedef vec3 triangle[3];

typedef struct rasterstate
{
    mat4 proj;
    mat4 view;
    mat4 model;

    obj_t obj;
    tex_t tex;
} rasterstate_t;

static rasterstate_t state = {0};

// NOTE : Should this be in raster state?
// probably better to move it into a "common.h" file
static float zbuffer[GRAFIKA_SCREEN_HEIGHT * GRAFIKA_SCREEN_WIDTH];

// static inline float edgefunc(const vec3 a, const vec3 b, const vec3 c)
//{
//     return (c[0] - a[0]) * (b[1] - a[1]) - (c[1] - a[1]) * (b[0] - a[0]);
// }

static inline float edgefunc(const vec3 A, const vec3 B, const vec3 P)
{
    return (B[0] - A[0]) * (P[1] - A[1]) - (B[1] - A[1]) * (P[0] - A[0]);
}

void drawtriangle(triangle t, vec3 texcoords[3], mat4 MVP)
{
    // convert to clip space
    vec4 clipspace[3] = {0};
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
    vec3 ndc[3] = {0};
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

    vec3 screenspace[3] = {0};
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

    vec3 uv[3];
    v3div(texcoords[0], clipspace[0][3], uv[0]);
    v3div(texcoords[1], clipspace[1][3], uv[1]);
    v3div(texcoords[2], clipspace[2][3], uv[2]);

    float A01 = screenspace[2][1] - screenspace[1][1], B01 = screenspace[1][0] - screenspace[2][0];
    float A12 = screenspace[0][1] - screenspace[2][1], B12 = screenspace[2][0] - screenspace[0][0];
    float A20 = screenspace[1][1] - screenspace[0][1], B20 = screenspace[0][0] - screenspace[1][0];

    float C01 = (screenspace[2][0] * screenspace[1][1]) - (screenspace[2][1] * screenspace[1][0]);
    float C12 = (screenspace[0][0] * screenspace[2][1]) - (screenspace[0][1] * screenspace[2][0]);
    float C20 = (screenspace[1][0] * screenspace[0][1]) - (screenspace[1][1] * screenspace[0][0]);

    const float inv_area = 1.0f / (A20 * B12 - A12 * B20);

    /* Should be ints? */
    // const vec3 P     = {minX + 0.5f, minY + 0.5f, 0.0f};
    float alpha = (A01 * ((float)minX + 0.5f)) + (B01 * ((float)minY + 0.5f)) + C01;
    float beta  = (A12 * ((float)minX + 0.5f)) + (B12 * ((float)minY + 0.5f)) + C12;
    float gamma = (A20 * ((float)minX + 0.5f)) + (B20 * ((float)minY + 0.5f)) + C20;

    // Rasterize
    triangle tri = {
        {screenspace[0][0], screenspace[0][1], 1.0f / clipspace[0][3]},
        {screenspace[1][0], screenspace[1][1], 1.0f / clipspace[1][3]},
        {screenspace[2][0], screenspace[2][1], 1.0f / clipspace[2][3]},
    };

    for (int y = minY; y <= maxY; ++y)
    {
        // Barycentric coordinates at start of row
        float w0 = alpha;
        float w1 = beta;
        float w2 = gamma;

        for (int x = minX; x <= maxX; ++x,
                 w0 += A01,
                 w1 += A12,
                 w2 += A20) // One step to the right
        {
            if (w0 < 0.0f || w1 < 0.0f || w2 < 0.0f)
                continue;

            const float bw0 = w0 * inv_area;
            const float bw1 = w1 * inv_area;
            const float bw2 = w2 * inv_area;

            const int   index = (x * GRAFIKA_SCREEN_WIDTH) + y;
            const float newZ  = bw0 * screenspace[0][2] + bw1 * screenspace[1][2] + bw2 * screenspace[2][2];
            const float invZ  = 1.0f / newZ;

            float *oldZ = &zbuffer[index];

            if (invZ > *oldZ)
                continue;

            *oldZ = invZ;

            // Weights, see what i did there ;)
            const float wait0 = bw0 * tri[0][2];
            const float wait1 = bw1 * tri[1][2];
            const float wait2 = bw2 * tri[2][2];

            // correction factor
            const float cf = 1.0f / (wait0 + wait1 + wait2);

            float u = (texcoords[0][0] * wait0 + texcoords[1][0] * wait1 + texcoords[2][0] * wait2) * cf;
            float v = (texcoords[0][1] * wait0 + texcoords[1][1] * wait1 + texcoords[2][1] * wait2) * cf;

            u = fabsf(u);
            v = fabsf(v);

            u *= (float)state.tex.w - 1;
            v *= (float)state.tex.h - 1;

            unsigned char *texcolour   = state.tex.data + (((int)u + state.tex.w * (int)v) * state.tex.bpp);
            uint32_t       pixelcolour = (0xFF << 24) + (texcolour[2] << 16) + (texcolour[1] << 8) + (texcolour[0] << 0);

            grafika_setpixel(x, y, pixelcolour);
        }
        // One row step
        alpha += B01;
        beta += B12;
        gamma += B20;
    }
}

#endif // __STEPPING_H__