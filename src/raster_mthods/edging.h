#ifndef __EDGING_H__
#define __EDGING_H__

#include "common.h"

static void drawtriangle(const triangle t, const vec3 texcoords[3], const mat4 MVP)
{
    vec3 screenspace[3] = {0};
    vec4 clipspace[3]   = {0};
    vec3 ndc[3]         = {0};

    // convert to clip space
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

    vec3 col1, col2, col3;
    v3div((vec3){1.0f, 0.0f, 0.0f}, clipspace[0][3], col1);
    v3div((vec3){0.0f, 1.0f, 0.0f}, clipspace[1][3], col2);
    v3div((vec3){0.0f, 0.0f, 1.0f}, clipspace[2][3], col3);

    vec3 uv[3];
    v3div(texcoords[0], clipspace[0][3], uv[0]);
    v3div(texcoords[1], clipspace[1][3], uv[1]);
    v3div(texcoords[2], clipspace[2][3], uv[2]);

    // Rasterize
    triangle tri = {
        {screenspace[0][0], screenspace[0][1], 1.0f / clipspace[0][3]},
        {screenspace[1][0], screenspace[1][1], 1.0f / clipspace[1][3]},
        {screenspace[2][0], screenspace[2][1], 1.0f / clipspace[2][3]},
    };

    float area     = edgefunc(tri[0], tri[1], tri[2]);
    float inv_area = 1.0f / area;

    for (int y = minY; y <= maxY; ++y)
    {
        for (int x = minX; x <= maxX; ++x)
        {
            vec3 point = {0.5f + (float)x, 0.5f + (float)y, 0.0f};

            float w0 = edgefunc(tri[1], tri[2], point);
            float w1 = edgefunc(tri[2], tri[0], point);
            float w2 = edgefunc(tri[0], tri[1], point);

            if (w0 < 0.0f || w1 < 0.0f || w2 < 0.0f)
                continue;

            w0 *= inv_area;
            w1 *= inv_area;
            w2 *= inv_area;

            // correction factor
            float cf = 1.0f / (w0 * tri[0][2] + w1 * tri[1][2] + w2 * tri[2][2]);

            const int   index = (x * GRAFIKA_SCREEN_WIDTH) + y;
            const float newZ  = w0 * screenspace[0][2] + w1 * screenspace[1][2] + w2 * screenspace[2][2];
            const float oldZ  = zbuffer[index];
            if (newZ < oldZ)
                continue;

            zbuffer[index] = newZ;

            float u = (uv[0][0] * w0 + uv[1][0] * w1 + uv[2][0] * w2) * cf;
            float v = (uv[0][1] * w0 + uv[1][1] * w1 + uv[2][1] * w2) * cf;

            u = fabsf(u);
            v = fabsf(v);

            u *= (float)state.tex.w - 1;
            v *= (float)state.tex.h - 1;

            unsigned char *texcolour   = state.tex.data + (((int)u + state.tex.w * (int)v) * state.tex.bpp);
            uint32_t       pixelcolour = (0xFF << 24) + (texcolour[2] << 16) + (texcolour[1] << 8) + (texcolour[0] << 0);

            // Draw pixel
            // 0xFF 00 00 00
            // uint32_t pixelcolour = (0xFF << 24) + (red << 16) + (blu << 8) + (gre << 0);
            grafika_setpixel(x, y, pixelcolour);
        }
    }
}

#endif // __EDGING_H__