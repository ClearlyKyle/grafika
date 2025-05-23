#include "common.h"

struct rasterstate raster_state = {0};

static void draw_triangle(const struct triangle t)
{
    vec3 screen_space[3] = {0};
    vec4 clip_space[3]   = {0};

    // convert to clip space
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
    vec3 ndc[3] = {0}, w_vals = {0};
    for (size_t i = 0; i < 3; i++)
    {
        w_vals[i] = 1.0f / clip_space[i][3]; // 1.0f / w
        ndc[i][0] = clip_space[i][0] * w_vals[i];
        ndc[i][1] = clip_space[i][1] * w_vals[i];
        ndc[i][2] = clip_space[i][2] * w_vals[i];
    }

    // back face culling (surface normal)
    vec3 sub10 = {0}, sub20 = {0}, normal = {0};
    v3_sub(ndc[1], ndc[0], sub10);
    v3_sub(ndc[2], ndc[0], sub20);
    v3_cross(sub10, sub20, normal);
    if (normal[2] > 0.0f)
        return;

    for (int i = 0; i < 3; ++i)
    {
        screen_space[i][0] = (ndc[i][0] + 1.0f) * 0.5f * GRAFIKA_SCREEN_WIDTH;
        screen_space[i][1] = (ndc[i][1] + 1.0f) * 0.5f * GRAFIKA_SCREEN_HEIGHT;
        screen_space[i][2] = (ndc[i][2] + 1.0f) * 0.5f;
    }

    // calculate bounding rectangle
    int AABB[4] = {0};
    AABB_make(screen_space, AABB);

    float area     = v3_edgefunc(screen_space[0], screen_space[1], screen_space[2]);
    float inv_area = 1.0f / area;

    // pre fetch tex data
    const int      tex_w    = raster_state.tex.w;
    const int      tex_h    = raster_state.tex.h;
    const int      tex_bpp  = raster_state.tex.bpp;
    unsigned char *tex_data = raster_state.tex.data;

    float *pDepthBuffer = rend.depth_buffer;

    // rasterize
    for (int y = AABB[1]; y <= AABB[3]; ++y)
    {
        for (int x = AABB[0]; x <= AABB[2]; ++x)
        {
            vec3 point = {0.5f + (float)x, 0.5f + (float)y, 0.0f};

            float w0 = v3_edgefunc(screen_space[1], screen_space[2], point);
            float w1 = v3_edgefunc(screen_space[2], screen_space[0], point);
            float w2 = v3_edgefunc(screen_space[0], screen_space[1], point);

            if (w0 < 0.0f || w1 < 0.0f || w2 < 0.0f)
                continue;

            w0 *= inv_area;
            w1 *= inv_area;
            w2 *= inv_area;

            const int index = (x * GRAFIKA_SCREEN_WIDTH) + y;

            const float depth = w0 * screen_space[0][2] + w1 * screen_space[1][2] + w2 * screen_space[2][2];

            float      *oldZ = pDepthBuffer + index;
            const float newZ = depth;

            if (newZ > *oldZ)
                continue;

            *oldZ = newZ;

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

            uint32_t pixelcolour = ((uint32_t)tex_colour[0] << 16) |
                                   ((uint32_t)tex_colour[1] << 8) |
                                   (uint32_t)tex_colour[2];

            grafika_setpixel(x, y, pixelcolour);
        }
    }
}

static struct triangle *triangles = NULL;

void draw_object(struct arena *arena)
{
    obj_t          obj      = raster_state.obj;
    float         *pPos     = obj.pos;
    float         *pTex     = obj.texs;
    vertindices_t *pIndices = obj.indices;

#if 1
    if (triangles == NULL)
        triangles = arena_alloc_aligned(arena, sizeof(struct triangle) * obj.num_f_rows, 16);

    // #pragma omp parallel
    {
        // #pragma omp for
        //  gather triangles
        for (size_t i = 0; i < obj.num_f_rows; ++i)
        { // CW triangles (which blender uses)
            vertindices_t indices0   = pIndices[i * 3 + 2];
            const int     pos_index0 = indices0.v_idx;

            triangles[i].pos[0][0] = pPos[3 * pos_index0 + 0];
            triangles[i].pos[0][1] = pPos[3 * pos_index0 + 1];
            triangles[i].pos[0][2] = pPos[3 * pos_index0 + 2];

            const int tex_index0   = indices0.vt_idx;
            triangles[i].tex[0][0] = pTex[2 * tex_index0 + 0];
            triangles[i].tex[0][1] = pTex[2 * tex_index0 + 1];

            vertindices_t indices1   = pIndices[i * 3 + 1];
            const int     pos_index1 = indices1.v_idx;

            triangles[i].pos[1][0] = pPos[3 * pos_index1 + 0];
            triangles[i].pos[1][1] = pPos[3 * pos_index1 + 1];
            triangles[i].pos[1][2] = pPos[3 * pos_index1 + 2];

            const int tex_index1   = indices1.vt_idx;
            triangles[i].tex[1][0] = pTex[2 * tex_index1 + 0];
            triangles[i].tex[1][1] = pTex[2 * tex_index1 + 1];

            vertindices_t indices2   = pIndices[i * 3 + 0];
            const int     pos_index2 = indices2.v_idx;

            triangles[i].pos[2][0] = pPos[3 * pos_index2 + 0];
            triangles[i].pos[2][1] = pPos[3 * pos_index2 + 1];
            triangles[i].pos[2][2] = pPos[3 * pos_index2 + 2];

            const int tex_index2   = indices2.vt_idx;
            triangles[i].tex[2][0] = pTex[2 * tex_index2 + 0];
            triangles[i].tex[2][1] = pTex[2 * tex_index2 + 1];
        }

        // #pragma omp for
        for (size_t i = 0; i < obj.num_f_rows; ++i)
        {
            draw_triangle(triangles[i]);
        }
    }
#else
    #pragma omp parallel
    {
    #pragma omp for
        for (size_t i = 0; i < obj.num_f_rows; ++i)
        {
            struct triangle t;

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
#endif
}