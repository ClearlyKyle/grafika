#ifndef __PHONG_H__
#define __PHONG_H__

#include "common.h"

#define PHONG_AMBI_AMOUNT 0.2f   /* 0.0f -> 1.0f */
#define PHONG_SPEC_AMOUNT 0.4f   /* 0.0f -> 1.0f */
#define PHONG_SHININESS   128.0f /* pow of 2 pls */
#define LIGHT_POSITION \
    (vec3){            \
        0.5f, 0.5f, 5.0f}

typedef struct phong
{
    mat3 nrm_matrix;
} phong_t;

static phong_t phong_data           = {0};
static float  *transformed_vertices = NULL;

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
    if (edge_value > 0.0f)
        return true;
    if (edge_value < 0.0f)
        return false;

    return AB_test_result;
}

static inline float Interpolate_Values(vec3 littlef_values, vec3 attribute)
{
    return v3_dot(littlef_values, attribute);
}

static inline float calculate_specular_amount(vec3 L, vec3 E, vec3 N, const float shininess)
{
#if 1
    {
        // BLINN PHONG
        vec3 halfway_direction = {0};
        v3_add(L, E, halfway_direction);
        v3_norm(halfway_direction);

        const float specular_intensity = (const float)pow(fmax(v3_dot(N, halfway_direction), 0.0), shininess);
        return specular_intensity;
    }
#else
    {
        // Regular PHONG
        vec3 R = {0};
        v3_reflect(L, N, R);
        R[0] = -R[0];
        R[1] = -R[1];
        R[2] = -R[2];

        const float specular_intensity = powf(fmaxf(v3_dot(E, R), 0.0f), shininess);
        return specular_intensity;
    }
#endif
}

static inline void phong_lighting(vec3 diff_colour, vec3 frag_pos, vec3 N, vec3 FragColour)
{
    //  Normalise the Noraml - N
    v3_norm(N);

    // L - direction to the light source
    vec3 L = {0};
    v3_sub(LIGHT_POSITION, frag_pos, L);
    v3_norm(L);

    // E - view direction (x and y position is fixed)
    vec3 E = {0.0f, 0.0f, 1.0f};

    // Ambient Term:
    vec3  Iamb           = {0};
    float ambient_amount = PHONG_AMBI_AMOUNT;
    v3_scale(diff_colour, ambient_amount, Iamb);

    // Diffuse Term:
    vec3        Idiff          = {0};
    const float dot_product    = v3_dot(L, N);
    const float diffuse_amount = fmaxf(dot_product, 0.0f);
    v3_scale(diff_colour, diffuse_amount, Idiff);

    // Specular Term:
    vec3        Ispec           = {0};
    const float specular_amount = calculate_specular_amount(L, E, N, PHONG_SHININESS);
    v3_broadcast(Ispec, specular_amount * PHONG_SPEC_AMOUNT);

    // float shading_amount = ambient_amount + diffuse_amount + specular_amount;
    // shading_amount       = shading_amount > 1.0f ? 1.0f : shading_amount; // clamp to 1.0f max

    v3_add(Iamb, Idiff, FragColour);
    v3_add(FragColour, Ispec, FragColour);

    v3_clamp(FragColour, 0.0f, 1.0f);
}

static void draw_triangle(vec4 trans[3], vec3 raw[3], vec3 nrm[3], vec2 texcoord[3])
{
    /* From the paper, calculate out "A" matrix */
    const float a0 = (trans[1][1] * trans[2][3]) - (trans[2][1] * trans[1][3]);
    const float a1 = (trans[2][1] * trans[0][3]) - (trans[0][1] * trans[2][3]);
    const float a2 = (trans[0][1] * trans[1][3]) - (trans[1][1] * trans[0][3]);

    const float b0 = (trans[2][0] * trans[1][3]) - (trans[1][0] * trans[2][3]);
    const float b1 = (trans[0][0] * trans[2][3]) - (trans[2][0] * trans[0][3]);
    const float b2 = (trans[1][0] * trans[0][3]) - (trans[0][0] * trans[1][3]);

    const float c0 = (trans[1][0] * trans[2][1]) - (trans[2][0] * trans[1][1]);
    const float c1 = (trans[2][0] * trans[0][1]) - (trans[0][0] * trans[2][1]);
    const float c2 = (trans[0][0] * trans[1][1]) - (trans[1][0] * trans[0][1]);

    const float detM = (c0 * trans[0][3]) + (c1 * trans[1][3]) + (c2 * trans[2][3]);
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
        proj[i][0] = trans[i][0] / trans[i][3];
        proj[i][1] = trans[i][1] / trans[i][3];
        proj[i][2] = trans[i][2] / trans[i][3];
    }

    /* transform our normals using the "normals matrix" this is needed due to
        the model matrix having the ability to scale and scew the model, and our
        normal do not like this */
    // NOTE: This might need to be a m3 * v3...
    vec3 new_nrm[3] = {0};
    m3_mul_v3(phong_data.nrm_matrix, nrm[0], new_nrm[0]);
    m3_mul_v3(phong_data.nrm_matrix, nrm[1], new_nrm[1]);
    m3_mul_v3(phong_data.nrm_matrix, nrm[2], new_nrm[2]);

    // get world space vertex positions
    vec4 ws[3] = {0};
    m4_mul_v4(raster_state.model, (vec4){raw[0][0], raw[0][1], raw[0][2], 1.0f}, ws[0]);
    m4_mul_v4(raster_state.model, (vec4){raw[1][0], raw[1][1], raw[1][2], 1.0f}, ws[1]);
    m4_mul_v4(raster_state.model, (vec4){raw[2][0], raw[2][1], raw[2][2], 1.0f}, ws[2]);

    // calculate bounding rectangle
    int AABB[4] = {0};
    AABB_make(proj, AABB);

    // Evaluaate edge equation at first tile origin
    const float edgeFunc0 = (E0[0] * (float)AABB[0]) + (E0[1] * (float)AABB[1]) + E0[2];
    const float edgeFunc1 = (E1[0] * (float)AABB[0]) + (E1[1] * (float)AABB[1]) + E1[2];
    const float edgeFunc2 = (E2[0] * (float)AABB[0]) + (E2[1] * (float)AABB[1]) + E2[2];

    /* Pre compute some tie breaker test results */
    const bool pre_comp_tie_E0 = Tie_Breaker_AB_Test(E0);
    const bool pre_comp_tie_E1 = Tie_Breaker_AB_Test(E1);
    const bool pre_comp_tie_E2 = Tie_Breaker_AB_Test(E2);

    // pre fetch tex data
    const int      texw    = raster_state.tex.w;
    const int      texh    = raster_state.tex.h;
    const int      texbpp  = raster_state.tex.bpp;
    unsigned char *texdata = raster_state.tex.data;

    float *pDepthBuffer = rend.depth_buffer;

    // Start rasterizing by looping over pixels to output a per-pixel color
    for (int y = AABB[1], step_y = 0; y <= AABB[3]; y++, step_y++)
    {
        for (int x = AABB[0], step_x = 0; x <= AABB[2]; x++, step_x++)
        {
            // Step from edge function by multiples of the edge values
            // edgevalue * step_multiple

            // Evaluate edge functions at current fragment
            // E(x + s, y + t) = E(x, y) + sa + tb,
            const float edgeFuncTR0 = edgeFunc0 + ((E0[0] * (float)step_x) + (E0[1] * (float)step_y));
            const float edgeFuncTR1 = edgeFunc1 + ((E1[0] * (float)step_x) + (E1[1] * (float)step_y));
            const float edgeFuncTR2 = edgeFunc2 + ((E2[0] * (float)step_x) + (E2[1] * (float)step_y));

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

            const float depth = Interpolate_Values(littlef_values, (vec3){trans[0][2], trans[1][2], trans[2][2]});

            float *oldZ = pDepthBuffer + index;
            // const float invZ = 1.0f / depth;
            const float invZ = depth;

            // Perform depth test
            if (invZ > *oldZ)
                continue;

            // Depth test passed, update depth buffer value
            *oldZ = invZ;

            vec3 frag_pos;
            frag_pos[0] = Interpolate_Values(littlef_values, (vec3){ws[0][0], ws[1][0], ws[2][0]}); // x
            frag_pos[1] = Interpolate_Values(littlef_values, (vec3){ws[0][1], ws[1][1], ws[2][1]}); // y
            frag_pos[2] = Interpolate_Values(littlef_values, (vec3){ws[0][2], ws[1][2], ws[2][2]}); // z

            vec3 frag_nrm;
            frag_nrm[0] = Interpolate_Values(littlef_values, (vec3){new_nrm[0][0], new_nrm[1][0], new_nrm[2][0]}); // x
            frag_nrm[1] = Interpolate_Values(littlef_values, (vec3){new_nrm[0][1], new_nrm[1][1], new_nrm[2][1]}); // y
            frag_nrm[2] = Interpolate_Values(littlef_values, (vec3){new_nrm[0][2], new_nrm[1][2], new_nrm[2][2]}); // z

            // Interpolate texture coordinates
            float u = Interpolate_Values(littlef_values, (vec3){texcoord[0][0], texcoord[1][0], texcoord[2][0]});
            float v = Interpolate_Values(littlef_values, (vec3){texcoord[0][1], texcoord[1][1], texcoord[2][1]});

            u *= (float)texw - 1;
            v *= (float)texh - 1;

            unsigned char *texcolour = texdata + (((int)u + texw * (int)v) * texbpp);

            vec3 diffuse_colour;
            diffuse_colour[0] = texcolour[0] / 255.0f;
            diffuse_colour[1] = texcolour[1] / 255.0f;
            diffuse_colour[2] = texcolour[2] / 255.0f;

            // diffuse_colour[0] = (float)(texcolour[0] >> 8);
            // diffuse_colour[1] = (float)(texcolour[1] >> 8);
            // diffuse_colour[2] = (float)(texcolour[2] >> 8);

            vec3 frag_colour = {0};
            phong_lighting(diffuse_colour, frag_pos, frag_nrm, frag_colour);

            unsigned char red = (unsigned char)(frag_colour[0] * 255.0f);
            unsigned char gre = (unsigned char)(frag_colour[1] * 255.0f);
            unsigned char blu = (unsigned char)(frag_colour[2] * 255.0f);

            // uint32_t pixelcolour = ((uint32_t)0xFFU << 24) | ((uint32_t)blu << 16) | ((uint32_t)gre << 8) | (uint32_t)red;
            uint32_t pixel_colour = ((uint32_t)red << 16) | ((uint32_t)gre << 8) | (uint32_t)blu;

            grafika_setpixel(x, y, pixel_colour);
        }
    }
}

static void draw_onexit(void) {} /* Do nothing */

static void draw_object(struct arena *arena)
{
    mat4 cum_matrix = {0};
    m4_mul_m4(VP_MATRIX, raster_state.MVP, cum_matrix);

    mat3 tmp = {0};
    m3_from_m4(raster_state.model, tmp);
    m3_inv(tmp, tmp);
    m3_transpose(tmp, phong_data.nrm_matrix);

    const obj_t object = raster_state.obj;

    if (!transformed_vertices)
        transformed_vertices = arena_alloc_aligned(arena, sizeof(vec4) * object.num_pos, 16);

#pragma omp parallel
    {
        float *pos = object.pos;

#pragma omp for
        // transform all verts to screen space
        for (size_t i = 0; i < object.num_pos; i++)
        {
            vec4 vertex = {pos[i * 3 + 0],
                           pos[i * 3 + 1],
                           pos[i * 3 + 2],
                           1.0f};
            m4_mul_v4(cum_matrix, vertex, &transformed_vertices[i * 4]);
        }

        float         *pPos     = transformed_vertices;
        float         *pRawpos  = object.pos;
        float         *pNrm     = object.norms;
        float         *pTex     = object.texs;
        vertindices_t *pIndices = object.indices;

#pragma omp for
        for (size_t i = 0; i < object.num_f_rows; ++i)
        {
            vec4 trans[3]     = {0};
            vec3 raw[3]       = {0};
            vec3 nrm[3]       = {0};
            vec2 texcoords[3] = {0};

            for (size_t j = 0; j < 3; ++j)
            {
                vertindices_t indices = pIndices[i * 3 + j];

                const int posIndex = indices.v_idx;
                const int nrmIndex = indices.vn_idx * 3;
                const int texIndex = indices.vt_idx * 2;
#if 0
                const size_t storeIndex = j; // CCW triangles
#else
                const size_t storeIndex = 2 - j; // CW triangles (which blender uses)
#endif
                float *pP            = pPos + 4 * posIndex;
                float *pRP           = pRawpos + 3 * posIndex;
                trans[storeIndex][0] = pP[0];
                trans[storeIndex][1] = pP[1];
                trans[storeIndex][2] = pP[2];
                trans[storeIndex][3] = pP[3];

                raw[storeIndex][0] = pRP[0];
                raw[storeIndex][1] = pRP[1];
                raw[storeIndex][2] = pRP[2];

                nrm[storeIndex][0] = pNrm[nrmIndex + 0];
                nrm[storeIndex][1] = pNrm[nrmIndex + 1];
                nrm[storeIndex][2] = pNrm[nrmIndex + 2];

                texcoords[storeIndex][0] = pTex[texIndex + 0];
                texcoords[storeIndex][1] = pTex[texIndex + 1];
            }
            draw_triangle(trans, raw, nrm, texcoords);
        }
    }
}

#endif // __PHONG_H__