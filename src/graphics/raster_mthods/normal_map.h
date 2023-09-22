#ifndef __NORMAL_MAP_H__
#define __NORMAL_MAP_H__

#include "common.h"

#define PHONG_AMBI_AMOUNT 0.2f  /* 0.0f -> 1.0f */
#define PHONG_SPEC_AMOUNT 0.4f  /* 0.0f -> 1.0f */
#define PHONG_SHININESS   64.0f /* pow of 2 pls */
#define LIGHT_POSITION   \
    (vec3)               \
    {                    \
        0.5f, 0.5f, 5.0f \
    }

typedef struct nrm_map
{
    mat3  nrm_matrix;
    mat3  TBN;
    tex_t nrm_tex;
} nrm_map_t;

static nrm_map_t normal_mapping_data  = {0};
static float    *transformed_vertices = NULL;

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
#if 0
    {
        // BLINN PHONG
        vec3 halfway_direction;
        v3_add(L, E, halfway_direction);
        v3_norm(halfway_direction);

        const float specular_intensity = (const float)pow(fmax(v3_dot(N, halfway_direction), 0.0), shininess);
        return specular_intensity;
    }
#else
    {
        // Regular PHONG
        vec3 R;
        v3_reflect(L, N, R);
        R[0] = -R[0];
        R[1] = -R[1];
        R[2] = -R[2];

        const float specular_intensity = powf(fmaxf(v3_dot(E, R), 0.0f), shininess);
        return specular_intensity;
    }
#endif
}

static inline void phong_lighting(vec3 diff_colour, vec3 frag_pos, vec3 cam_pos, vec3 light_pos, vec3 N, vec3 frag_colour)
{
    //  Normalise the Noraml - N
    v3_norm(N);

    // L - direction to the light source
    vec3 L = {0};
    v3_sub(light_pos, frag_pos, L);
    v3_norm(L);

    // E - view direction
    vec3 E;
    v3_sub(cam_pos, frag_pos, E);
    v3_norm(E);

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

    v3_add(Iamb, Idiff, frag_colour);
    v3_add(frag_colour, Ispec, frag_colour);

    v3_clamp(frag_colour, 0.0f, 1.0f);
}

static void TBN_create(vec3 normal, vec3 raw[3], vec2 tex[3], mat3 ret_TBN)
{
    // the tangent could be precomputed for all vertices to save doing it every frame
    vec3 edge1, edge2;
    v3_sub(raw[1], raw[0], edge1);
    v3_sub(raw[2], raw[0], edge2);

    vec2 delta_uv1, delta_uv2;
    v2_sub(tex[1], tex[0], delta_uv1);
    v2_sub(tex[2], tex[0], delta_uv2);

    const float f = 1.0f / (delta_uv1[0] * delta_uv2[1] - delta_uv2[0] * delta_uv1[1]);

    vec3 tangent;
    tangent[0] = f * (delta_uv2[1] * edge1[0] - delta_uv1[1] * edge2[0]);
    tangent[1] = f * (delta_uv2[1] * edge1[1] - delta_uv1[1] * edge2[1]);
    tangent[2] = f * (delta_uv2[1] * edge1[2] - delta_uv1[1] * edge2[2]);

    vec3 N;
    m3_mul_v3(normal_mapping_data.nrm_matrix, normal, N); // Transfrom from tangent to normal space
    v3_norm(N);

    vec3 T;
    m3_mul_v3(normal_mapping_data.nrm_matrix, tangent, T);
    v3_norm(T);

    /* Gram-Schmidt orthogonalize : t = normalize(t - n * dot(n, t)) */
    vec3 rhs;
    v3_scale(N, v3_dot(T, N), rhs);
    v3_sub(T, rhs, T);
    v3_norm(T);

    vec3 B;
    v3_cross(N, T, B);
    // v3_norm(B);

    /* Handedness */
    vec3 tmp;
    v3_cross(N, T, tmp);
    if (v3_dot(tmp, B) < 0.0f)
    {
        T[0] = -T[0];
        T[1] = -T[1];
        T[2] = -T[2];
    }

    //  Now our matrix is ready to take world coordinates, and put them into tangent space
    //  Transpose to make a TBN that can convert vertices into Tangent space
    //  Dont Transpose, and this lest us take vertices from Tangent space into world space

#if 1
    ret_TBN[0][0] = T[0];
    ret_TBN[1][0] = T[1];
    ret_TBN[2][0] = T[2];

    ret_TBN[0][1] = B[0];
    ret_TBN[1][1] = B[1];
    ret_TBN[2][1] = B[2];

    ret_TBN[0][2] = N[0];
    ret_TBN[1][2] = N[1];
    ret_TBN[2][2] = N[2];
#else
    ret_TBN[0][0] = T[0];
    ret_TBN[0][1] = T[1];
    ret_TBN[0][2] = T[2];

    ret_TBN[1][0] = B[0];
    ret_TBN[1][1] = B[1];
    ret_TBN[1][2] = B[2];

    ret_TBN[2][0] = N[0];
    ret_TBN[2][1] = N[1];
    ret_TBN[2][2] = N[2];
#endif
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

    // get world space vertex positions
    vec4 ws[3];
    m4_mul_v4(state.model, (vec4){raw[0][0], raw[0][1], raw[0][2], 1.0f}, ws[0]);
    m4_mul_v4(state.model, (vec4){raw[1][0], raw[1][1], raw[1][2], 1.0f}, ws[1]);
    m4_mul_v4(state.model, (vec4){raw[2][0], raw[2][1], raw[2][2], 1.0f}, ws[2]);

    // calculate bounding rectangle
    int AABB[4];
    AABB_make(proj, AABB);

    // normal mapping setup
    vec3 tangent_verts[3], tangent_light_pos[3], tangent_cam_pos[3];
    for (size_t i = 0; i < 3; i++) // Do we really do this for each vert?
    {
        mat3 TBN = {0};
        TBN_create(nrm[i], raw, texcoord, TBN);

        m3_mul_v3(TBN, ws[i], tangent_verts[i]);

        m3_mul_v3(TBN, LIGHT_POSITION, tangent_light_pos[i]);
        m3_mul_v3(TBN, (vec3){0.0f, 0.0f, 2.0f}, tangent_cam_pos[i]);
    }

    // Evaluaate edge equation at first tile origin
    const float edgeFunc0 = (E0[0] * (float)AABB[0]) + (E0[1] * (float)AABB[1]) + E0[2];
    const float edgeFunc1 = (E1[0] * (float)AABB[0]) + (E1[1] * (float)AABB[1]) + E1[2];
    const float edgeFunc2 = (E2[0] * (float)AABB[0]) + (E2[1] * (float)AABB[1]) + E2[2];

    /* Pre compute some tie breaker test results */
    const bool pre_comp_tie_E0 = Tie_Breaker_AB_Test(E0);
    const bool pre_comp_tie_E1 = Tie_Breaker_AB_Test(E1);
    const bool pre_comp_tie_E2 = Tie_Breaker_AB_Test(E2);

    // pre fetch tex data
    const int      texw    = state.tex.w;
    const int      texh    = state.tex.h;
    const int      texbpp  = state.tex.bpp;
    unsigned char *texdata = state.tex.data;

    const int      nrmw    = normal_mapping_data.nrm_tex.w;
    const int      nrmbpp  = normal_mapping_data.nrm_tex.bpp;
    unsigned char *nrmdata = normal_mapping_data.nrm_tex.data;

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

            const float depth = Interpolate_Values(littlef_values, (vec3){trans[0][2], trans[1][2], trans[2][2]});

            float *oldZ = pDepthBuffer + index;
            // const float invZ = 1.0f / depth;
            const float invZ = depth;

            // Perform depth test
            if (invZ > *oldZ)
                continue;

            // Depth test passed, update depth buffer value
            *oldZ = invZ;

            // Interpolate texture coordinates
            float u = Interpolate_Values(littlef_values, (vec3){texcoord[0][0], texcoord[1][0], texcoord[2][0]});
            float v = Interpolate_Values(littlef_values, (vec3){texcoord[0][1], texcoord[1][1], texcoord[2][1]});

            u *= (float)texw - 1;
            v *= (float)texh - 1;

            unsigned char *nrm_colour = nrmdata + (((int)u + nrmw * (int)v) * nrmbpp);
            vec3           frag_nrm;
            frag_nrm[0] = ((nrm_colour[0] / 255.0f) * 2.0f) - 1.0f; // x
            frag_nrm[1] = ((nrm_colour[1] / 255.0f) * 2.0f) - 1.0f; // y
            frag_nrm[2] = ((nrm_colour[2] / 255.0f) * 2.0f) - 1.0f; // z

            vec3 frag_pos, light_pos, cam_pos;
            for (size_t i = 0; i < 3; i++)
            {
                frag_pos[i]  = Interpolate_Values(littlef_values, (vec3){tangent_verts[0][i], tangent_verts[1][i], tangent_verts[2][i]});
                light_pos[i] = Interpolate_Values(littlef_values, (vec3){tangent_light_pos[0][i], tangent_light_pos[1][i], tangent_light_pos[2][i]});
                cam_pos[i]   = Interpolate_Values(littlef_values, (vec3){tangent_cam_pos[0][i], tangent_cam_pos[1][i], tangent_cam_pos[2][i]});
            }

            unsigned char *texcolour = texdata + (((int)u + texw * (int)v) * texbpp);

            vec3 diffuse_colour;
            diffuse_colour[0] = texcolour[0] / 255.0f;
            diffuse_colour[1] = texcolour[1] / 255.0f;
            diffuse_colour[2] = texcolour[2] / 255.0f;

            vec3 frag_colour;
            phong_lighting(diffuse_colour, frag_pos, cam_pos, light_pos, frag_nrm, frag_colour);

            unsigned char red = (unsigned char)(frag_colour[0] * 255.0f);
            unsigned char gre = (unsigned char)(frag_colour[1] * 255.0f);
            unsigned char blu = (unsigned char)(frag_colour[2] * 255.0f);

            uint32_t pixelcolour = ((uint32_t)0xFFU << 24) | ((uint32_t)blu << 16) | ((uint32_t)gre << 8) | (uint32_t)red;
            grafika_setpixel(x, y, pixelcolour);
        }
    }
}

static void draw_onexit(void)
{
    SAFE_FREE(transformed_vertices);
    tex_destroy(&normal_mapping_data.nrm_tex);
}

static void draw_object(void)
{
    mat4 cum_matrix;
    m4_mul_m4(VP_MATRIX, state.MVP, cum_matrix);

    mat4 tmp, nrm_mat;
    m4_inv(state.model, tmp);
    m4_transpose(tmp, nrm_mat);
    m3_from_m4(nrm_mat, normal_mapping_data.nrm_matrix);

    const obj_t object = state.obj;

    if (!transformed_vertices)
    {
        transformed_vertices        = malloc(sizeof(vec4) * object.num_pos);
        normal_mapping_data.nrm_tex = tex_load(object.mats[0].map_bump, true);
    }

    float *pos = object.pos;
#pragma omp parallel for
    // transform all verts to screen space
    for (size_t i = 0; i < object.num_pos; i++)
    {
        vec4 vertex = {pos[i * 3 + 0],
                       pos[i * 3 + 1],
                       pos[i * 3 + 2],
                       1.0f};
        m4_mul_v4(cum_matrix, vertex, &transformed_vertices[i * 4]);
    }

#pragma omp parallel
    {
        float         *pPos     = transformed_vertices;
        float         *pRawpos  = object.pos;
        float         *pNrm     = object.norms;
        float         *pTex     = object.texs;
        vertindices_t *pIndices = object.indices;

#pragma omp for
        for (size_t i = 0; i < object.num_f_rows; ++i)
        {
            vec4 trans[3];
            vec3 raw[3];
            vec3 nrm[3];
            vec2 texcoords[3];

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
#endif // __NORMAL_MAP_H__