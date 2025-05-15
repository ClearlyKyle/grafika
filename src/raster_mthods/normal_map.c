#include "common.h"

#define PHONG_AMBI_AMOUNT (0.2f)  // 0.0f -> 1.0f
#define PHONG_SPEC_AMOUNT (0.4f)  // 0.0f -> 1.0f
#define PHONG_SHININESS   (64.0f) // pow of 2 pls
#define LIGHT_POSITION    ((vec3){0.5f, 0.5f, 5.0f})

struct nrm_map
{
    mat3  nrm_matrix;
    mat3  TBN;
    tex_t nrm_tex;
};

static struct nrm_map normal_mapping_data  = {0};
struct rasterstate    raster_state         = {0};
static float         *transformed_vertices = NULL;

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

static inline bool tie_breaker_ab_test(const vec3 E)
{
    return (E[0] > 0.0f) || (E[0] == 0.0f && E[1] >= 0.0f);
}

static inline bool edge_tie_breaker(const float edge_value, const bool ab_test_result)
{
    return (edge_value > 0.0f) || (edge_value == 0.0f && ab_test_result);
}

static inline float interpolate_values(vec3 little_f_values, vec3 attribute)
{
    return v3_dot(little_f_values, attribute);
}

static float calculate_specular_amount(vec3 L, vec3 E, vec3 N, const float shininess)
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
        // PHONG
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

static void phong_lighting(vec3 diff_colour, vec3 frag_pos, vec3 cam_pos, vec3 light_pos, vec3 N, vec3 frag_colour)
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
    m3_mul_v3(normal_mapping_data.nrm_matrix, normal, N); // transfrom from tangent to normal space
    v3_norm(N);

    vec3 T;
    m3_mul_v3(normal_mapping_data.nrm_matrix, tangent, T);
    v3_norm(T);

    // Gram-Schmidt orthogonalize : t = normalize(t - n * dot(n, t))
    vec3 rhs;
    v3_scale(N, v3_dot(T, N), rhs);
    v3_sub(T, rhs, T);
    v3_norm(T);

    vec3 B;
    v3_cross(N, T, B);
    // v3_norm(B);

    // handedness
    vec3 tmp;
    v3_cross(N, T, tmp);
    if (v3_dot(tmp, B) < 0.0f)
    {
        T[0] = -T[0];
        T[1] = -T[1];
        T[2] = -T[2];
    }

    //  now our matrix is ready to take world coordinates, and put them into tangent space
    //  transpose to make a TBN that can convert vertices into Tangent space
    //  dont Transpose, and this lest us take vertices from Tangent space into world space

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
        proj[i][0] = trans[i][0] / trans[i][3];
        proj[i][1] = trans[i][1] / trans[i][3];
        proj[i][2] = trans[i][2] / trans[i][3];
    }

    // get world space vertex positions
    vec4 ws[3] = {0};
    m4_mul_v4(raster_state.model, (vec4){raw[0][0], raw[0][1], raw[0][2], 1.0f}, ws[0]);
    m4_mul_v4(raster_state.model, (vec4){raw[1][0], raw[1][1], raw[1][2], 1.0f}, ws[1]);
    m4_mul_v4(raster_state.model, (vec4){raw[2][0], raw[2][1], raw[2][2], 1.0f}, ws[2]);

    // calculate bounding rectangle
    int AABB[4] = {0};
    AABB_make(proj, AABB);

    // normal mapping setup
    vec3 tangent_verts[3], tangent_light_pos[3], tangent_cam_pos[3];
    for (size_t i = 0; i < 3; i++) // NOTE : do we really do this for each vert?
    {
        mat3 TBN = {0};
        TBN_create(nrm[i], raw, texcoord, TBN);

        m3_mul_v3(TBN, ws[i], tangent_verts[i]);

        m3_mul_v3(TBN, LIGHT_POSITION, tangent_light_pos[i]);
        m3_mul_v3(TBN, raster_state.cam_pos, tangent_cam_pos[i]);
    }

    // evaluaate edge equation at first tile origin
    const float edge_func_0 = (E0[0] * (float)AABB[0]) + (E0[1] * (float)AABB[1]) + E0[2];
    const float edge_func_1 = (E1[0] * (float)AABB[0]) + (E1[1] * (float)AABB[1]) + E1[2];
    const float edge_func_2 = (E2[0] * (float)AABB[0]) + (E2[1] * (float)AABB[1]) + E2[2];

    // pre compute some tie breaker test results
    const bool pre_comp_tie_E0 = tie_breaker_ab_test(E0);
    const bool pre_comp_tie_E1 = tie_breaker_ab_test(E1);
    const bool pre_comp_tie_E2 = tie_breaker_ab_test(E2);

    // pre fetch tex data
    const int      tex_w    = raster_state.tex.w;
    const int      tex_h    = raster_state.tex.h;
    const int      tex_bpp  = raster_state.tex.bpp;
    unsigned char *tex_data = raster_state.tex.data;

    const int      nrm_w    = normal_mapping_data.nrm_tex.w;
    const int      nrm_bpp  = normal_mapping_data.nrm_tex.bpp;
    unsigned char *nrm_data = normal_mapping_data.nrm_tex.data;

    vec3 u_values = {texcoord[0][0], texcoord[1][0], texcoord[2][0]};
    vec3 v_values = {texcoord[0][1], texcoord[1][1], texcoord[2][1]};

    vec3   z_values     = {trans[0][2], trans[1][2], trans[2][2]};
    float *depth_buffer = rend.depth_buffer;

    for (int y = AABB[1], step_y = 0; y <= AABB[3]; y++, step_y++)
    {
        for (int x = AABB[0], step_x = 0; x <= AABB[2]; x++, step_x++)
        {
            // E(x + s, y + t) = E(x, y) + sa + tb,
            const float eval_edge_func_0 = edge_func_0 + ((E0[0] * (float)step_x) + (E0[1] * (float)step_y));
            const float eval_edge_func_1 = edge_func_1 + ((E1[0] * (float)step_x) + (E1[1] * (float)step_y));
            const float eval_edge_func_2 = edge_func_2 + ((E2[0] * (float)step_x) + (E2[1] * (float)step_y));

            if (edge_tie_breaker(eval_edge_func_0, pre_comp_tie_E0)) continue;
            if (edge_tie_breaker(eval_edge_func_1, pre_comp_tie_E1)) continue;
            if (edge_tie_breaker(eval_edge_func_2, pre_comp_tie_E2)) continue;

            const float F0 = eval_edge_func_0;
            const float F1 = eval_edge_func_1;
            const float F2 = eval_edge_func_2;

            const float r = 1.0f / (F0 + F1 + F2);

            vec3 littlef_values = {eval_edge_func_0 * r, eval_edge_func_1 * r, eval_edge_func_2 * r};

            // interpolate depth value
            const int index = (x * GRAFIKA_SCREEN_WIDTH) + y;

            const float depth = interpolate_values(littlef_values, z_values);

            float      *oldZ = depth_buffer + index;
            const float invZ = depth;

            // depth test
            if (invZ > *oldZ) continue;
            *oldZ = invZ;

            float u = interpolate_values(littlef_values, u_values);
            float v = interpolate_values(littlef_values, v_values);

            u = SDL_clamp(u, 0.0f, 1.0f);
            v = SDL_clamp(v, 0.0f, 1.0f);

            u *= (float)tex_w - 1;
            v *= (float)tex_h - 1;

            unsigned char *nrm_colour = nrm_data + (((int)u + nrm_w * (int)v) * nrm_bpp);
            vec3           frag_nrm;
            const float    inv_255 = 1.0f / 255.0f;
            frag_nrm[0]            = ((nrm_colour[0] * inv_255) * 2.0f) - 1.0f; // x
            frag_nrm[1]            = ((nrm_colour[1] * inv_255) * 2.0f) - 1.0f; // y
            frag_nrm[2]            = ((nrm_colour[2] * inv_255) * 2.0f) - 1.0f; // z

            vec3 frag_pos, light_pos, cam_pos;
            for (size_t i = 0; i < 3; i++)
            {
                frag_pos[i]  = interpolate_values(littlef_values, (vec3){tangent_verts[0][i], tangent_verts[1][i], tangent_verts[2][i]});
                light_pos[i] = interpolate_values(littlef_values, (vec3){tangent_light_pos[0][i], tangent_light_pos[1][i], tangent_light_pos[2][i]});
                cam_pos[i]   = interpolate_values(littlef_values, (vec3){tangent_cam_pos[0][i], tangent_cam_pos[1][i], tangent_cam_pos[2][i]});
            }

            unsigned char *tex_colour = tex_data + (((int)u + tex_w * (int)v) * tex_bpp);

            vec3 diffuse_colour;
            diffuse_colour[0] = tex_colour[0] * inv_255;
            diffuse_colour[1] = tex_colour[1] * inv_255;
            diffuse_colour[2] = tex_colour[2] * inv_255;

            vec3 frag_colour = {0};
            phong_lighting(diffuse_colour, frag_pos, cam_pos, light_pos, frag_nrm, frag_colour);

            unsigned char red = (unsigned char)(frag_colour[0] * 255.0f);
            unsigned char gre = (unsigned char)(frag_colour[1] * 255.0f);
            unsigned char blu = (unsigned char)(frag_colour[2] * 255.0f);

            uint32_t pixel_colour = ((uint32_t)red << 16) | ((uint32_t)gre << 8) | (uint32_t)blu;

            grafika_setpixel(x, y, pixel_colour);
        }
    }
}

void draw_object(struct arena *arena)
{
    mat4 cum_matrix = {0};
    m4_mul_m4(VP_MATRIX, raster_state.MVP, cum_matrix);

    mat4 tmp = {0}, nrm_mat = {0};
    m4_inv(raster_state.model, tmp);
    m4_transpose(tmp, nrm_mat);
    m3_from_m4(nrm_mat, normal_mapping_data.nrm_matrix);

    const struct obj obj = raster_state.obj;

    if (!transformed_vertices)
    {
        transformed_vertices = arena_alloc_aligned(arena, sizeof(vec4) * obj.num_pos, 16);

        ASSERT(obj.mats[0].map_bump, "Object must have a normal map!\n");
        normal_mapping_data.nrm_tex = tex_load(obj.mats[0].map_bump);
    }

#pragma omp parallel
    {
        float         *tran_pos    = transformed_vertices;
        float         *obj_pos     = obj.pos;
        float         *obj_nrm     = obj.norms;
        float         *obj_tex     = obj.texs;
        struct vertindices *obj_indices = obj.indices;

#pragma omp for
        // transform all verts to screen space
        for (size_t i = 0; i < obj.num_pos; i++)
        {
            vec4 vertex = {obj_pos[i * 3 + 0],
                           obj_pos[i * 3 + 1],
                           obj_pos[i * 3 + 2],
                           1.0f};
            m4_mul_v4(cum_matrix, vertex, &tran_pos[i * 4]);
        }

#pragma omp for
        for (size_t i = 0; i < obj.num_f_rows; ++i)
        {
            vec4 trans[3]     = {0};
            vec3 raw[3]       = {0};
            vec3 nrm[3]       = {0};
            vec2 texcoords[3] = {0};

            for (size_t j = 0; j < 3; ++j)
            {
                struct vertindices indices = obj_indices[i * 3 + j];

                const int pos_index = indices.v_idx;
                const int nrm_index = indices.vn_idx * 3;
                const int tex_index = indices.vt_idx * 2;
#if 0
                const size_t store_index = j; // CCW triangles
#else
                const size_t store_index = 2 - j; // CW triangles (which blender uses)
#endif
                float *p              = tran_pos + 4 * pos_index;
                float *raw_p          = obj_pos + 3 * pos_index;
                trans[store_index][0] = p[0];
                trans[store_index][1] = p[1];
                trans[store_index][2] = p[2];
                trans[store_index][3] = p[3];

                raw[store_index][0] = raw_p[0];
                raw[store_index][1] = raw_p[1];
                raw[store_index][2] = raw_p[2];

                nrm[store_index][0] = obj_nrm[nrm_index + 0];
                nrm[store_index][1] = obj_nrm[nrm_index + 1];
                nrm[store_index][2] = obj_nrm[nrm_index + 2];

                texcoords[store_index][0] = obj_tex[tex_index + 0];
                texcoords[store_index][1] = obj_tex[tex_index + 1];
            }
            draw_triangle(trans, raw, nrm, texcoords);
        }
    }
}

void draw_onexit(void)
{
    tex_destroy(&normal_mapping_data.nrm_tex);
}