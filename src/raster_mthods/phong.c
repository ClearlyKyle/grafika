#include "common.h"

#define PHONG_AMBI_AMOUNT (0.2f)   // 0.0f -> 1.0f
#define PHONG_SPEC_AMOUNT (0.4f)   // 0.0f -> 1.0f
#define PHONG_SHININESS   (128.0f) // pow of 2 pls
#define LIGHT_POSITION    ((vec3){0.5f, 0.5f, 5.0f})

struct phong
{
    mat3 nrm_matrix;
};

static struct phong phong_data           = {0};
struct rasterstate  raster_state         = {0};
static float       *transformed_vertices = NULL;

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

static void phong_lighting(vec3 diff_colour, vec3 frag_pos, vec3 N, vec3 frag_colour)
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

    v3_add(Iamb, Idiff, frag_colour);
    v3_add(frag_colour, Ispec, frag_colour);

    v3_clamp(frag_colour, 0.0f, 1.0f);
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

    // transform our normals using the "normals matrix" this is needed due to
    // the model matrix having the ability to scale and scew the model, and our
    // normal do not like this
    vec3 new_nrm[3] = {0};
    m3_mul_v3(phong_data.nrm_matrix, nrm[0], new_nrm[0]);
    m3_mul_v3(phong_data.nrm_matrix, nrm[1], new_nrm[1]);
    m3_mul_v3(phong_data.nrm_matrix, nrm[2], new_nrm[2]);

    vec3 nrm_x = {new_nrm[0][0], new_nrm[1][0], new_nrm[2][0]};
    vec3 nrm_y = {new_nrm[0][1], new_nrm[1][1], new_nrm[2][1]};
    vec3 nrm_z = {new_nrm[0][2], new_nrm[1][2], new_nrm[2][2]};

    // get world space vertex positions
    vec4 ws[3] = {0};
    m4_mul_v4(raster_state.model, (vec4){raw[0][0], raw[0][1], raw[0][2], 1.0f}, ws[0]);
    m4_mul_v4(raster_state.model, (vec4){raw[1][0], raw[1][1], raw[1][2], 1.0f}, ws[1]);
    m4_mul_v4(raster_state.model, (vec4){raw[2][0], raw[2][1], raw[2][2], 1.0f}, ws[2]);

    vec3 ws_x = {ws[0][0], ws[1][0], ws[2][0]};
    vec3 ws_y = {ws[0][1], ws[1][1], ws[2][1]};
    vec3 ws_z = {ws[0][2], ws[1][2], ws[2][2]};

    // calculate bounding rectangle
    int AABB[4] = {0};
    AABB_make(proj, AABB);

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

    vec3 u_values = {texcoord[0][0], texcoord[1][0], texcoord[2][0]};
    vec3 v_values = {texcoord[0][1], texcoord[1][1], texcoord[2][1]};

    vec3   z_values     = {trans[0][2], trans[1][2], trans[2][2]};
    float *depth_buffer = rend.depth_buffer;

    // Start rasterizing by looping over pixels to output a per-pixel color
    for (int y = AABB[1], step_y = 0; y <= AABB[3]; y++, step_y++)
    {
        for (int x = AABB[0], step_x = 0; x <= AABB[2]; x++, step_x++)
        {
            // step from edge function by multiples of the edge values
            // edge_value * step_multiple

            // evaluate edge functions at current fragment
            // E(x + s, y + t) = E(x, y) + sa + tb,
            const float eval_edge_func_0 = edge_func_0 + ((E0[0] * (float)step_x) + (E0[1] * (float)step_y));
            const float eval_edge_func_1 = edge_func_1 + ((E1[0] * (float)step_x) + (E1[1] * (float)step_y));
            const float eval_edge_func_2 = edge_func_2 + ((E2[0] * (float)step_x) + (E2[1] * (float)step_y));

            // check if the current point is inside a traingle using Tie breaker rules
            if (edge_tie_breaker(eval_edge_func_0, pre_comp_tie_E0)) continue;
            if (edge_tie_breaker(eval_edge_func_1, pre_comp_tie_E1)) continue;
            if (edge_tie_breaker(eval_edge_func_2, pre_comp_tie_E2)) continue;

            // calculate interpolation coefficients
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

            vec3 frag_pos;
            frag_pos[0] = interpolate_values(littlef_values, ws_x); // x
            frag_pos[1] = interpolate_values(littlef_values, ws_y); // y
            frag_pos[2] = interpolate_values(littlef_values, ws_z); // z

            vec3 frag_nrm;
            frag_nrm[0] = interpolate_values(littlef_values, nrm_x); // x
            frag_nrm[1] = interpolate_values(littlef_values, nrm_y); // y
            frag_nrm[2] = interpolate_values(littlef_values, nrm_z); // z

            float u = interpolate_values(littlef_values, u_values);
            float v = interpolate_values(littlef_values, v_values);

            u = SDL_clamp(u, 0.0f, 1.0f);
            v = SDL_clamp(v, 0.0f, 1.0f);

            u *= (float)tex_w - 1;
            v *= (float)tex_h - 1;

            unsigned char *tex_colour = tex_data + (((int)u + tex_w * (int)v) * tex_bpp);

            vec3        diffuse_colour;
            const float inv_255 = 1.0f / 255.0f;
            diffuse_colour[0]   = tex_colour[0] * inv_255;
            diffuse_colour[1]   = tex_colour[1] * inv_255;
            diffuse_colour[2]   = tex_colour[2] * inv_255;

            vec3 frag_colour = {0};
            phong_lighting(diffuse_colour, frag_pos, frag_nrm, frag_colour);

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

    mat3 tmp = {0};
    m3_from_m4(raster_state.model, tmp);
    m3_inv(tmp, tmp);
    m3_transpose(tmp, phong_data.nrm_matrix);

    const struct obj obj = raster_state.obj;

    if (!transformed_vertices)
        transformed_vertices = arena_alloc_aligned(arena, sizeof(vec4) * obj.num_pos, 16);

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