#define GRAFIKA_TITLE ("grafika - matrix")
#include "common.h"

struct triangle
{
    vec3 pos[3];
    vec2 tex[3];
};

struct matrix
{
    mat4       model;
    mat4       MVP;
    vec3       cam_pos;
    struct tex tex;
    struct obj obj;
};

struct matrix raster_state = {0};

// use this to precompute the tie braker edge conditions, use the E value before looping through the pixels
static inline bool tie_breaker_ab_test(const vec3 E)
{
    return (E[0] < 0.0f) || (E[0] == 0.0f && E[1] <= 0.0f);
}

// apply tie-breaking rules on shared vertices in order to avoid double-shading fragments
static inline bool edge_tie_breaker(const float edge_value, const bool ab_test_result)
{
    return (edge_value < 0.0f) || (edge_value == 0.0f && ab_test_result);
}

static inline float interpolate_values(vec3 little_f_values, vec3 attribute)
{
    return v3_dot(little_f_values, attribute);
}

static void draw_AABB_outline(int AABB[4], uint32_t colour)
{
    for (int x = AABB[0]; x <= AABB[2]; ++x)
    {
        grafika_setpixel((uint32_t)x, (uint32_t)AABB[1], colour); // top
        grafika_setpixel((uint32_t)x, (uint32_t)AABB[3], colour); // bottom
    }

    for (int y = AABB[1]; y <= AABB[3]; ++y)
    {
        grafika_setpixel((uint32_t)AABB[0], (uint32_t)y, colour); // left
        grafika_setpixel((uint32_t)AABB[2], (uint32_t)y, colour); // right
    }
}

static void draw_triangle(const struct triangle t)
{
    // > Local space (or object space, loaded from OBJ file)
    //      + Model Matrix
    // > World space
    //      + View Matrix
    // > View space
    //      + Projection Matrix
    // > Clip space

    // we can first multiply together all the matricies and jump straight from "local space"
    // directly to "clip space"

    // "local space" to "clip space"
    vec4 clip_space[3] = {0};
    for (int i = 0; i < 3; ++i)
    {
        vec4 pos = {t.pos[i][0], t.pos[i][1], t.pos[i][2], 1.0f};
        m4_mul_v4(raster_state.MVP, pos, clip_space[i]);

        // clipping
        // -w ≤ x ≤ w
        // -w ≤ y ≤ w
        // -w ≤ z ≤ w
        // if (clip_space[i][0] < -clip_space[i][3] || clip_space[i][0] > clip_space[i][3] ||
        //    clip_space[i][1] < -clip_space[i][3] || clip_space[i][1] > clip_space[i][3] ||
        //    clip_space[i][2] < -clip_space[i][3] || clip_space[i][2] > clip_space[i][3])
        //{
        //    return; // triangle is outside
        //}
    }

    // after clipping has taken place, our x, y, z values in the range of [-w, w]
    // NOTE : normally here, we would perform "perspective-divison" (x/w, y/w, z/w) to get
    //  them into "ndc space" which is the range [-1, 1]

    // instead, in this "matrix method", we avoid the perspective division

    // so what is going on here?
    // normally, after perspective division, we convert from ndc to screen space.
    // that is : ndc [-1, 1] -> screen space [0, width|height]
    // we do this by doing the following
    //     screen_space.x = (ndc.x + 1.0f) * 0.5f * GRAFIKA_SCREEN_WIDTH;
    //     screen_space.y = (ndc.y + 1.0f) * 0.5f * GRAFIKA_SCREEN_HEIGHT;
    // instead we are going directly from clip space into screen space

    vec4 screen_space[3];
    for (int i = 0; i < 3; ++i)
    {
        screen_space[i][0] = GRAFIKA_SCREEN_WIDTH * (clip_space[i][0] + clip_space[i][3]) * 0.5f;
        screen_space[i][1] = GRAFIKA_SCREEN_HEIGHT * (clip_space[i][3] - clip_space[i][1]) * 0.5f;
        screen_space[i][2] = clip_space[i][2];
        screen_space[i][3] = clip_space[i][3];
    }

    // from the paper, calculate out "A" matrix
    const float a0 = (screen_space[1][1] * screen_space[2][3]) - (screen_space[2][1] * screen_space[1][3]);
    const float a1 = (screen_space[2][1] * screen_space[0][3]) - (screen_space[0][1] * screen_space[2][3]);
    const float a2 = (screen_space[0][1] * screen_space[1][3]) - (screen_space[1][1] * screen_space[0][3]);

    const float b0 = (screen_space[2][0] * screen_space[1][3]) - (screen_space[1][0] * screen_space[2][3]);
    const float b1 = (screen_space[0][0] * screen_space[2][3]) - (screen_space[2][0] * screen_space[0][3]);
    const float b2 = (screen_space[1][0] * screen_space[0][3]) - (screen_space[0][0] * screen_space[1][3]);

    const float c0 = (screen_space[1][0] * screen_space[2][1]) - (screen_space[2][0] * screen_space[1][1]);
    const float c1 = (screen_space[2][0] * screen_space[0][1]) - (screen_space[0][0] * screen_space[2][1]);
    const float c2 = (screen_space[0][0] * screen_space[1][1]) - (screen_space[1][0] * screen_space[0][1]);

    const float detM = (c0 * screen_space[0][3]) + (c1 * screen_space[1][3]) + (c2 * screen_space[2][3]);

    // The sign of the determinant gives
    // the orientation of the triangle: a CCW order gives a positive determinant
    // if (detM >= 0.0f) return; // reject CCW triangles
    if (detM < 0.0f) return; // reject CW triangles

    // set up edge functions (this is A = adj(M))
    vec3 E0 = {a0, b0, c0};
    vec3 E1 = {a1, b1, c1};
    vec3 E2 = {a2, b2, c2};

    vec3 bounds[3];
    for (int i = 0; i < 3; ++i)
    {
        bounds[i][0] = screen_space[i][0] / screen_space[i][3];
        bounds[i][1] = screen_space[i][1] / screen_space[i][3];
    }

    // calculate bounding rectangle
    int AABB[4] = {0};
    AABB_make(bounds, AABB);

    // draw_AABB_outline(AABB, 0x00FF00FF);
    // LOG("AABB : (%d, %d) to (%d, %d)\n", AABB[0], AABB[1], AABB[2], AABB[3]);

    // evaluaate edge equation at first tile origin
    const float edge_func_origin_0 = (E0[0] * (float)AABB[0]) + (E0[1] * (float)AABB[1]) + E0[2];
    const float edge_func_origin_1 = (E1[0] * (float)AABB[0]) + (E1[1] * (float)AABB[1]) + E1[2];
    const float edge_func_origin_2 = (E2[0] * (float)AABB[0]) + (E2[1] * (float)AABB[1]) + E2[2];

    // pre compute some tie breaker test results
    const bool pre_comp_tie_E0 = tie_breaker_ab_test(E0);
    const bool pre_comp_tie_E1 = tie_breaker_ab_test(E1);
    const bool pre_comp_tie_E2 = tie_breaker_ab_test(E2);

    // pre fetch tex data
    const int      tex_w    = raster_state.tex.w;
    const int      tex_h    = raster_state.tex.h;
    const int      tex_bpp  = raster_state.tex.bpp;
    unsigned char *tex_data = raster_state.tex.data;

    vec3 u_values = {t.tex[0][0], t.tex[1][0], t.tex[2][0]};
    vec3 v_values = {t.tex[0][1], t.tex[1][1], t.tex[2][1]};

    // NOTE : do these Z values need to be divided by the W?
    vec3   z_values     = {screen_space[0][2], screen_space[1][2], screen_space[2][2]};
    float *depth_buffer = rend.depth_buffer;

    for (int y = AABB[1], step_y = 0; y <= AABB[3]; y++, step_y++)
    {
        for (int x = AABB[0], step_x = 0; x <= AABB[2]; x++, step_x++)
        {
            // step from edge function by multiples of the edge values
            // edge_value * step_multiple

            // evaluate edge functions at current fragment
            // E(x + s, y + t) = E(x, y) + sa + tb,
            const float eval_edge_func_0 = edge_func_origin_0 + ((E0[0] * (float)(step_x)) + (E0[1] * (float)(step_y)));
            const float eval_edge_func_1 = edge_func_origin_1 + ((E1[0] * (float)(step_x)) + (E1[1] * (float)(step_y)));
            const float eval_edge_func_2 = edge_func_origin_2 + ((E2[0] * (float)(step_x)) + (E2[1] * (float)(step_y)));

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
            const int index = (x * GRAFIKA_SCREEN_WIDTH) + y; // NOTE : is this wrong?

            const float depth = interpolate_values(littlef_values, z_values);

            float      *oldZ = depth_buffer + index;
            const float invZ = depth;

            // depth test
            if (invZ > *oldZ) continue;
            *oldZ = invZ;

            // interpolate texture coordinates
            float u = interpolate_values(littlef_values, u_values);
            float v = interpolate_values(littlef_values, v_values);

            u = SDL_clamp(u, 0.0f, 1.0f);
            v = SDL_clamp(v, 0.0f, 1.0f);

            u *= (float)tex_w - 1;
            v *= (float)tex_h - 1;

            unsigned char *tex_colour = tex_data + (((int)u + tex_w * (int)v) * tex_bpp);

            uint32_t pixel_colour = ((uint32_t)tex_colour[0] << 16) |
                                    ((uint32_t)tex_colour[1] << 8) |
                                    (uint32_t)tex_colour[2];

            grafika_setpixel((uint32_t)x, (uint32_t)y, pixel_colour);
        }
    }
}

void draw_onstart(struct arena *arena)
{
    UNUSED(arena);

    struct obj obj = raster_state.obj;

    ASSERT(obj.mats, "Object must have atleast a diffuse texture\n");
    raster_state.tex = tex_load(obj.mats[0].map_Kd);
}

void draw_object(struct arena *arena)
{
    UNUSED(arena);

#pragma omp parallel
    {
        struct obj          obj         = raster_state.obj;
        float              *obj_pos     = obj.pos;
        float              *obj_tex     = obj.texs;
        struct vertindices *obj_indices = obj.indices;

#pragma omp for
        for (size_t i = 0; i < obj.num_f_rows; ++i)
        {
            struct triangle t = {0};

            for (size_t j = 0; j < 3; ++j)
            {
                struct vertindices indices = obj_indices[i * 3 + j];

                const int pos_index = 3 * indices.v_idx;
                const int tex_index = 2 * indices.vt_idx;

#if 1
                const size_t store_index = j; // CCW triangles
#else
                const size_t store_index = 2 - j; // CW triangles (which blender uses)
#endif
                t.pos[store_index][0] = obj_pos[pos_index + 0];
                t.pos[store_index][1] = obj_pos[pos_index + 1];
                t.pos[store_index][2] = obj_pos[pos_index + 2];

                t.tex[store_index][0] = obj_tex[tex_index + 0];
                t.tex[store_index][1] = obj_tex[tex_index + 1];
            }

            draw_triangle(t);
        }
    }
}

void draw_onexit(void)
{
    tex_destroy(&raster_state.tex);
    // raster_state.obj       // cleaned up by arena
}