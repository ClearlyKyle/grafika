#define GRAFIKA_TITLE ("grafika - binning_simd_jobs")
#include "common.h"

#define RABOTY_IMPLEMENTATION
#include "../raboty.h"

#define MAX_TRIANGLES_PER_TILE (2048)

#define TILE_SIZE_32x32   (32u)
#define TILE_SIZE_64x64   (64u)
#define TILE_SIZE_128x128 (128u)
#define TILE_SIZE_256x256 (256u)
#define TILE_SIZE_512x512 (512u)

#define RENDER_TILE_SIZE     TILE_SIZE_64x64
#define RENDER_TILE_SIZE_INV (1.0f / RENDER_TILE_SIZE) // RENDER_TILE_SIZE must be power of 2

#define RENDER_NUM_TILE_X (GRAFIKA_SCREEN_WIDTH / RENDER_TILE_SIZE)
#define RENDER_NUM_TILE_Y (GRAFIKA_SCREEN_HEIGHT / RENDER_TILE_SIZE)

#define RENDER_TILE_COUNT (RENDER_NUM_TILE_X * RENDER_NUM_TILE_Y)

#define GRAFIKA_SCREEN_HALF_WIDTH  (0.5f * GRAFIKA_SCREEN_WIDTH)
#define GRAFIKA_SCREEN_HALF_HEIGHT (0.5f * GRAFIKA_SCREEN_HEIGHT)

struct aabb
{
    float minx, miny, maxx, maxy;
};

enum raster_type
{
    RASTER_FULL_TILE,
    RASTER_PARTIAL_TILE,
};

struct triangle_ref
{
    enum raster_type type;
    uint32_t         triangle_index;
};

struct tile
{
    struct triangle_ref triangles[MAX_TRIANGLES_PER_TILE]; // TODO : should be dynamic?
    uint32_t            triangle_count;
    uint8_t             has_been_queued;
};

struct tile_coords
{
    uint8_t tile_x, tile_y;
};

struct tile_coords tile_coords[RENDER_TILE_COUNT]     = {0};
struct tile        tile_array[RENDER_TILE_COUNT]      = {0};
uint32_t           tiles_with_bins[RENDER_TILE_COUNT] = {0};
uint32_t           tiles_queue_index                  = 0;

struct binning_simd_jobs
{
    mat4       model;
    mat4       MVP;
    vec3       cam_pos;
    struct tex tex;
    struct obj obj;
};

struct binning_simd_jobs raster_state = {0};

static inline __m128 interpolate_values(const __m128 F[3], const __m128 V[3])
{
    return _mm_add_ps(
        _mm_mul_ps(V[2], F[2]),
        _mm_add_ps(
            _mm_mul_ps(V[1], F[1]),
            _mm_mul_ps(V[0], F[0])));
}

static inline void tile_add_data(uint32_t tile_index, struct triangle_ref ref)
{
    struct tile *tile = &tile_array[tile_index];

    uint32_t index = tile->triangle_count++;
    ASSERT(index < MAX_TRIANGLES_PER_TILE, "Too many triangles in this bin\n");

    tile->triangles[index] = ref;

    if (tile->has_been_queued) return;

    tile->has_been_queued                = 1;
    tiles_with_bins[tiles_queue_index++] = tile_index;
}

static bool _triangle_clipping(vec4 pos[3])
{
    enum
    {
        PLANE_LEFT   = 1 << 0,
        PLANE_RIGHT  = 1 << 1,
        PLANE_BOTTOM = 1 << 2,
        PLANE_TOP    = 1 << 3,
        PLANE_NEAR   = 1 << 4,
        PLANE_FAR    = 1 << 5
    };

    uint8_t mask0 = 0;
    if (pos[0][0] < -pos[0][3]) mask0 |= PLANE_LEFT;
    if (pos[0][0] > pos[0][3]) mask0 |= PLANE_RIGHT;
    if (pos[0][1] < -pos[0][3]) mask0 |= PLANE_BOTTOM;
    if (pos[0][1] > pos[0][3]) mask0 |= PLANE_TOP;
    if (pos[0][2] < 0.0f) mask0 |= PLANE_NEAR;
    if (pos[0][2] > pos[0][3]) mask0 |= PLANE_FAR;

    uint8_t mask1 = 0;
    if (pos[1][0] < -pos[1][3]) mask1 |= PLANE_LEFT;
    if (pos[1][0] > pos[1][3]) mask1 |= PLANE_RIGHT;
    if (pos[1][1] < -pos[1][3]) mask1 |= PLANE_BOTTOM;
    if (pos[1][1] > pos[1][3]) mask1 |= PLANE_TOP;
    if (pos[1][2] < 0.0f) mask1 |= PLANE_NEAR;
    if (pos[1][2] > pos[1][3]) mask1 |= PLANE_FAR;

    uint8_t mask2 = 0;
    if (pos[2][0] < -pos[2][3]) mask2 |= PLANE_LEFT;
    if (pos[2][0] > pos[2][3]) mask2 |= PLANE_RIGHT;
    if (pos[2][1] < -pos[2][3]) mask2 |= PLANE_BOTTOM;
    if (pos[2][1] > pos[2][3]) mask2 |= PLANE_TOP;
    if (pos[2][2] < 0.0f) mask2 |= PLANE_NEAR;
    if (pos[2][2] > pos[2][3]) mask2 |= PLANE_FAR;

    uint8_t and_mask = mask0 & mask1 & mask2;
    uint8_t or_mask  = mask0 | mask1 | mask2;

    if (and_mask != 0)
        return false; // CLIP_OUTSIDE : all vertices outside same plane
    else if (or_mask == 0)
        return true; // CLIP_INSIDE : all vertices inside all planes
    else
        return true; // CLIP_PARTIAL : partial overlap, needs clipping
}

#define GRAFIKA_SCREEN_HALF_WIDTH  (0.5f * GRAFIKA_SCREEN_WIDTH)
#define GRAFIKA_SCREEN_HALF_HEIGHT (0.5f * GRAFIKA_SCREEN_HEIGHT)

static bool _triangle_setup_and_cull(vec4 pos[3], vec3 E[3], vec3 Z, struct aabb *aabb)
{
    // first, transform clip-space (x, y, z, w) vertices to device-space 2D homogeneous coordinates (x, y, w)
    vec4 homo[3];
    for (size_t i = 0; i < 3; i++)
    {
        homo[i][0] = (pos[i][0] + pos[i][3]) * GRAFIKA_SCREEN_HALF_WIDTH;
        homo[i][1] = (pos[i][3] - pos[i][1]) * GRAFIKA_SCREEN_HALF_HEIGHT; // origin Top Left
        homo[i][2] = pos[i][2];
        homo[i][3] = pos[i][3];
    }

    // from the paper, calculate out "A" matrix
    const float a0 = (homo[1][1] * homo[2][3]) - (homo[2][1] * homo[1][3]);
    const float a1 = (homo[2][1] * homo[0][3]) - (homo[0][1] * homo[2][3]);
    const float a2 = (homo[0][1] * homo[1][3]) - (homo[1][1] * homo[0][3]);

    const float b0 = (homo[2][0] * homo[1][3]) - (homo[1][0] * homo[2][3]);
    const float b1 = (homo[0][0] * homo[2][3]) - (homo[2][0] * homo[0][3]);
    const float b2 = (homo[1][0] * homo[0][3]) - (homo[0][0] * homo[1][3]);

    const float c0 = (homo[1][0] * homo[2][1]) - (homo[2][0] * homo[1][1]);
    const float c1 = (homo[2][0] * homo[0][1]) - (homo[0][0] * homo[2][1]);
    const float c2 = (homo[0][0] * homo[1][1]) - (homo[1][0] * homo[0][1]);

    // det(M) <= 0  -> back-facing triangle
    float detM = (c0 * homo[0][3]) + (c1 * homo[1][3]) + (c2 * homo[2][3]);

    if (detM <= 0.0f) return false; // reject CW triangles

    vec2 bounds[3];
    for (int i = 0; i < 3; ++i)
    {
        bounds[i][0] = homo[i][0] / homo[i][3];
        bounds[i][1] = homo[i][1] / homo[i][3];
    }

    const float screen_width  = GRAFIKA_SCREEN_WIDTH - 1;
    const float screen_height = GRAFIKA_SCREEN_HEIGHT - 1;

    aabb->minx = max(0.0f, min(bounds[0][0], min(bounds[1][0], bounds[2][0])));
    aabb->maxx = min(screen_width, max(bounds[0][0], max(bounds[1][0], bounds[2][0])));

    aabb->miny = max(0.0f, min(bounds[0][1], min(bounds[1][1], bounds[2][1])));
    aabb->maxy = min(screen_height, max(bounds[0][1], max(bounds[1][1], bounds[2][1])));

    // set up edge functions (this is A = adj(M))
    v3_copy(E[0], (vec3){a0, b0, c0});
    v3_copy(E[1], (vec3){a1, b1, c1});
    v3_copy(E[2], (vec3){a2, b2, c2});

    v3_copy(Z, (vec3){(pos[0][2] - pos[2][2]), (pos[1][2] - pos[2][2]), pos[2][2]});

    return true;
    // return (detM > 0.0f);
}

static inline uint32_t _get_tile_index(uint32_t tile_x, uint32_t tile_y)
{
    return (tile_x + tile_y * RENDER_NUM_TILE_X);
}

static inline void _tile_index_to_coords(uint32_t index, uint32_t *tile_x, uint32_t *tile_y)
{
    *tile_y = index / RENDER_NUM_TILE_X;
    *tile_x = index % RENDER_NUM_TILE_X;
}

static void _binner(uint32_t tri_index, struct aabb *aabb, vec3 E[3])
{
    uint32_t min_tile_x = (uint32_t)(aabb->minx * RENDER_TILE_SIZE_INV);
    uint32_t min_tile_y = (uint32_t)(aabb->miny * RENDER_TILE_SIZE_INV);

    uint32_t max_tile_x = (uint32_t)ceilf(aabb->maxx * RENDER_TILE_SIZE_INV);
    uint32_t max_tile_y = (uint32_t)ceilf(aabb->maxy * RENDER_TILE_SIZE_INV);

    vec3 nrm_E0, nrm_E1, nrm_E2;
    v3_div(E[0], (fabsf(E[0][0]) + fabsf(E[0][1])), nrm_E0);
    v3_div(E[1], (fabsf(E[1][0]) + fabsf(E[1][1])), nrm_E1);
    v3_div(E[2], (fabsf(E[2][0]) + fabsf(E[2][1])), nrm_E2);

    static const vec2 tile_corner_offsets[] =
        {
            {0.0f, 0.0f},                        // LL (origin)
            {RENDER_TILE_SIZE, 0.0f},            // LR
            {0.0f, RENDER_TILE_SIZE},            // UL
            {RENDER_TILE_SIZE, RENDER_TILE_SIZE} // UR
        };

    const uint8_t TR_corner0 = (nrm_E0[1] >= 0.0f) ? ((nrm_E0[0] >= 0.0f) ? 3u : 2u)
                                                   : ((nrm_E0[0] >= 0.0f) ? 1u : 0u);
    const uint8_t TR_corner1 = (nrm_E1[1] >= 0.0f) ? ((nrm_E1[0] >= 0.0f) ? 3u : 2u)
                                                   : ((nrm_E1[0] >= 0.0f) ? 1u : 0u);
    const uint8_t TR_corner2 = (nrm_E2[1] >= 0.0f) ? ((nrm_E2[0] >= 0.0f) ? 3u : 2u)
                                                   : ((nrm_E2[0] >= 0.0f) ? 1u : 0u);

    const uint8_t TA_corner0 = 3u - TR_corner0;
    const uint8_t TA_corner1 = 3u - TR_corner1;
    const uint8_t TA_corner2 = 3u - TR_corner2;

    const float tile_pos_x = min(GRAFIKA_SCREEN_WIDTH, (float)(min_tile_x * RENDER_TILE_SIZE));
    const float tile_pos_y = min(GRAFIKA_SCREEN_HEIGHT, (float)(min_tile_y * RENDER_TILE_SIZE));

    float TR_edge_func_start_0 = ((nrm_E0[0] * tile_pos_x) + (nrm_E0[1] * tile_pos_y)) + nrm_E0[2];
    float TR_edge_func_start_1 = ((nrm_E1[0] * tile_pos_x) + (nrm_E1[1] * tile_pos_y)) + nrm_E1[2];
    float TR_edge_func_start_2 = ((nrm_E2[0] * tile_pos_x) + (nrm_E2[1] * tile_pos_y)) + nrm_E2[2];

    const float step0_x = (nrm_E0[0] * (float)(RENDER_TILE_SIZE)); // stepping 1 step
    const float step1_x = (nrm_E1[0] * (float)(RENDER_TILE_SIZE));
    const float step2_x = (nrm_E2[0] * (float)(RENDER_TILE_SIZE));

    const float step0_y = (nrm_E0[1] * (float)(RENDER_TILE_SIZE));
    const float step1_y = (nrm_E1[1] * (float)(RENDER_TILE_SIZE));
    const float step2_y = (nrm_E2[1] * (float)(RENDER_TILE_SIZE));

    const float _tr_contants0 = (nrm_E0[0] * (tile_corner_offsets[TR_corner0][0])) + (nrm_E0[1] * (tile_corner_offsets[TR_corner0][1]));
    const float _tr_contants1 = (nrm_E1[0] * (tile_corner_offsets[TR_corner1][0])) + (nrm_E1[1] * (tile_corner_offsets[TR_corner1][1]));
    const float _tr_contants2 = (nrm_E2[0] * (tile_corner_offsets[TR_corner2][0])) + (nrm_E2[1] * (tile_corner_offsets[TR_corner2][1]));

    const float _ta_contants0 = (nrm_E0[0] * (tile_corner_offsets[TA_corner0][0])) + (nrm_E0[1] * (tile_corner_offsets[TA_corner0][1]));
    const float _ta_contants1 = (nrm_E1[0] * (tile_corner_offsets[TA_corner1][0])) + (nrm_E1[1] * (tile_corner_offsets[TA_corner1][1]));
    const float _ta_contants2 = (nrm_E2[0] * (tile_corner_offsets[TA_corner2][0])) + (nrm_E2[1] * (tile_corner_offsets[TA_corner2][1]));

    for (uint32_t tile_y = min_tile_y, y = 0; tile_y < max_tile_y; tile_y++, y++)
    {
        const float y_step0 = step0_y * (float)y + TR_edge_func_start_0;
        const float y_step1 = step1_y * (float)y + TR_edge_func_start_1;
        const float y_step2 = step2_y * (float)y + TR_edge_func_start_2;

        float TR_edge_func_0 = _tr_contants0 + y_step0;
        float TR_edge_func_1 = _tr_contants1 + y_step1;
        float TR_edge_func_2 = _tr_contants2 + y_step2;

        float TA_edge_func_0 = _ta_contants0 + y_step0;
        float TA_edge_func_1 = _ta_contants1 + y_step1;
        float TA_edge_func_2 = _ta_contants2 + y_step2;

        for (uint32_t tile_x = min_tile_x; tile_x < max_tile_x; tile_x++,
                      TR_edge_func_0 += step0_x,
                      TR_edge_func_1 += step1_x,
                      TR_edge_func_2 += step2_x,

                      TA_edge_func_0 += step0_x,
                      TA_edge_func_1 += step1_x,
                      TA_edge_func_2 += step2_x)
        {
            if ((TR_edge_func_0 < 0.0f) || (TR_edge_func_1 < 0.0f) || (TR_edge_func_2 < 0.0f))
            {
                continue;
            }

            struct triangle_ref ref;
            ref.triangle_index = tri_index;

            uint32_t tile_index = _get_tile_index(tile_x, tile_y);

            if ((TA_edge_func_0 >= 0.0f) && (TA_edge_func_1 >= 0.0f) && (TA_edge_func_2 >= 0.0f))
            {
                //_draw_square_filled(tx * RENDER_TILE_SIZE, ty * RENDER_TILE_SIZE, RENDER_TILE_SIZE, 0xFF00FF00);
                ref.type = RASTER_FULL_TILE;
            }
            else
            {
                //_draw_square_outline(tx * RENDER_TILE_SIZE, ty * RENDER_TILE_SIZE, RENDER_TILE_SIZE, 0xFF00FF00);
                ref.type = RASTER_PARTIAL_TILE;
            }
            tile_add_data(tile_index, ref);
        }
    }
}

static void _shade_full_tile(uint32_t pos_x, uint32_t pos_y, vec3 E0, vec3 E1, vec3 E2, vec2 t0, vec2 t1, vec2 t2, vec3 Z)
{
    ASSERT(raster_state.tex.bpp == 4, "Texture bpp must be 4\n");

    const uint32_t start_x = pos_x * RENDER_TILE_SIZE;
    const uint32_t end_x   = start_x + RENDER_TILE_SIZE;

    const uint32_t start_y = pos_y * RENDER_TILE_SIZE;
    const uint32_t end_y   = start_y + RENDER_TILE_SIZE;

    const __m128 mm_E0_a = _mm_set1_ps(E0[0]);
    const __m128 mm_E1_a = _mm_set1_ps(E1[0]);
    const __m128 mm_E2_a = _mm_set1_ps(E2[0]);

    const __m128 mm_E0_b = _mm_set1_ps(E0[1]);
    const __m128 mm_E1_b = _mm_set1_ps(E1[1]);
    const __m128 mm_E2_b = _mm_set1_ps(E2[1]);

    const __m128 mm_E0_c = _mm_set1_ps(E0[2]);
    const __m128 mm_E1_c = _mm_set1_ps(E1[2]);
    const __m128 mm_E2_c = _mm_set1_ps(E2[2]);

    const __m128 mm_start_x = _mm_setr_ps((float)start_x + 0.5f,
                                          (float)start_x + 1.5f,
                                          (float)start_x + 2.5f,
                                          (float)start_x + 3.5f);

    __m128 mm_edge_start_0 = _mm_add_ps(_mm_add_ps(_mm_mul_ps(mm_E0_a, mm_start_x), _mm_mul_ps(mm_E0_b, _mm_set1_ps((float)start_y))), mm_E0_c);
    __m128 mm_edge_start_1 = _mm_add_ps(_mm_add_ps(_mm_mul_ps(mm_E1_a, mm_start_x), _mm_mul_ps(mm_E1_b, _mm_set1_ps((float)start_y))), mm_E1_c);
    __m128 mm_edge_start_2 = _mm_add_ps(_mm_add_ps(_mm_mul_ps(mm_E2_a, mm_start_x), _mm_mul_ps(mm_E2_b, _mm_set1_ps((float)start_y))), mm_E2_c);

    __m128       mm_4        = _mm_set1_ps(4.0f);
    const __m128 mm_step_x_0 = _mm_mul_ps(mm_E0_a, mm_4); // jump 4 pixels at a time
    const __m128 mm_step_x_1 = _mm_mul_ps(mm_E1_a, mm_4);
    const __m128 mm_step_x_2 = _mm_mul_ps(mm_E2_a, mm_4);

    const __m128 mm_step_y_0 = mm_E0_b;
    const __m128 mm_step_y_1 = mm_E1_b;
    const __m128 mm_step_y_2 = mm_E2_b;

    // pre fetch tex data
    const int      tex_w    = raster_state.tex.w;
    const int      tex_h    = raster_state.tex.h;
    unsigned char *tex_data = raster_state.tex.data;

    __m128 U[3], V[3], mm_Z[3];
    U[0] = _mm_set1_ps(t0[0]);
    U[1] = _mm_set1_ps(t1[0]);
    U[2] = _mm_set1_ps(t2[0]);

    V[0] = _mm_set1_ps(t0[1]);
    V[1] = _mm_set1_ps(t1[1]);
    V[2] = _mm_set1_ps(t2[1]);

    mm_Z[0] = _mm_set1_ps(Z[0]);
    mm_Z[1] = _mm_set1_ps(Z[1]);
    mm_Z[2] = _mm_set1_ps(Z[2]);

    float *depth_buffer = rend.depth_buffer;

    __m128 mm_tex_w = _mm_set1_ps((float)tex_w - 1);
    __m128 mm_tex_h = _mm_set1_ps((float)tex_h - 1);

    // end_x = end_x & ~3;       // round down end_x to the nearest multiple of 4
    // end_y = (end_y + 3) & ~3; // round up end_y to the next multiple of 4 - 1

    for (uint32_t py              = start_y; py < end_y; py++,
                  mm_edge_start_0 = _mm_add_ps(mm_edge_start_0, mm_step_y_0),
                  mm_edge_start_1 = _mm_add_ps(mm_edge_start_1, mm_step_y_1),
                  mm_edge_start_2 = _mm_add_ps(mm_edge_start_2, mm_step_y_2))
    {
        __m128 mm_edge_0 = mm_edge_start_0;
        __m128 mm_edge_1 = mm_edge_start_1;
        __m128 mm_edge_2 = mm_edge_start_2;

        for (uint32_t px        = start_x; px < end_x; px += 4,
                      mm_edge_0 = _mm_add_ps(mm_edge_0, mm_step_x_0),
                      mm_edge_1 = _mm_add_ps(mm_edge_1, mm_step_x_1),
                      mm_edge_2 = _mm_add_ps(mm_edge_2, mm_step_x_2))
        {
            // calculate interpolation coefficients
            __m128 F[3] = {0};
            F[0]        = mm_edge_0;
            F[1]        = mm_edge_1;
            F[2]        = mm_edge_2;

            __m128 r = _mm_add_ps(F[0], _mm_add_ps(F[1], F[2]));
            r        = _mm_rcp_ps(r);

            F[0] = _mm_mul_ps(F[0], r);
            F[1] = _mm_mul_ps(F[1], r);
            F[2] = _mm_mul_ps(F[2], r);

            const size_t pixel_index = (const size_t)(py * GRAFIKA_SCREEN_WIDTH + px);

            // interpolate depth value
            // float new_z = Z[0] * littlef_values[0] + Z[1] * littlef_values[1] + Z[2];
            const __m128 mm_pix_Z = _mm_add_ps(_mm_add_ps(_mm_mul_ps(mm_Z[0], F[0]),
                                                          _mm_mul_ps(mm_Z[1], F[1])),
                                               mm_Z[2]);

            float       *depth_location   = &depth_buffer[pixel_index];
            const __m128 mm_current_depth = _mm_load_ps(depth_location);

            const __m128 mm_depth_test_result = _mm_cmple_ps(mm_pix_Z, mm_current_depth);

            if ((uint16_t)_mm_movemask_ps(mm_depth_test_result) == 0x0) continue;

            const __m128 mm_write_mask = mm_depth_test_result;

            // write interpolated Z values
            _mm_maskmoveu_si128(
                _mm_castps_si128(mm_pix_Z),
                _mm_castps_si128(mm_write_mask),
                (char *)depth_location);

            // interpolate texture coordinates
            __m128 pixU = interpolate_values(F, U);
            __m128 pixV = interpolate_values(F, V);

            // clamp the vector to the range [0.0f, 1.0f]
            pixU = _mm_max_ps(_mm_min_ps(pixU, _mm_set1_ps(1.0f)), _mm_setzero_ps());
            pixV = _mm_max_ps(_mm_min_ps(pixV, _mm_set1_ps(1.0f)), _mm_setzero_ps());

            // TODO : simd the tex values
            pixU = _mm_mul_ps(pixU, mm_tex_w);
            pixV = _mm_mul_ps(pixV, mm_tex_h);

            __m128i mm_vw = _mm_mullo_epi32(_mm_cvtps_epi32(pixV), _mm_set1_epi32(tex_w));

            __m128i mm_withou_bpp = _mm_add_epi32(_mm_cvtps_epi32(pixU), mm_vw);
            __m128i tex_offsets   = _mm_slli_epi32(mm_withou_bpp, 2); // multiply each element by 4

            uint32_t *offset = (uint32_t *)&tex_offsets; // I think this is unsafe

            char *sample[4] = {0};
            sample[0]       = (char *)(tex_data + offset[3]);
            sample[1]       = (char *)(tex_data + offset[2]);
            sample[2]       = (char *)(tex_data + offset[1]);
            sample[3]       = (char *)(tex_data + offset[0]);

            __m128i final_colour = _mm_set_epi8(-1, sample[0][0], sample[0][1], sample[0][2],
                                                -1, sample[1][0], sample[1][1], sample[1][2],
                                                -1, sample[2][0], sample[2][1], sample[2][2],
                                                -1, sample[3][0], sample[3][1], sample[3][2]);

            uint32_t *pixel_location = &rend.pixels[pixel_index];

            const __m128i write_mask    = _mm_castps_si128(mm_write_mask);
            const __m128i masked_output = _mm_or_si128(_mm_and_si128(write_mask, final_colour),
                                                       _mm_andnot_si128(write_mask, *(__m128i *)pixel_location));

            _mm_store_si128((__m128i *)pixel_location, masked_output);
        }
    }
}

static void _shade_partial_tile(uint32_t pos_x, uint32_t pos_y, vec3 E0, vec3 E1, vec3 E2, vec2 t0, vec2 t1, vec2 t2, vec3 Z)
{
    ASSERT(raster_state.tex.bpp == 4, "Texture bpp must be 4\n");

    const uint32_t start_x = pos_x * RENDER_TILE_SIZE;
    const uint32_t end_x   = start_x + RENDER_TILE_SIZE;

    const uint32_t start_y = pos_y * RENDER_TILE_SIZE;
    const uint32_t end_y   = start_y + RENDER_TILE_SIZE;

    const __m128 mm_E0_a = _mm_set1_ps(E0[0]);
    const __m128 mm_E1_a = _mm_set1_ps(E1[0]);
    const __m128 mm_E2_a = _mm_set1_ps(E2[0]);

    const __m128 mm_E0_b = _mm_set1_ps(E0[1]);
    const __m128 mm_E1_b = _mm_set1_ps(E1[1]);
    const __m128 mm_E2_b = _mm_set1_ps(E2[1]);

    const __m128 mm_E0_c = _mm_set1_ps(E0[2]);
    const __m128 mm_E1_c = _mm_set1_ps(E1[2]);
    const __m128 mm_E2_c = _mm_set1_ps(E2[2]);

    const __m128 mm_start_x = _mm_setr_ps((float)start_x + 0.5f,
                                          (float)start_x + 1.5f,
                                          (float)start_x + 2.5f,
                                          (float)start_x + 3.5f);

    __m128 mm_edge_start_0 = _mm_add_ps(_mm_add_ps(_mm_mul_ps(mm_E0_a, mm_start_x), _mm_mul_ps(mm_E0_b, _mm_set1_ps((float)start_y))), mm_E0_c);
    __m128 mm_edge_start_1 = _mm_add_ps(_mm_add_ps(_mm_mul_ps(mm_E1_a, mm_start_x), _mm_mul_ps(mm_E1_b, _mm_set1_ps((float)start_y))), mm_E1_c);
    __m128 mm_edge_start_2 = _mm_add_ps(_mm_add_ps(_mm_mul_ps(mm_E2_a, mm_start_x), _mm_mul_ps(mm_E2_b, _mm_set1_ps((float)start_y))), mm_E2_c);

    __m128       mm_4        = _mm_set1_ps(4.0f);
    const __m128 mm_step_x_0 = _mm_mul_ps(mm_E0_a, mm_4); // jump 4 pixels at a time
    const __m128 mm_step_x_1 = _mm_mul_ps(mm_E1_a, mm_4);
    const __m128 mm_step_x_2 = _mm_mul_ps(mm_E2_a, mm_4);

    const __m128 mm_step_y_0 = mm_E0_b;
    const __m128 mm_step_y_1 = mm_E1_b;
    const __m128 mm_step_y_2 = mm_E2_b;

    const __m128 mm_edge0_AB_tie_breaker_result = _mm_or_ps(_mm_cmpgt_ps(mm_E0_a, _mm_setzero_ps()),
                                                            _mm_and_ps(_mm_cmpeq_ps(mm_E0_a, _mm_setzero_ps()),
                                                                       _mm_cmpge_ps(mm_E0_b, _mm_setzero_ps())));

    const __m128 mm_edge1_AB_tie_breaker_result = _mm_or_ps(_mm_cmpgt_ps(mm_E1_a, _mm_setzero_ps()),
                                                            _mm_and_ps(_mm_cmpeq_ps(mm_E1_a, _mm_setzero_ps()),
                                                                       _mm_cmpge_ps(mm_E1_b, _mm_setzero_ps())));

    const __m128 mm_edge2_AB_tie_breaker_result = _mm_or_ps(_mm_cmpgt_ps(mm_E2_a, _mm_setzero_ps()),
                                                            _mm_and_ps(_mm_cmpeq_ps(mm_E2_a, _mm_setzero_ps()),
                                                                       _mm_cmpge_ps(mm_E2_b, _mm_setzero_ps())));

    // pre fetch tex data
    const int      tex_w    = raster_state.tex.w;
    const int      tex_h    = raster_state.tex.h;
    unsigned char *tex_data = raster_state.tex.data;

    __m128 U[3], V[3], mm_Z[3];
    U[0] = _mm_set1_ps(t0[0]);
    U[1] = _mm_set1_ps(t1[0]);
    U[2] = _mm_set1_ps(t2[0]);

    V[0] = _mm_set1_ps(t0[1]);
    V[1] = _mm_set1_ps(t1[1]);
    V[2] = _mm_set1_ps(t2[1]);

    mm_Z[0] = _mm_set1_ps(Z[0]);
    mm_Z[1] = _mm_set1_ps(Z[1]);
    mm_Z[2] = _mm_set1_ps(Z[2]);

    float *depth_buffer = rend.depth_buffer;

    __m128 mm_tex_w = _mm_set1_ps((float)tex_w - 1);
    __m128 mm_tex_h = _mm_set1_ps((float)tex_h - 1);

    // end_x = end_x & ~3;       // round down end_x to the nearest multiple of 4
    // end_y = (end_y + 3) & ~3; // round up end_y to the next multiple of 4 - 1

    for (uint32_t py              = start_y; py < end_y; py++,
                  mm_edge_start_0 = _mm_add_ps(mm_edge_start_0, mm_step_y_0),
                  mm_edge_start_1 = _mm_add_ps(mm_edge_start_1, mm_step_y_1),
                  mm_edge_start_2 = _mm_add_ps(mm_edge_start_2, mm_step_y_2))
    {
        __m128 mm_edge_0 = mm_edge_start_0;
        __m128 mm_edge_1 = mm_edge_start_1;
        __m128 mm_edge_2 = mm_edge_start_2;

        for (uint32_t px        = start_x; px < end_x; px += 4,
                      mm_edge_0 = _mm_add_ps(mm_edge_0, mm_step_x_0),
                      mm_edge_1 = _mm_add_ps(mm_edge_1, mm_step_x_1),
                      mm_edge_2 = _mm_add_ps(mm_edge_2, mm_step_x_2))
        {

            //     (edge_value < 0.0f) || (edge_value == 0.0f && ab_test_result);
            __m128 _edge_positive      = _mm_cmpgt_ps(mm_edge_0, _mm_setzero_ps());
            __m128 _edge_negative      = _mm_cmplt_ps(mm_edge_0, _mm_setzero_ps());
            __m128 mm_edge_0_test_mask = _mm_or_ps(_edge_positive, _mm_andnot_ps(_edge_negative, mm_edge0_AB_tie_breaker_result));

            _edge_positive             = _mm_cmpgt_ps(mm_edge_1, _mm_setzero_ps());
            _edge_negative             = _mm_cmplt_ps(mm_edge_1, _mm_setzero_ps());
            __m128 mm_edge_1_test_mask = _mm_or_ps(_edge_positive, _mm_andnot_ps(_edge_negative, mm_edge1_AB_tie_breaker_result));

            _edge_positive             = _mm_cmpgt_ps(mm_edge_2, _mm_setzero_ps());
            _edge_negative             = _mm_cmplt_ps(mm_edge_2, _mm_setzero_ps());
            __m128 mm_edge_2_test_mask = _mm_or_ps(_edge_positive, _mm_andnot_ps(_edge_negative, mm_edge2_AB_tie_breaker_result));

            //__m128 mm_edge_tie_breaker = _mm_or_ps(mm_edge_0_test_mask, _mm_or_ps(mm_edge_1_test_mask, mm_edge_2_test_mask));
            __m128 mm_edge_tie_breaker = _mm_and_ps(mm_edge_0_test_mask, _mm_and_ps(mm_edge_1_test_mask, mm_edge_2_test_mask));

            uint16_t tie_breaker_result = (uint16_t)_mm_movemask_ps(mm_edge_tie_breaker);

            if (tie_breaker_result == 0x0) continue;

            // calculate interpolation coefficients
            __m128 F[3] = {0};
            F[0]        = mm_edge_0;
            F[1]        = mm_edge_1;
            F[2]        = mm_edge_2;

            // R(x, y) = F0(x, y) + F1(x, y) + F2(x, y)
            // r = 1 / (F0(x, y) + F1(x, y) + F2(x, y))
            __m128 r = _mm_add_ps(F[0], _mm_add_ps(F[1], F[2]));
            r        = _mm_rcp_ps(r);

            F[0] = _mm_mul_ps(F[0], r);
            F[1] = _mm_mul_ps(F[1], r);
            F[2] = _mm_mul_ps(F[2], r);

            const size_t pixel_index = (const size_t)(py * GRAFIKA_SCREEN_WIDTH + px);

            // interpolate depth value
            // float new_z = Z[0] * littlef_values[0] + Z[1] * littlef_values[1] + Z[2];
            const __m128 mm_pix_Z = _mm_add_ps(_mm_add_ps(_mm_mul_ps(mm_Z[0], F[0]),
                                                          _mm_mul_ps(mm_Z[1], F[1])),
                                               mm_Z[2]);

            float       *depth_location   = &depth_buffer[pixel_index];
            const __m128 mm_current_depth = _mm_load_ps(depth_location);

            const __m128 mm_depth_test_result = _mm_cmple_ps(mm_pix_Z, mm_current_depth);

            if ((uint16_t)_mm_movemask_ps(mm_depth_test_result) == 0x0) continue;

            // AND depth mask & coverage mask for quads of fragments
            const __m128 mm_write_mask = _mm_and_ps(mm_depth_test_result, mm_edge_tie_breaker);

            // write interpolated Z values
            _mm_maskmoveu_si128(
                _mm_castps_si128(mm_pix_Z),
                _mm_castps_si128(mm_write_mask),
                (char *)depth_location);

            // interpolate texture coordinates
            __m128 pixU = interpolate_values(F, U);
            __m128 pixV = interpolate_values(F, V);

            // clamp the vector to the range [0.0f, 1.0f]
            pixU = _mm_max_ps(_mm_min_ps(pixU, _mm_set1_ps(1.0f)), _mm_setzero_ps());
            pixV = _mm_max_ps(_mm_min_ps(pixV, _mm_set1_ps(1.0f)), _mm_setzero_ps());

            pixU = _mm_mul_ps(pixU, mm_tex_w);
            pixV = _mm_mul_ps(pixV, mm_tex_h);

            __m128i mm_vw = _mm_mullo_epi32(_mm_cvtps_epi32(pixV), _mm_set1_epi32(tex_w));

            __m128i mm_withou_bpp = _mm_add_epi32(_mm_cvtps_epi32(pixU), mm_vw);
            __m128i tex_offsets   = _mm_slli_epi32(mm_withou_bpp, 2); // multiply each element by 4

            uint32_t *offset = (uint32_t *)&tex_offsets; // I think this is unsafe

            char *sample[4] = {0};
            sample[0]       = (char *)(tex_data + offset[3]);
            sample[1]       = (char *)(tex_data + offset[2]);
            sample[2]       = (char *)(tex_data + offset[1]);
            sample[3]       = (char *)(tex_data + offset[0]);

            __m128i final_colour = _mm_set_epi8(-1, sample[0][0], sample[0][1], sample[0][2],
                                                -1, sample[1][0], sample[1][1], sample[1][2],
                                                -1, sample[2][0], sample[2][1], sample[2][2],
                                                -1, sample[3][0], sample[3][1], sample[3][2]);

            uint32_t *pixel_location = &rend.pixels[pixel_index];
#if 0
                const __m128i write_mask    = _mm_castps_si128(mm_write_mask);
                const __m128i masked_output = _mm_or_si128(_mm_and_si128(write_mask, final_colour),
                _mm_andnot_si128(write_mask, *(__m128i *)pixel_location));
                
                _mm_store_si128((__m128i *)pixel_location, masked_output);
#else
            _mm_maskstore_epi32(
                (int *)pixel_location,
                _mm_castps_si128(mm_write_mask),
                final_colour);
#endif
        }
    }
}

struct triangle_data
{
    vec4 pos[3];
    vec2 tex[3];
    vec3 E[3];
    vec3 Z;
};

struct raster_job_data
{
    struct triangle_data *triangles;
    uint32_t              tile_lookup_index;
};

static void raster_job(void *data)
{
    struct raster_job_data *d = (struct raster_job_data *)data;
    ASSERT(d, "Raster job data is inavlid\n");

    struct triangle_data *triangles = d->triangles;

    struct tile        *t_data   = &tile_array[d->tile_lookup_index];
    struct tile_coords *t_coords = &tile_coords[d->tile_lookup_index];

    uint32_t tile_x = t_coords->tile_x;
    uint32_t tile_y = t_coords->tile_y;

    for (uint32_t b = 0; b < t_data->triangle_count; b++)
    {
        struct triangle_ref t = t_data->triangles[b];

        struct triangle_data *td = triangles + t.triangle_index;

        switch (t.type)
        {
            case RASTER_FULL_TILE:
            {
                _shade_full_tile(tile_x, tile_y, td->E[0], td->E[1], td->E[2], td->tex[0], td->tex[1], td->tex[2], td->Z);
                break;
            }
            case RASTER_PARTIAL_TILE:
            {
                _shade_partial_tile(tile_x, tile_y, td->E[0], td->E[1], td->E[2], td->tex[0], td->tex[1], td->tex[2], td->Z);
                break;
            }
            default: ASSERT(false, "Invalid tile rendering\n");
        }
    }
}

void draw_onstart(struct arena *arena)
{
    UNUSED(arena);

    struct obj obj = raster_state.obj;

    ASSERT(obj.mats, "Object must have a material\n");
    ASSERT(obj.mats[0].map_Kd, "Object must have a diffuse texture\n");

    raster_state.tex = tex_load(obj.mats[0].map_Kd);

    LOG("Processor count : %u\n", jobs_get_proccessor_count());

    jobs_init();

    for (uint32_t i = 0; i < RENDER_TILE_COUNT; i++)
    {
        uint32_t tile_x, tile_y;
        _tile_index_to_coords(i, &tile_x, &tile_y);

        tile_coords[i].tile_x = cast_u32_to_u8(tile_x);
        tile_coords[i].tile_y = cast_u32_to_u8(tile_y);
    }
}

void draw_object(struct arena *arena)
{
    struct arena tmp_arena = arena_create_tmp(arena);

    struct obj          obj         = raster_state.obj;
    float              *obj_pos     = obj.pos;
    float              *obj_tex     = obj.texs;
    struct vertindices *obj_indices = obj.indices;

    // pass triangles through our "Vertex Shader"
    float *transformed_vertices = arena_alloc_aligned(&tmp_arena, sizeof(vec4) * obj.num_pos, 16);

    // TODO : multi thread this loop
    for (size_t i = 0; i < obj.num_pos; i++)
    {
        vec4 vertex = {obj_pos[i * 3 + 0],
                       obj_pos[i * 3 + 1],
                       obj_pos[i * 3 + 2],
                       1.0f};

        // vertex pos * MVP
        m4_mul_v4(raster_state.MVP, vertex, &transformed_vertices[i * 4]);
    }

    // now are verts are in CLIP SPACE
    const size_t          triangle_count = obj.num_f_rows;
    struct triangle_data *triangles      = ma_push_size(&tmp_arena, sizeof(struct triangle_data) * triangle_count);

    // TODO : multi thread this loop
    for (size_t i = 0; i < triangle_count; ++i)
    {
        for (size_t j = 0; j < 3; ++j)
        {
            struct vertindices indices = obj_indices[i * 3 + j];

            const int pos_index = indices.v_idx;
            const int tex_index = indices.vt_idx;

#ifdef CCW_TRIANGLES
            const size_t store_index = j; // CCW triangles
#else
            const size_t store_index = (2 - j); // CW triangles (which blender uses)
#endif

            float *p = transformed_vertices + 4 * pos_index;
            float *t = obj_tex + 2 * tex_index;

            triangles[i].pos[store_index][0] = p[0];
            triangles[i].pos[store_index][1] = p[1];
            triangles[i].pos[store_index][2] = p[2];
            triangles[i].pos[store_index][3] = p[3];

            // TODO : why are we doing this every draw call?
            triangles[i].tex[store_index][0] = t[0];
            triangles[i].tex[store_index][1] = t[1];
        }
        // TODO : we have a triangle, why not bin it now?
    }

    // clean out the bins
    memset(tile_array, 0, sizeof(tile_array));
    memset(tiles_with_bins, -1, sizeof(tiles_with_bins));
    tiles_queue_index = 0;

    // now we have all the triangles, in CLIP SPACE
    // TODO : multi thread this loop
    for (size_t i = 0; i < triangle_count; i++)
    {
        struct triangle_data *t = triangles + i;

        struct aabb aabb;

        // CLIPPER
        if (!_triangle_clipping(t->pos))
        {
            continue; // triangle clipped, move to next triangle
        }

        // TRIANGLE SETUP & CULL
        if (!_triangle_setup_and_cull(t->pos, t->E, t->Z, &aabb))
        {
            continue; // triangle culled, move to next triangle
        }

        // BINNER
        _binner((uint32_t)i, &aabb, t->E);
    }

    // NOTE : here we have collected all "tiles" are in contact with a triangle,
    //          we could go a step further and basically find all pixels that are
    //          touching a triangle, but for now, I am not doing this

    for (uint32_t i = 0; i < tiles_queue_index; i++)
    {
        uint32_t tile_lookup_index = tiles_with_bins[i];
        // if (tile_lookup_index > RENDER_TILE_COUNT) continue;

        ASSERT(tile_lookup_index < RENDER_TILE_COUNT, "Lookup tile index is out of range\n");

        struct raster_job_data *rd = ma_push_struct(&tmp_arena, struct raster_job_data);

        // NOTE : tile_array are shared between threads, doesnt change, so can store index instead
        // NOTE : triangles are shared between threads, doesnt change
        rd->tile_lookup_index = tile_lookup_index;
        rd->triangles         = triangles;

        job_submit((struct job){.arguments = rd, .function = raster_job});
    }

    jobs_complete_all_work();
}

void draw_onexit(void)
{
    tex_destroy(&raster_state.tex);
    // raster_state.obj       // cleaned up by arena
}
