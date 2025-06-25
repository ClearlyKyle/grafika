// https://tayfunkayhan.wordpress.com/2019/07/26/chasing-triangles-in-a-tile-based-rasterizer/

// 1 - Triangles are evenly spread between threads
// 2 - Each triangle is checked againts which "tile" it interacts with
//      save data for the triangle that is in a "tile"
// 3 - Now each tile knows what triangle it interacts with we can raster each tile
// 4 - we can raster a triangle in small tiles, think 8x8 for AVX SIMD

#include "common.h"

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

struct aabb
{
    float minx, miny, maxx, maxy;
};

enum raster_type
{
    RASTER_FULL_TILE,
    RASTER_PARTIAL_TILE,
    RASTER_QUAD,
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
    bool                has_been_queued;
};

struct tile tile_array[RENDER_TILE_COUNT];
uint32_t    tiles_with_bins[RENDER_TILE_COUNT];
uint32_t    tiles_queue_index = 0;

struct binning
{
    mat4       model;
    mat4       MVP;
    vec3       cam_pos;
    struct tex tex;
    struct obj obj;
};

struct binning raster_state = {0};

static inline bool tie_breaker_ab_test(const vec3 E)
{
    return (E[0] > 0.0f) || (E[0] == 0.0f && E[1] >= 0.0f);
}

static inline bool edge_tie_breaker(const float edge_value, const bool ab_test_result)
{
    // NOTE : cation with this less than direction and triangle windings!!
    //          and with the matrix a, b, c value setup!!
    return (edge_value < 0.0f) || (edge_value == 0.0f && ab_test_result);
}

static inline float interpolate_values(vec3 little_f_values, vec3 attribute)
{
    return v3_dot(little_f_values, attribute);
}

static inline void tile_add_data(uint32_t tile_index, struct triangle_ref ref)
{
    struct tile *tile = &tile_array[tile_index];

#pragma omp atomic capture
    uint32_t index = tile->triangle_count++;

    ASSERT(index < MAX_TRIANGLES_PER_TILE, "Too many triangles in this bin\n");

    tile->triangles[index] = ref;

    // prevent double-queueing tiles
#pragma omp critical(tiles_queue_lock)
    {
        if (!tile->has_been_queued)
        {
            tile->has_been_queued                = true;
            tiles_with_bins[tiles_queue_index++] = tile_index;
        }
    }
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
        // homo[i][1] = (pos[i][3] + pos[i][1]) * GRAFIKA_SCREEN_HALF_HEIGHT; // origin Bottom Left
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
}

static inline uint32_t _get_tile_index(uint32_t tile_x, uint32_t tile_y)
{
    return (tile_x + tile_y * RENDER_NUM_TILE_X);
}

// TODO : maybe store these somewhere?
static inline void _tile_index_to_coords(uint32_t index, uint32_t *tile_x, uint32_t *tile_y)
{
    *tile_y = index / RENDER_NUM_TILE_X;
    *tile_x = index % RENDER_NUM_TILE_X;
}

static void draw_aabb_outline(struct aabb *box)
{
    for (int x = (int)box->minx; x <= (int)box->maxx; ++x)
    {
        grafika_setpixel((uint32_t)x, (uint32_t)box->miny, 0xFF00FF00); // Top
        grafika_setpixel((uint32_t)x, (uint32_t)box->maxy, 0xFF00FF00); // Bottom
    }

    for (int y = (int)box->miny; y <= (int)box->maxy; ++y)
    {
        grafika_setpixel((uint32_t)box->minx, (uint32_t)y, 0xFF00FF00); // Left
        grafika_setpixel((uint32_t)box->maxx, (uint32_t)y, 0xFF00FF00); // Right
    }
}

static void _draw_square_outline(uint32_t x, uint32_t y, uint32_t length, uint32_t colour)
{
    for (uint32_t i = 0; i < length; i++)
    {
        grafika_setpixel(x + i, y, colour);              // Top edge
        grafika_setpixel(x + i, y + length - 1, colour); // Bottom edge
        grafika_setpixel(x, y + i, colour);              // Left edge
        grafika_setpixel(x + length - 1, y + i, colour); // Right edge
    }
}

static void _draw_square_filled(uint32_t x, uint32_t y, uint32_t length, uint32_t colour)
{
    for (uint32_t j = 0; j < length; j++)
    {
        for (uint32_t i = 0; i < length; i++)
        {
            grafika_setpixel(x + i, y + j, colour);
        }
    }
}

static inline uint8_t _compute_corner_index(const float x, const float y)
{
    const uint8_t x_bit = (x >= 0.0f) << 1;
    const uint8_t y_bit = (y >= 0.0f) << 0;
    return (y_bit | x_bit) ^ (y_bit >> 1); // returns : 0,1,2,3
}

#if 1
// stepped binner
static void _binner(uint32_t tri_index, struct aabb *aabb, vec3 E[3])
{
    ASSERT(IS_POWER_OF_2(RENDER_TILE_SIZE), "RENDER_TILE_SIZE must be a power of 2!\n");

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

    // (TA) Trivial-Accept: tile is inside of triangles bounding box
    // (TR) Trivial-Reject: tile is outside of triangles bounding box
    // (O)  Overlap: tile overlaps the triangles bounding box

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
                //_draw_square_outline(tile_x * RENDER_TILE_SIZE, tile_y * RENDER_TILE_SIZE, RENDER_TILE_SIZE, 0x00FFFFFF);
                ref.type = RASTER_FULL_TILE;
            }
            else
            {
                //_draw_square_outline(tile_x * RENDER_TILE_SIZE, tile_y * RENDER_TILE_SIZE, RENDER_TILE_SIZE, 0xFF00FF00);
                ref.type = RASTER_PARTIAL_TILE;
            }
            tile_add_data(tile_index, ref);
        }
    }
}
#else
static void _binner(uint32_t tri_index, struct aabb *aabb, vec3 E[3])
{
    ASSERT((aabb->minx >= 0.0f) && (aabb->maxx >= 0.0f) && (aabb->miny >= 0.0f) && (aabb->maxy >= 0.0f), "clipper must have clamped aabb to screen extents!");
    ASSERT((aabb->minx <= aabb->maxx) && (aabb->miny <= aabb->maxy), "clipper must have clamped aabb to screen extents!");

    // given the tile size and window dimensions, we compute the min/max tile range
    // that intersects the triangles AABB

    ASSERT(IS_POWER_OF_2(RENDER_TILE_SIZE), "RENDER_TILE_SIZE must be a power of 2!\n");

    uint32_t min_tile_x = (uint32_t)(aabb->minx * RENDER_TILE_SIZE_INV);
    uint32_t min_tile_y = (uint32_t)(aabb->miny * RENDER_TILE_SIZE_INV);

    uint32_t max_tile_x = (uint32_t)ceilf(aabb->maxx * RENDER_TILE_SIZE_INV);
    uint32_t max_tile_y = (uint32_t)ceilf(aabb->maxy * RENDER_TILE_SIZE_INV);

    ASSERT((min_tile_x <= max_tile_x) && (max_tile_x <= RENDER_NUM_TILE_X), "tile_x error\n");
    ASSERT((min_tile_y <= max_tile_y) && (max_tile_y <= RENDER_NUM_TILE_Y), "tile_y error\n");

    // normalize edge functions
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

    // (x, y) -> sample location | (a, b, c) -> edge equation coefficients
    // Edge function ->             E(x, y) = (a * x) + (b * y) + c
    // Stepping edge function ->    E(x + s, y + t) = E(x, y) + (a * s) + (b * t)

    // (TA) Trivial-Accept: tile is inside of triangles bounding box
    // (TR) Trivial-Reject: tile is outside of triangles bounding box
    // (O)  Overlap: tile overlaps the triangles bounding box

    const uint8_t TR_corner0 = (nrm_E0[1] >= 0.0f) ? ((nrm_E0[0] >= 0.0f) ? 3u : 2u)
                                                   : ((nrm_E0[0] >= 0.0f) ? 1u : 0u);
    const uint8_t TR_corner1 = (nrm_E1[1] >= 0.0f) ? ((nrm_E1[0] >= 0.0f) ? 3u : 2u)
                                                   : ((nrm_E1[0] >= 0.0f) ? 1u : 0u);
    const uint8_t TR_corner2 = (nrm_E2[1] >= 0.0f) ? ((nrm_E2[0] >= 0.0f) ? 3u : 2u)
                                                   : ((nrm_E2[0] >= 0.0f) ? 1u : 0u);

    const uint8_t TA_corner0 = 3u - TR_corner0;
    const uint8_t TA_corner1 = 3u - TR_corner1;
    const uint8_t TA_corner2 = 3u - TR_corner2;

    // evaluate edge function for the first tile incontact with the AABB region, then
    // step it over the entire AABB
    const float tile_pos_x = min(GRAFIKA_SCREEN_WIDTH, (float)(min_tile_x * RENDER_TILE_SIZE));
    const float tile_pos_y = min(GRAFIKA_SCREEN_HEIGHT, (float)(min_tile_y * RENDER_TILE_SIZE));

    float edge_func_start_0 = ((nrm_E0[0] * tile_pos_x) + (nrm_E0[1] * tile_pos_y)) + nrm_E0[2];
    float edge_func_start_1 = ((nrm_E1[0] * tile_pos_x) + (nrm_E1[1] * tile_pos_y)) + nrm_E1[2];
    float edge_func_start_2 = ((nrm_E2[0] * tile_pos_x) + (nrm_E2[1] * tile_pos_y)) + nrm_E2[2];

    // (tile_y, tile_x) tile index
    // (x, y) pixel location relative to tile origin
    for (uint32_t tile_y = min_tile_y, y = 0; tile_y < max_tile_y; tile_y++, y++)
    {
        for (uint32_t tile_x = min_tile_x, x = 0; tile_x < max_tile_x; tile_x++, x++)
        {
            //_draw_square_outline(tile_x * RENDER_TILE_SIZE, tile_y * RENDER_TILE_SIZE, RENDER_TILE_SIZE, 0xFFFFFFFF);

            const float x_step = (float)(x * RENDER_TILE_SIZE);
            const float y_step = (float)(y * RENDER_TILE_SIZE);

            float TR_edge_func_0 = edge_func_start_0 + ((nrm_E0[0] * (tile_corner_offsets[TR_corner0][0] + x_step)) + (nrm_E0[1] * (tile_corner_offsets[TR_corner0][1] + y_step)));
            float TR_edge_func_1 = edge_func_start_1 + ((nrm_E1[0] * (tile_corner_offsets[TR_corner1][0] + x_step)) + (nrm_E1[1] * (tile_corner_offsets[TR_corner1][1] + y_step)));
            float TR_edge_func_2 = edge_func_start_2 + ((nrm_E2[0] * (tile_corner_offsets[TR_corner2][0] + x_step)) + (nrm_E2[1] * (tile_corner_offsets[TR_corner2][1] + y_step)));

            if ((TR_edge_func_0 < 0.0f) || (TR_edge_func_1 < 0.0f) || (TR_edge_func_2 < 0.0f))
            {
                // Trivial Rejected - tile is completely outside of the triangle
                // LOG("TR Tile %u\n", _get_tile_index(tx, ty));
                //_draw_square_outline(x * RENDER_TILE_SIZE, y * RENDER_TILE_SIZE, RENDER_TILE_SIZE, 0xFF0000FF);
                continue;
            }

            // tile is either partially or completely inside of the triangle
            struct triangle_ref ref;
            ref.triangle_index = tri_index;

            uint32_t tile_index = _get_tile_index(tile_x, tile_y); // NOTE : 2D tile (x, y), into 1D array value

            float TA_edge_func_0 = edge_func_start_0 + ((nrm_E0[0] * (tile_corner_offsets[TA_corner0][0] + x_step)) + (nrm_E0[1] * (tile_corner_offsets[TA_corner0][1] + y_step)));
            float TA_edge_func_1 = edge_func_start_1 + ((nrm_E1[0] * (tile_corner_offsets[TA_corner1][0] + x_step)) + (nrm_E1[1] * (tile_corner_offsets[TA_corner1][1] + y_step)));
            float TA_edge_func_2 = edge_func_start_2 + ((nrm_E2[0] * (tile_corner_offsets[TA_corner2][0] + x_step)) + (nrm_E2[1] * (tile_corner_offsets[TA_corner2][1] + y_step)));

            if ((TA_edge_func_0 >= 0.0f) && (TA_edge_func_1 >= 0.0f) && (TA_edge_func_2 >= 0.0f))
            {
                // Trivial Accepted - tile is completely inside of the triangle

                // LOG("TA Tile %u\n", _get_tile_index(tx, ty));
                //_draw_square_outline(x * RENDER_TILE_SIZE, y * RENDER_TILE_SIZE, RENDER_TILE_SIZE, 0xFFFF00FF);
                ref.type = RASTER_FULL_TILE;
            }
            else
            {
                // Overlap - tile is partially covered by the triangle, bin the triangle for the tile
                // LOG("O Tile %u\n", _get_tile_index(tx, ty));
                //_draw_square_outline(x * RENDER_TILE_SIZE, y * RENDER_TILE_SIZE, RENDER_TILE_SIZE, 0xFF0000FF);
                ref.type = RASTER_PARTIAL_TILE;
            }
            tile_add_data(tile_index, ref);
        }
    }
}
#endif

static void _shade_plain_full_tile(uint32_t pos_x, uint32_t pos_y, vec3 E0, vec3 E1, vec3 E2, vec2 t0, vec2 t1, vec2 t2, vec3 Z)
{
    const uint32_t start_x = pos_x * RENDER_TILE_SIZE;
    const uint32_t end_x   = start_x + RENDER_TILE_SIZE;

    const uint32_t start_y = pos_y * RENDER_TILE_SIZE;
    const uint32_t end_y   = start_y + RENDER_TILE_SIZE;

    const float step0_x = E0[0];
    const float step1_x = E1[0];
    const float step2_x = E2[0];

    const float step0_y = E0[1];
    const float step1_y = E1[1];
    const float step2_y = E2[1];

    float edge_func_start_0 = (E0[0] * (float)start_x) + (E0[1] * (float)start_y) + E0[2];
    float edge_func_start_1 = (E1[0] * (float)start_x) + (E1[1] * (float)start_y) + E1[2];
    float edge_func_start_2 = (E2[0] * (float)start_x) + (E2[1] * (float)start_y) + E2[2];

    // pre fetch tex data
    const int      tex_w    = raster_state.tex.w;
    const int      tex_h    = raster_state.tex.h;
    const int      tex_bpp  = raster_state.tex.bpp;
    unsigned char *tex_data = raster_state.tex.data;

    vec3 u_values = {t0[0], t1[0], t2[0]};
    vec3 v_values = {t0[1], t1[1], t2[1]};

    float *depth_buffer = rend.depth_buffer;

    for (uint32_t py = start_y; py < end_y; py++,
                  edge_func_start_0 += step0_y,
                  edge_func_start_1 += step1_y,
                  edge_func_start_2 += step2_y

    )
    {
        float edge_func_0 = edge_func_start_0;
        float edge_func_1 = edge_func_start_1;
        float edge_func_2 = edge_func_start_2;

        for (uint32_t px = start_x; px < end_x; px++,
                      edge_func_0 += step0_x,
                      edge_func_1 += step1_x,
                      edge_func_2 += step2_x

        )
        {
            // calculate interpolation coefficients
            const float F0 = edge_func_0;
            const float F1 = edge_func_1;
            const float F2 = edge_func_2;

            const float r = 1.0f / (F0 + F1 + F2);

            vec3 littlef_values = {F0 * r, F1 * r, F2 * r};

            // if we have deltas for out Z values, use this form, else, interpolate_values
            float new_z = Z[0] * littlef_values[0] + Z[1] * littlef_values[1] + Z[2];
            // float  new_z = interpolate_values(littlef_values, Z);

            float *old_z = depth_buffer + ((GRAFIKA_SCREEN_WIDTH * py) + px);

            // depth test
            if (new_z >= *old_z) continue;
            *old_z = new_z;

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

            grafika_setpixel(px, py, pixel_colour);
        }
    }
}

static void _shade_full_tile(uint32_t pos_x, uint32_t pos_y, vec3 E0, vec3 E1, vec3 E2, vec2 t0, vec2 t1, vec2 t2, vec3 Z)
{
    const uint32_t start_x = pos_x * RENDER_TILE_SIZE;
    const uint32_t end_x   = start_x + RENDER_TILE_SIZE;

    const uint32_t start_y = pos_y * RENDER_TILE_SIZE;
    const uint32_t end_y   = start_y + RENDER_TILE_SIZE;

    const float step0_x = E0[0];
    const float step1_x = E1[0];
    const float step2_x = E2[0];

    const float step0_y = E0[1];
    const float step1_y = E1[1];
    const float step2_y = E2[1];

    float edge_func_start_0 = (E0[0] * (float)start_x) + (E0[1] * (float)start_y) + E0[2];
    float edge_func_start_1 = (E1[0] * (float)start_x) + (E1[1] * (float)start_y) + E1[2];
    float edge_func_start_2 = (E2[0] * (float)start_x) + (E2[1] * (float)start_y) + E2[2];

    // pre fetch tex data
    const int      tex_w    = raster_state.tex.w;
    const int      tex_h    = raster_state.tex.h;
    const int      tex_bpp  = raster_state.tex.bpp;
    unsigned char *tex_data = raster_state.tex.data;

    vec3 u_values = {t0[0], t1[0], t2[0]};
    vec3 v_values = {t0[1], t1[1], t2[1]};

    float *depth_buffer = rend.depth_buffer;

    for (uint32_t py = start_y; py < end_y; py++,
                  edge_func_start_0 += step0_y,
                  edge_func_start_1 += step1_y,
                  edge_func_start_2 += step2_y

    )
    {
        float edge_func_0 = edge_func_start_0;
        float edge_func_1 = edge_func_start_1;
        float edge_func_2 = edge_func_start_2;

        for (uint32_t px = start_x; px < end_x; px++,
                      edge_func_0 += step0_x,
                      edge_func_1 += step1_x,
                      edge_func_2 += step2_x

        )
        {
            // calculate interpolation coefficients
            const float F0 = edge_func_0;
            const float F1 = edge_func_1;
            const float F2 = edge_func_2;

            const float r = 1.0f / (F0 + F1 + F2);

            vec3 littlef_values = {F0 * r, F1 * r, F2 * r};

            // if we have deltas for out Z values, use this form, else, interpolate_values
            float new_z = Z[0] * littlef_values[0] + Z[1] * littlef_values[1] + Z[2];
            // float  new_z = interpolate_values(littlef_values, Z);

            float *old_z = depth_buffer + ((GRAFIKA_SCREEN_WIDTH * py) + px);

            // depth test
            if (new_z > *old_z) continue;
            *old_z = new_z;

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

            grafika_setpixel(px, py, pixel_colour);
        }
    }
}

static void _shade_partial_tile(uint32_t pos_x, uint32_t pos_y, vec3 E0, vec3 E1, vec3 E2, vec2 t0, vec2 t1, vec2 t2, vec3 Z)
{
    const uint32_t start_x = pos_x * RENDER_TILE_SIZE;
    const uint32_t end_x   = start_x + RENDER_TILE_SIZE;

    const uint32_t start_y = pos_y * RENDER_TILE_SIZE;
    const uint32_t end_y   = start_y + RENDER_TILE_SIZE;

    const bool pre_comp_tie_E0 = tie_breaker_ab_test(E0);
    const bool pre_comp_tie_E1 = tie_breaker_ab_test(E1);
    const bool pre_comp_tie_E2 = tie_breaker_ab_test(E2);

    const float step0_x = E0[0] * (float)1;
    const float step1_x = E1[0] * (float)1;
    const float step2_x = E2[0] * (float)1;

    const float step0_y = E0[1] * (float)1;
    const float step1_y = E1[1] * (float)1;
    const float step2_y = E2[1] * (float)1;

    float edge_func_start_0 = (E0[0] * (float)start_x) + (E0[1] * (float)start_y) + E0[2];
    float edge_func_start_1 = (E1[0] * (float)start_x) + (E1[1] * (float)start_y) + E1[2];
    float edge_func_start_2 = (E2[0] * (float)start_x) + (E2[1] * (float)start_y) + E2[2];

    // pre fetch tex data
    const int      tex_w    = raster_state.tex.w;
    const int      tex_h    = raster_state.tex.h;
    const int      tex_bpp  = raster_state.tex.bpp;
    unsigned char *tex_data = raster_state.tex.data;

    vec3 u_values = {t0[0], t1[0], t2[0]};
    vec3 v_values = {t0[1], t1[1], t2[1]};

    float *depth_buffer = rend.depth_buffer;

    for (uint32_t py = start_y; py < end_y; py++,
                  edge_func_start_0 += step0_y,
                  edge_func_start_1 += step1_y,
                  edge_func_start_2 += step2_y

    )
    {
        float edge_func_0 = edge_func_start_0;
        float edge_func_1 = edge_func_start_1;
        float edge_func_2 = edge_func_start_2;

        for (uint32_t px = start_x; px < end_x; px++,
                      edge_func_0 += step0_x,
                      edge_func_1 += step1_x,
                      edge_func_2 += step2_x

        )
        {
            if (edge_tie_breaker(edge_func_0, pre_comp_tie_E0)) continue;
            if (edge_tie_breaker(edge_func_1, pre_comp_tie_E1)) continue;
            if (edge_tie_breaker(edge_func_2, pre_comp_tie_E2)) continue;

            // calculate interpolation coefficients
            const float F0 = edge_func_0;
            const float F1 = edge_func_1;
            const float F2 = edge_func_2;

            const float r = 1.0f / (F0 + F1 + F2);

            vec3 littlef_values = {F0 * r, F1 * r, F2 * r};

            float new_z = Z[0] * littlef_values[0] + Z[1] * littlef_values[1] + Z[2];

            // float  new_z = interpolate_values(littlef_values, Z);
            float *old_z = depth_buffer + ((GRAFIKA_SCREEN_WIDTH * py) + px);

            // depth test
            if (new_z > *old_z) continue;
            *old_z = new_z;

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

            grafika_setpixel(px, py, pixel_colour);
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

void draw_onstart(struct arena *arena)
{
    UNUSED(arena);

    struct obj obj = raster_state.obj;

    ASSERT(obj.mats, "Object must have a material\n");
    ASSERT(obj.mats[0].map_Kd, "Object must have a diffuse texture\n");

    raster_state.tex = tex_load(obj.mats[0].map_Kd);
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
#pragma omp parallel for
    for (size_t i = 0; i < obj.num_pos; i++)
    {
        size_t pos_index = 3 * i;
        vec4   vertex    = {obj_pos[pos_index + 0], obj_pos[pos_index + 1], obj_pos[pos_index + 2], 1.0f};

        m4_mul_v4(raster_state.MVP, vertex, &transformed_vertices[i * 4]);
    }

    // now are verts are in CLIP SPACE
    const size_t          triangle_count = obj.num_f_rows;
    struct triangle_data *triangles      = ma_push_size(&tmp_arena, sizeof(struct triangle_data) * triangle_count);

#if 0
    for (size_t i = 0; i < triangle_count; ++i)
    {
        for (size_t j = 0; j < 3; ++j)
        {
            struct vertindices indices = obj_indices[i * 3 + j];

            const int pos_index = indices.v_idx;
            const int tex_index = indices.vt_idx;
#if 0
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

    // LOG("Triangle count : %u\n", triangle_count);

    // clean out the bins
    memset(tile_array, 0, sizeof(tile_array));
    memset(tiles_with_bins, -1, sizeof(tiles_with_bins));
    tiles_queue_index = 0;



    // now we have all the triangles, in CLIP SPACE
    for (size_t i = 0; i < triangle_count; i++)
    {
        struct triangle_data *t = triangles + i;

        // LOG("v0 : %f, %f, %f, %f\n", v0[0], v0[1], v0[2], v0[3]);
        // LOG("v1 : %f, %f, %f, %f\n", v1[0], v1[1], v1[2], v1[3]);
        // LOG("v2 : %f, %f, %f, %f\n", v2[0], v2[1], v0[2], v0[3]);

        // LOG("t0 : %f, %f\n", t0[0], t0[1]);
        // LOG("t1 : %f, %f\n", t1[0], t1[1]);
        // LOG("t2 : %f, %f\n", t2[0], t2[1]);

        struct aabb aabb;

        // CLIPPER
        if (!_triangle_clipping(t->pos, &aabb))
        {
            // LOG("Triangle %u CLIPPED\n", i);
            //  LOG("v0 : %f, %f, %f\n", v0[0], v0[1], v0[2]);
            //  LOG("v1 : %f, %f, %f\n", v1[0], v1[1], v1[2]);
            //  LOG("v2 : %f, %f, %f\n", v2[0], v2[1], v0[2]);

            continue; // triangle clipped, move to next triangle
        }

        // TRIANGLE SETUP & CULL
        if (!_triangle_setup_and_cull(t->pos, t->E, t->Z))
        {
            // LOG("Triangle %u CULLED\n", i);
            //  LOG("v0 : %f, %f, %f\n", v0[0], v0[1], v0[2]);
            //  LOG("v1 : %f, %f, %f\n", v1[0], v1[1], v1[2]);
            //  LOG("v2 : %f, %f, %f\n", v2[0], v2[1], v0[2]);

            continue; // triangle culled, move to next triangle
        }

        // LOG("Triangle %u will be BINNED\n", i);
        // LOG("v0 : %f, %f, %f\n", v0[0], v0[1], v0[2]);
        // LOG("v1 : %f, %f, %f\n", v1[0], v1[1], v1[2]);
        // LOG("v2 : %f, %f, %f\n", v2[0], v2[1], v0[2]);

        // LOG("AABB : (%f, %f), (%f, %f)\n", aabb.minx, aabb.miny, aabb.maxx, aabb.maxy);

        // BINNER
        // LOG("Tri index : %u\n", i);
        _binner((uint32_t)i, &aabb, t->E);
    }

#endif

    // clean out the bins
    memset(tile_array, 0, sizeof(tile_array));
    memset(tiles_with_bins, -1, sizeof(tiles_with_bins));
    tiles_queue_index = 0;

    TIMED_BLOCK_BEGIN(binner);

#pragma omp parallel for
    for (size_t i = 0; i < triangle_count; ++i)
    {
        struct triangle_data *t = triangles + i;
        struct aabb           aabb;

        struct vertindices indices0 = obj_indices[i * 3 + 0];
        struct vertindices indices1 = obj_indices[i * 3 + 1];
        struct vertindices indices2 = obj_indices[i * 3 + 2];

#ifdef CCW_TRIANGLES
        float *p1    = transformed_vertices + (4 * indices0.v_idx);
        t->pos[0][0] = p1[0];
        t->pos[0][1] = p1[1];
        t->pos[0][2] = p1[2];
        t->pos[0][3] = p1[3];

        float *p2    = transformed_vertices + (4 * indices1.v_idx);
        t->pos[1][0] = p2[0];
        t->pos[1][1] = p2[1];
        t->pos[1][2] = p2[2];
        t->pos[1][3] = p2[3];

        float *p3    = transformed_vertices + (4 * indices2.v_idx);
        t->pos[2][0] = p3[0];
        t->pos[2][1] = p3[1];
        t->pos[2][2] = p3[2];
        t->pos[2][3] = p3[3];
#endif
        // CLIPPER
        if (!_triangle_clipping(t->pos))
        {
            // LOG("Triangle %u CLIPPED\n", i);
            //  LOG("v0 : %f, %f, %f\n", v0[0], v0[1], v0[2]);
            //  LOG("v1 : %f, %f, %f\n", v1[0], v1[1], v1[2]);
            //  LOG("v2 : %f, %f, %f\n", v2[0], v2[1], v0[2]);

            continue; // triangle clipped, move to next triangle
        }

        // TRIANGLE SETUP & CULL
        if (!_triangle_setup_and_cull(t->pos, t->E, t->Z, &aabb))
        {
            // LOG("Triangle %u CULLED\n", i);
            //  LOG("v0 : %f, %f, %f\n", v0[0], v0[1], v0[2]);
            //  LOG("v1 : %f, %f, %f\n", v1[0], v1[1], v1[2]);
            //  LOG("v2 : %f, %f, %f\n", v2[0], v2[1], v0[2]);

            continue; // triangle culled, move to next triangle
        }

#ifdef CCW_TRIANGLES
        v2_copy(t->tex[0], &obj_tex[indices0.vt_idx * 2]);
        v2_copy(t->tex[1], &obj_tex[indices1.vt_idx * 2]);
        v2_copy(t->tex[2], &obj_tex[indices2.vt_idx * 2]);
#endif

        _binner((uint32_t)i, &aabb, t->E);
    }

    TIMED_BLOCK_END(binner);

    // NOTE : here we have collected all "tiles" are in contact with a triangle,
    //          we could go a step further and basically find all pixels that are
    //          touching a triangle, but for now, I am not doing this
#pragma omp parallel for
    for (uint32_t i = 0; i < tiles_queue_index; i++)
    {
        uint32_t tile_lookup_index = tiles_with_bins[i];

        if (tile_lookup_index >= RENDER_TILE_COUNT) continue;

        ASSERT(tile_lookup_index < RENDER_TILE_COUNT, "Lookup tile index is out of range\n");

        struct tile tile_data = tile_array[tile_lookup_index];

        // LOG("Queue index %u, has tile %u, with %u bins\n", i, tile_lookup_index, tile_data.triangle_count);

        // TODO : setup the tiles on startup
        uint32_t tile_x, tile_y;
        _tile_index_to_coords(tile_lookup_index, &tile_x, &tile_y);
        // LOG("Tile coords : %u, %u\n", tile_x, tile_y);

        // #pragma omp parallel for
        for (uint32_t b = 0; b < tile_data.triangle_count; b++)
        {
            struct triangle_ref t = tile_data.triangles[b];

            struct triangle_data *td = triangles + t.triangle_index;

            switch (t.type)
            {
                case RASTER_FULL_TILE:
                {
                    _shade_plain_full_tile(tile_x, tile_y, td->E[0], td->E[1], td->E[2], td->tex[0], td->tex[1], td->tex[2], td->Z);
                    break;
                }
                case RASTER_PARTIAL_TILE:
                {
                    _shade_partial_tile(tile_x, tile_y, td->E[0], td->E[1], td->E[2], td->tex[0], td->tex[1], td->tex[2], td->Z);
                    break;
                }
                case RASTER_QUAD:
                {
                    ASSERT(false, "RASTER_QUAD not implemented\n");
                    break;
                }
                default:
                    ASSERT(false, "Invalid tile rendering\n");
            }
        }
    }
}

void draw_onexit(void)
{
    tex_destroy(&raster_state.tex);
    // raster_state.obj       // cleaned up by arena
}
