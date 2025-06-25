// TODO :
// 5 - colours will be stored in each vertex texture coordinates
//      v1 - x, y, z + u, v     (u here can be Red)
//      v2 - x, y, z + u, v     (u here can be Blue)
//      v3 - x, y, z + u, v     (u here can be Green)
//    v can store surface normals?
#include "raster_mthods/common.h"

#define BENCH

//
// CAMERA
//
struct camera
{
    vec3  target, up, position;
    float yaw, pitch, speed;
};

static void camera_view_matrix(struct camera *cam, mat4 view)
{
    // https://learnopengl.com/Getting-started/Camera

    cam->target[0] = cosf(DEG2RAD(cam->yaw)) * cosf(DEG2RAD(cam->pitch));
    cam->target[1] = sinf(DEG2RAD(cam->pitch));
    cam->target[2] = sinf(DEG2RAD(cam->yaw)) * cosf(DEG2RAD(cam->pitch));
    v3_norm(cam->target);

    vec3 pos_front;
    v3_add(cam->position, cam->target, pos_front);
    m4_lookat(cam->position, pos_front, cam->up, view);
}

//
// TERRAIN
//
struct chunk
{
    int  x, z;
    bool should_render;
};

struct terrain_data
{
    struct chunk chunk_info;

    float *vertices;
    size_t vertex_count;

    int   *indices;
    size_t indices_count;

    size_t triangle_count;
};

struct terrain_data terrain_generate(struct arena *arena, int points, float scale, float offset_x, float offset_z);
static void         terrain_raster(struct terrain_data *terrain, mat4 MVP, vec3 cam_pos);

int main(int argc, char *argv[])
{
    UNUSED(argc), UNUSED(argv);

    struct arena arena = arena_create(MEGABYRES(16));

    grafika_startup(&arena);
    text_startup(rend.surface, 12);

    mat4 proj;
    m4_proj(DEG2RAD(60.0f), (float)GRAFIKA_SCREEN_WIDTH / (float)GRAFIKA_SCREEN_HEIGHT, 0.1f, 1000.0f, proj);

    mat4 trans;
    m4_make_trans(0.0f, 0.0f, 0.0f, trans);

    mat4 scale = {0};
    m4_make_scale(1.0f, 1.0f, 1.0f, scale);

    struct tymer frame_timer = TIMER_INIT;
    double       avg_time    = 0.0;

    struct camera camera = {0};
    camera.yaw           = -90.0f;
    camera.pitch         = 0.0f;
    camera.speed         = 1.5f;
    vec3_set(camera.position, 37.0f, 92.0f, 15.0f);
    vec3_set(camera.target, 0.0f, 0.0f, 1.0f);
    vec3_set(camera.up, 0.0f, -1.0f, 0.0f);

    const float chunk_point_scale    = 32.0f;
    const int   chunk_point_count    = 16;
    const float chunk_total_size     = (float)(chunk_point_count - 1) * chunk_point_scale;
    const float inv_chunk_total_size = 1.0f / chunk_total_size;

    const uint32_t       terrain_list_max   = 64;
    uint32_t             terrain_list_count = 0;
    struct terrain_data *terrain_list       = ma_push_size(&arena, sizeof(struct terrain_data) * terrain_list_max);

    for (bool running = true; running; /* blank */)
    {
        SDL_Event event;
        while (SDL_PollEvent(&event))
        {
            if (SDL_QUIT == event.type) running = false;

            if (event.type == SDL_MOUSEMOTION && (event.motion.state & SDL_BUTTON(SDL_BUTTON_LEFT)))
            {
                camera.yaw += (float)event.motion.xrel;
                camera.pitch += (float)event.motion.yrel;

                camera.pitch = SDL_clamp(camera.pitch, -89.0f, 89.0f);
            }

            if (event.type == SDL_KEYDOWN)
            {
                if (event.key.keysym.scancode == SDL_SCANCODE_W)
                {
                    vec3 tmp;
                    v3_scale(camera.target, camera.speed, tmp);
                    v3_add(camera.position, tmp, camera.position);
                }
                else if (event.key.keysym.scancode == SDL_SCANCODE_S)
                {
                    vec3 tmp;
                    v3_scale(camera.target, camera.speed, tmp);
                    v3_sub(camera.position, tmp, camera.position);
                }
                else if (event.key.keysym.scancode == SDL_SCANCODE_A)
                {
                    //   cameraPos -= glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;
                    vec3 cross;
                    v3_cross(camera.target, camera.up, cross);
                    v3_norm(cross);
                    v3_scale(cross, camera.speed, cross);

                    v3_sub(camera.position, cross, camera.position);
                }
                else if (event.key.keysym.scancode == SDL_SCANCODE_D)
                {
                    //   cameraPos += glm::normalize(glm::cross(cameraFront, cameraUp)) * cameraSpeed;
                    vec3 cross;
                    v3_cross(camera.target, camera.up, cross);
                    v3_norm(cross);
                    v3_scale(cross, camera.speed, cross);

                    v3_add(camera.position, cross, camera.position);
                }
            }
        }

        mat4 view;
        camera_view_matrix(&camera, view);

        // terrain generation
        int cam_chunk_x     = (int)floorf(camera.position[0] * inv_chunk_total_size);
        int cam_chunk_z     = (int)floorf(camera.position[2] * inv_chunk_total_size);
        int render_distance = 3;

        // mark all chunks as not rendered this frame
        for (uint32_t i = 0; i < terrain_list_count; i++)
        {
            terrain_list[i].chunk_info.should_render = false;
        }

        // loop through chunks in range of the camera
        for (int dz = -render_distance; dz <= render_distance; dz++)
        {
            for (int dx = -render_distance; dx <= render_distance; dx++)
            {
                int cx = cam_chunk_x + dx;
                int cz = cam_chunk_z + dz;

                bool found = false;

                // search for existing chunk
                for (uint32_t i = 0; i < terrain_list_count; i++)
                {
                    if (terrain_list[i].chunk_info.x == cx && terrain_list[i].chunk_info.z == cz)
                    {
                        // LOG("Reusing OLD Chunk (%d, %d)\n", cx, cz);

                        terrain_list[i].chunk_info.should_render = true;
                        found                                    = true;
                        break;
                    }
                }

                // if not found, generate it
                if (!found && terrain_list_count < terrain_list_max)
                {
                    // LOG("Generating NEW Chunk (%d, %d)\n", cx, cz);
                    // LOG("total chunks made : %u\n", terrain_list_count + 1);

                    float offset_x = (float)cx * chunk_total_size;
                    float offset_z = (float)cz * chunk_total_size;

                    terrain_list[terrain_list_count] = terrain_generate(&arena, chunk_point_count, chunk_point_scale, offset_x, offset_z);

                    terrain_list[terrain_list_count].chunk_info.x             = cx;
                    terrain_list[terrain_list_count].chunk_info.z             = cz;
                    terrain_list[terrain_list_count].chunk_info.should_render = true;

                    terrain_list_count++;
                }
            }
        }

        mat4 MVP   = MAT4_IDENTITY;
        mat4 model = MAT4_IDENTITY;

        m4_mul_m4(view, model, MVP);
        m4_mul_m4(proj, MVP, MVP);

        grafika_clear();

        size_t vert_counter = 0;

        for (uint32_t chunk_index = 0; chunk_index < terrain_list_count; chunk_index++)
        {
            struct terrain_data terrain = terrain_list[chunk_index];

            vert_counter += terrain.vertex_count;

            terrain_raster(&terrain, MVP, camera.position);
        }

#if 0 // Cumulative Moving Average (CMA)
        double current_time = TIMER_ELAPSED_MS(frame_timer);
        avg_time            = (avg_time * frame_counter + current_time) / (frame_counter + 1);
        frame_counter       = (frame_counter + 1) % 64; // Reset periodically to avoid overflow
#else //  Exponential Moving Average (EMA)
        double current_time = TIMER_ELAPSED_MS(frame_timer);
        avg_time            = 0.1 * current_time + (1.0f - 0.1) * avg_time;
#endif

        text_write(2, 2, "frame : %0.2fms", avg_time);
        text_write(2, 14, "verts : %zu", vert_counter);
        text_write(2, 26, "cam: (%0.2f, %0.2f, %0.2f)", camera.position[0], camera.position[1], camera.position[2]);

#ifdef BENCH
        debug_frame_end();
        debug_render_info();
#endif

        grafika_present();

        TIMER_UPDATE(frame_timer);
    }

    text_shutdown();
    grafika_shutdown();

    arena_destroy(&arena);

    LOG("Exiting...\n");
    return 0;
}

static float noise2d(int x, int y, int seed)
{
    int n = x + y * 57 + seed * 131;
    n     = (n << 13) ^ n;
    return (1.0f - (float)((n * (n * n * 15731 + 789221) + 1376312589) & 0x7fffffff) / 1073741824.0f);
}

// smooth interpolation between noise points
static float interpolate(float a, float b, float x)
{
    float ft = x * 3.1415927f;
    float f  = (1.0f - cosf(ft)) * 0.5f;
    return a * (1.0f - f) + b * f;
}

static float perlin_noise(float x, float y, int seed)
{
    int   x_int  = (int)x;
    int   y_int  = (int)y;
    float x_frac = x - (float)x_int;
    float y_frac = y - (float)y_int;

    float v1 = noise2d(x_int, y_int, seed);
    float v2 = noise2d(x_int + 1, y_int, seed);
    float v3 = noise2d(x_int, y_int + 1, seed);
    float v4 = noise2d(x_int + 1, y_int + 1, seed);

    float i1 = interpolate(v1, v2, x_frac);
    float i2 = interpolate(v3, v4, x_frac);

    return interpolate(i1, i2, y_frac);
}

static float calculate_elevation(float pos_x, float pos_y)
{
    const int   layer_count       = 8;
    const float lacunarity        = 2.0f;
    const float presistence       = 1.5f;
    const int   ridge_layer_start = 2;

    float frequency = 0.002f;
    float amplitude = 1.0f;
    float elevation = 0.0f;

    for (int i = 0; i < layer_count; i++)
    {
        float noise = perlin_noise(pos_x * frequency, pos_y * frequency, 12345);
        // float noise = noise_gen_evaluate(pos_x * frequency, pos_y * frequency);
        if (i >= ridge_layer_start) noise = 0.5f - fabsf(noise);

        elevation += noise * amplitude;
        amplitude -= presistence;
        frequency *= lacunarity;
    }

    return elevation * 15.0f;
}

static void calculate_jiggle(float in[2], float out[2])
{
    const float jiggle = 0.7f;

    float x = perlin_noise(in[0], in[1], 12345) * jiggle;
    float y = perlin_noise(in[0] - 10000.0f, in[1] - 10000.0f, 12345) * jiggle;

    out[0] = x;
    out[1] = y;
}

struct terrain_data terrain_generate(struct arena *arena, int points, float scale, float offset_x, float offset_z)
{
    size_t total_verts    = (size_t)(points * points);
    size_t indicies_count = (size_t)((points - 1) * (points - 1) * 6);

    struct terrain_data terrain = {0};

    float *vertices = ma_push_size(arena, sizeof(float) * 3 * total_verts);
    int   *indices  = ma_push_size(arena, sizeof(int) * indicies_count);

    int vertex_idx = 0;
    for (int x = 0; x < points; x++)
    {
        for (int y = 0; y < points; y++)
        {
            float x_pos = (float)x * scale + offset_x;
            float y_pos = (float)y * scale + offset_z;

            // calculate_jiggle(point, point);

            float elevation = fmaxf(0.0f, calculate_elevation(x_pos, y_pos));
            // float elevation = 1.0f;

            vertices[vertex_idx++] = x_pos;     // X
            vertices[vertex_idx++] = elevation; // Y
            vertices[vertex_idx++] = y_pos;     // Z
        }
    }

    // generate indices
    size_t trianlge_count = 0;
    int    index_idx      = 0;
    for (int z = 0; z < points - 1; z++)
    {
        for (int x = 0; x < points - 1; x++)
        {
            // Get indices of the 4 corners of the quad
            const int top_left     = z * points + x;
            const int top_right    = top_left + 1;
            const int bottom_left  = (z + 1) * points + x;
            const int bottom_right = bottom_left + 1;

#if 0 // CW
      // Triangle 1
            trianlge_count++;
            indices[index_idx].vt_idx  = top_left;
            indices[index_idx++].v_idx = top_left;
            
            indices[index_idx].vt_idx  = bottom_left;
            indices[index_idx++].v_idx = bottom_left;
            
            indices[index_idx].vt_idx  = top_right;
            indices[index_idx++].v_idx = top_right;
            
            trianlge_count++;
            // Triangle 2
            indices[index_idx].vt_idx  = top_right;
            indices[index_idx++].v_idx = top_right;
            
            indices[index_idx].vt_idx  = bottom_left;
            indices[index_idx++].v_idx = bottom_left;
            
            indices[index_idx].vt_idx  = bottom_right;
            indices[index_idx++].v_idx = bottom_right;
#else // CCW
      // Triangle 1
            trianlge_count++;

            indices[index_idx++] = top_left;
            indices[index_idx++] = top_right;
            indices[index_idx++] = bottom_left;

            trianlge_count++;

            // Triangle 2
            indices[index_idx++] = top_right;
            indices[index_idx++] = bottom_right;
            indices[index_idx++] = bottom_left;
#endif
        }
    }

    terrain.vertices     = vertices;
    terrain.vertex_count = total_verts * 3;

    terrain.indices       = indices;
    terrain.indices_count = indicies_count;

    terrain.triangle_count = trianlge_count;

    return terrain;
}

//
// RASTER TERRAIN
//

static void draw_AABB_outline(int AABB[4], uint32_t colour)
{
    int minX = AABB[0];
    int minY = AABB[1];
    int maxX = AABB[2];
    int maxY = AABB[3];

    for (int x = minX; x <= maxX; ++x)
    {
        grafika_setpixel((uint32_t)x, (uint32_t)minY, colour); // Top
        grafika_setpixel((uint32_t)x, (uint32_t)maxY, colour); // Bottom
    }

    for (int y = minY; y <= maxY; ++y)
    {
        grafika_setpixel((uint32_t)minX, (uint32_t)y, colour); // Left
        grafika_setpixel((uint32_t)maxX, (uint32_t)y, colour); // Right
    }
}

static void draw_triangle(mat4 MVP, vec3 verts[3], vec3 cam_pos)
{
    float max_y = max(verts[0][1], max(verts[1][1], verts[2][1]));

    int colour[3] = {255, 255, 255}; // snow tops
    if (max_y < 1.5f)
    {
        //  water
        colour[0] = 51;  // R
        colour[1] = 125; // G
        colour[2] = 245; // B
    }
    else if (max_y < 10.0f)
    {
        //  grass
        colour[0] = 51;  // R
        colour[1] = 125; // G
        colour[2] = 25;  // B
    }
    else if (max_y < 55.0f)
    {
        // rocky
        colour[0] = 139;
        colour[1] = 126;
        colour[2] = 102;
    }

    // convert to clip space
    vec4 clip_space[3];
    for (int i = 0; i < 3; ++i)
    {
        vec4 pos = {verts[i][0], verts[i][1], verts[i][2], 1.0f};
        m4_mul_v4(MVP, pos, clip_space[i]);
    }

    // clipping (after all vertices transformed)
    bool outside = false;
    for (int plane = 0; plane < 6; ++plane)
    {
        int count = 0;
        for (int j = 0; j < 3; ++j)
        {
            switch (plane)
            {
                case 0:
                    if (clip_space[j][0] < -clip_space[j][3]) count++;
                    break; // left
                case 1:
                    if (clip_space[j][0] > clip_space[j][3]) count++;
                    break; // right
                case 2:
                    if (clip_space[j][1] < -clip_space[j][3]) count++;
                    break; // bottom
                case 3:
                    if (clip_space[j][1] > clip_space[j][3]) count++;
                    break; // top
                case 4:
                    if (clip_space[j][2] < -clip_space[j][3]) count++;
                    break; // near
                case 5:
                    if (clip_space[j][2] > clip_space[j][3]) count++;
                    break; // far
                default: break;
            }
        }
        if (count == 3)
        {
            outside = true;
            break;
        }
    }
    if (outside) return;

    // perspective division (clip to ndc)
    vec3 ndc[3], w_vals;
    for (size_t i = 0; i < 3; i++)
    {
        w_vals[i] = 1.0f / clip_space[i][3]; // 1.0f / w
        ndc[i][0] = clip_space[i][0] * w_vals[i];
        ndc[i][1] = clip_space[i][1] * w_vals[i];
        ndc[i][2] = clip_space[i][2] * w_vals[i];
    }

    vec3 screen_space[3];
    for (int i = 0; i < 3; ++i)
    {
        screen_space[i][0] = (ndc[i][0] + 1.0f) * 0.5f * GRAFIKA_SCREEN_WIDTH;
        screen_space[i][1] = (ndc[i][1] + 1.0f) * 0.5f * GRAFIKA_SCREEN_HEIGHT;
        screen_space[i][2] = (ndc[i][2] + 1.0f) * 0.5f;
    }

    float dx1     = screen_space[1][0] - screen_space[0][0];
    float dy1     = screen_space[1][1] - screen_space[0][1];
    float dx2     = screen_space[2][0] - screen_space[0][0];
    float dy2     = screen_space[2][1] - screen_space[0][1];
    float cross_z = dx1 * dy2 - dy1 * dx2;

    if (cross_z >= 0.0f) return; // backface

    // compute triangle normal
    vec3 ab, ac, normal;
    v3_sub(verts[2], verts[0], ab);
    v3_sub(verts[1], verts[0], ac);
    v3_cross(ab, ac, normal);
    v3_norm(normal);

    // vector from triangle to camera
    vec3 to_camera;
    v3_sub(cam_pos, verts[0], to_camera); // use triangle center if you prefer
    v3_norm(to_camera);

    // Lambertian brightness
    float brightness = max(0.0f, v3_dot(normal, to_camera));

    float ambient          = 0.2f;
    float final_brightness = SDL_clamp(ambient + (1.0f - ambient) * brightness, 0.0f, 1.0f);

    int shaded_colour[3];
    for (int i = 0; i < 3; ++i)
    {
        shaded_colour[i] = (int)((float)colour[i] * final_brightness);

        shaded_colour[i] = SDL_clamp(shaded_colour[i], 0, 255);
    }

    int AABB[4] = {0};
    AABB_make(screen_space, AABB);

    // draw_AABB_outline(AABB, 0xFFFFFFFF);

    const float dY0 = screen_space[2][1] - screen_space[1][1], dX0 = screen_space[1][0] - screen_space[2][0];
    const float dY1 = screen_space[0][1] - screen_space[2][1], dX1 = screen_space[2][0] - screen_space[0][0];
    const float dY2 = screen_space[1][1] - screen_space[0][1], dX2 = screen_space[0][0] - screen_space[1][0];

    const float C0 = (screen_space[2][0] * screen_space[1][1]) - (screen_space[2][1] * screen_space[1][0]);
    const float C1 = (screen_space[0][0] * screen_space[2][1]) - (screen_space[0][1] * screen_space[2][0]);
    const float C2 = (screen_space[1][0] * screen_space[0][1]) - (screen_space[1][1] * screen_space[0][0]);

    const float inv_area = 1.0f / (dX1 * dY2 - dY1 * dX2);

    float alpha = (dY0 * ((float)AABB[0])) + (dX0 * ((float)AABB[1])) + C0;
    float betaa = (dY1 * ((float)AABB[0])) + (dX1 * ((float)AABB[1])) + C1;
    float gamma = (dY2 * ((float)AABB[0])) + (dX2 * ((float)AABB[1])) + C2;

    screen_space[1][2] = (screen_space[1][2] - screen_space[0][2]) * inv_area;
    screen_space[2][2] = (screen_space[2][2] - screen_space[0][2]) * inv_area;

    const float zstep = dY1 * screen_space[1][2] + dY2 * screen_space[2][2];

    float *depth_buffer = rend.depth_buffer;

    // brightness = max(0.0f, dot(normal, view_dir));

    for (int y = AABB[1]; y <= AABB[3]; ++y)
    {
        // barycentric coordinates at start of row
        float w0 = alpha;
        float w1 = betaa;
        float w2 = gamma;

        float depth = screen_space[0][2] + (screen_space[1][2] * betaa) + (screen_space[2][2] * gamma);

        for (int x = AABB[0]; x <= AABB[2]; ++x,            //
                                            w0 += dY0,      //
                                            w1 += dY1,      //
                                            w2 += dY2,      // one step to the right
                                            depth += zstep) // step the depth
        {
            if (w0 < 0.0f || w1 < 0.0f || w2 < 0.0f)
                continue;

            const int index = (y * GRAFIKA_SCREEN_WIDTH) + x;

            float      *oldZ = depth_buffer + index;
            const float invZ = depth;

            if (invZ > *oldZ) continue;
            *oldZ = invZ;

            uint32_t pixel_colour = ((uint32_t)shaded_colour[0] << 16) |
                                    ((uint32_t)shaded_colour[1] << 8) |
                                    (uint32_t)shaded_colour[2];

            grafika_setpixel((uint32_t)x, (uint32_t)y, pixel_colour);
        }

        // step one row
        alpha += dX0;
        betaa += dX1;
        gamma += dX2;
    }
}

static void terrain_raster(struct terrain_data *terrain, mat4 MVP, vec3 cam_pos)
{
    // #pragma omp parallel
    {
        size_t triangle_count = terrain->triangle_count;

        // #pragma omp for
        for (size_t i = 0; i < triangle_count; ++i)
        {
            vec3 verts[3];

            for (size_t j = 0; j < 3; ++j)
            {
                int pos_index = terrain->indices[i * 3 + j];
#if 0
                const size_t store_index = j; // CCW triangles
#else
                const size_t store_index = 2 - j; // CW triangles (which blender uses)
#endif
                verts[store_index][0] = terrain->vertices[3 * pos_index + 0];
                verts[store_index][1] = terrain->vertices[3 * pos_index + 1];
                verts[store_index][2] = terrain->vertices[3 * pos_index + 2];
            }
            draw_triangle(MVP, verts, cam_pos);
        }
    }
}

//
// DEBUG
//
#define BENCH_IMPLEMENTATION
#include "bench.h"