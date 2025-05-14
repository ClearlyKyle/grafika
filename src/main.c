#include <stdio.h>
#include <stdlib.h>

// TODO : move the debug stuff to another file
// TODO : run through all files and make sure all modes work
// TODO : add the job system
// TODO : make sure its all camel case
// TODO : the todos in the obj loader
// TODO : switch the obj loader code to a single header style
// TODO : move the timer code to the bench marking

// TODO : graphika and shrifty need to be .c files?

#define BENCH

// #include "raster_mthods/edging.c"
// #include "raster_mthods/stepping.c"
// #include "raster_mthods/stepping_simd.c"
// #include "raster_mthods/matrix.c"
// #include "raster_mthods/matrix_simd.c"
// #include "raster_mthods/phong.c"
// #include "raster_mthods/normal_map.c"
#include "raster_mthods/parallax.c"

static void debug_frame_end(void);
static void debug_render_info(void);

int main(int argc, char *argv[])
{
    UNUSED(argc), UNUSED(argv);

    struct arena arena = arena_create(MEGABYRES(8));

    grafika_startup(&arena);
    text_startup(rend.surface, 12);

    // edging, stepping, matrix, avx, phong
    // raster_state.obj = obj_load("res/Cube/cube.obj");
    // raster_state.obj = obj_load("res/Lorry/lorry.obj");
    // raster_state.obj = obj_load("res/Camera/Camera.obj");
    // raster_state.obj = obj_load("res/Plane/Plane.obj");

    // edging, stepping, matrix, avx, phong, normal mapping
    // raster_state.obj = obj_load("res/Wooden Box/wooden crate.obj");
    // raster_state.obj = obj_load("res/Dog House/Doghouse.obj");

    // edging, stepping, matrix, avx, phong, normal, parallax
    raster_state.obj = obj_load("res/Square/square.obj");

    // TODO : move this to the drawing functions?
    ASSERT(raster_state.obj.mats, "Object must have atleast a diffuse texture\n");
    raster_state.tex = tex_load(raster_state.obj.mats[0].map_Kd);

    mat4 proj; // projection matrix
    m4_proj(DEG2RAD(60.0f), (float)GRAFIKA_SCREEN_WIDTH / (float)GRAFIKA_SCREEN_HEIGHT, 0.1f, 100.0f, proj);

    // calculate model height
    BoundingBox_t bbox    = raster_state.obj.bbox;
    float         centerx = -(bbox.min[0] + bbox.max[0]) * 0.5f;
    float         centery = -(bbox.min[1] + bbox.max[1]) * 0.5f;
    float         centerz = -(bbox.min[2] + bbox.max[2]) * 0.5f;

    mat4 trans; // translation matrix
    m4_make_trans(centerx, centery, centerz, trans);

    // scale model height to 1
    float model_scale = 1.0f;
#ifdef LH_COORDINATE_SYSTEM
    const float model_scale = -1.0f / (bbox.max[1] - bbox.min[1]);
#else
    if ((centerx + centery + centerz) != 0.0f)
        model_scale = 1.0f / (bbox.max[1] - bbox.min[1]);
#endif

    mat4 scale = {0}; // scale matrix
    m4_make_scale(model_scale, model_scale, model_scale, scale);

    bool  view_update = true;
    float rot_angle_x = 0.0f, rot_angle_y = 0.0f;
    float scroll_amount = 4.0f;

    raster_state.cam_pos[0] = 0.0f;
    raster_state.cam_pos[1] = 0.0f;
    raster_state.cam_pos[2] = scroll_amount;

    struct tymer frame_timer = TIMER_INIT;

    int    frame_counter          = 0;
    double frame_time_accumulated = 0.0, avg_time = 0.0; // TODO : remove avg_time

    mat4 view = {0};

    for (bool running = true; running; /* blank */)
    {
        TIMED_BLOCK_NAMED("poll_events")
        {
            SDL_Event event;
            while (SDL_PollEvent(&event))
            {
                if (SDL_QUIT == event.type)
                {
                    running = false;
                }

                if (event.type == SDL_MOUSEMOTION && (event.motion.state & SDL_BUTTON(SDL_BUTTON_LEFT)))
                {
                    view_update = true;
                    rot_angle_y += (float)event.motion.xrel * 0.4f; // Cap these angles?
                    rot_angle_x += (float)event.motion.yrel * 0.4f;
                }

                if (event.type == SDL_MOUSEWHEEL)
                {
                    view_update = true;

                    scroll_amount -= (float)event.wheel.y * 0.1f;

                    if (scroll_amount < 1.5f) scroll_amount = 1.5f;
                }
            }
        }

        if (view_update)
        {
            TIMED_BLOCK_NAMED("new_rot_mats")
            {
                raster_state.cam_pos[2] = scroll_amount;
                vec3 target             = {0.0f, 0.0f, 0.0f};
                vec3 up                 = {0.0f, -1.0f, 0.0f};
                m4_lookat(raster_state.cam_pos, target, up, view);

                // Create rotation matrices
                mat4 rotx, roty;
                m4_make_rot(rotx, DEG2RAD(-rot_angle_x), 'x');
                m4_make_rot(roty, DEG2RAD(rot_angle_y), 'y');

                // Combine rotations by multiplying matrices (SRT)
                mat4 rot;
                m4_mul_m4(rotx, roty, rot);
                m4_mul_m4(rot, trans, raster_state.model);
                m4_mul_m4(scale, raster_state.model, raster_state.model);

                m4_mul_m4(view, raster_state.model, raster_state.MVP);
                m4_mul_m4(proj, raster_state.MVP, raster_state.MVP);

                view_update = false;
            }
        }

        TIMED_BLOCK_NAMED("grafika_clear")
        {
            grafika_clear();
        }

        draw_object(&arena);

        text_write(2, 2, "%0.2fms", avg_time);
        text_write(2, 14, "verts : %zu", raster_state.obj.num_verts);
        text_write(2, 26, "cam: (%0.2f, %0.2f, %0.2f)", raster_state.cam_pos[0], raster_state.cam_pos[1], raster_state.cam_pos[2]);

#ifdef BENCH
        debug_frame_end();
        debug_render_info();
#endif

        frame_time_accumulated += TIMER_ELAPSED_MS(frame_timer);

        if (frame_counter++ == 64)
        {
            avg_time               = frame_time_accumulated / 64.0;
            frame_counter          = 0;
            frame_time_accumulated = 0.0;
        }

        grafika_present();

        TIMER_UPDATE(frame_timer);
    }

    obj_destroy(&raster_state.obj);

    text_shutdown();
    grafika_shutdown();

    arena_destroy(&arena);

    LOG("Exiting...\n");
    return 0;

    return 0;
}

//
// DEBUG
//
#define BENCH_IMPLEMENTATION
#include "bench.h"