#include <stdio.h>
#include <stdlib.h>

// TODO : graphika and shrifty need to be .c files?
// TODO : rename the raster_mthods folder lol
// TODO : hot swap the raster methods
// TODO : fix the threading race condition for the depth testing due to overlapping
//      triangles, each thread needs a render region of the screen, applies to:
//          edging, stepping, phong, normal_map, parallax
// TODO : should raster options be prefixed with something?
//          RASTER_BENCH
//          RASTER_CCW_TRIANGLES
//          RASTER_LH_COOR_SYSTEM

#define BENCH
#define CCW_TRIANGLES
#define LH_COORDINATE_SYSTEM

#include "raster_mthods/edging.c"
// #include "raster_mthods/stepping.c"
// #include "raster_mthods/stepping_simd.c"
// #include "raster_mthods/stepping_simd_jobs.c"
// #include "raster_mthods/matrix.c"
// #include "raster_mthods/matrix_simd.c"
// #include "raster_mthods/phong.c"
// #include "raster_mthods/normal_map.c"
// #include "raster_mthods/parallax.c"
// #include "raster_mthods/binning.c"
// #include "raster_mthods/binning_simd.c"
// #include "raster_mthods/binning_simd_jobs.c"

int main(int argc, char *argv[])
{
    UNUSED(argc), UNUSED(argv);

    struct arena arena = arena_create(MEGABYRES(16));

    grafika_startup(&arena);
    text_startup(rend.surface, 12);

    // terrian created from terrain.c
    // raster_state.obj = obj_load("mountains.obj", &arena);

    // edging, stepping, matrix, avx, phong
    // raster_state.obj = obj_load("res/Cube/cube.obj", &arena);
    // raster_state.obj = obj_load("res/Lorry/lorry.obj", &arena);
    // raster_state.obj = obj_load("res/Camera/Camera.obj", &arena);
    // raster_state.obj = obj_load("res/Plane/Plane.obj", &arena);

    // edging, stepping, matrix, avx, phong, normal mapping
    raster_state.obj = obj_load("res/Wooden Box/wooden crate.obj", &arena);
    // raster_state.obj = obj_load("res/Dog House/Doghouse.obj", &arena);

    // edging, stepping, matrix, avx, phong, normal, parallax
    // raster_state.obj = obj_load("res/Square/square.obj", &arena);

    mat4 proj; // projection matrix
    m4_proj(DEG2RAD(60.0f), (float)GRAFIKA_SCREEN_WIDTH / (float)GRAFIKA_SCREEN_HEIGHT, 0.1f, 1000.0f, proj);

    // calculate model center
    struct bounding_box bbox    = raster_state.obj.bbox;
    float               centerx = -(bbox.min[0] + bbox.max[0]) * 0.5f;
    float               centery = -(bbox.min[1] + bbox.max[1]) * 0.5f;
    float               centerz = -(bbox.min[2] + bbox.max[2]) * 0.5f;
    LOG("BBOX : min(%f, %f, %f), max(%f, %f, %f)\n", bbox.min[0], bbox.min[1], bbox.min[2], bbox.max[0], bbox.max[1], bbox.max[2]);
    LOG("CENTER : %f, %f, %f\n", centerx, centery, centerz);

    // create translation matrix to center the model
    mat4 trans;
    m4_make_trans(centerx, centery, centerz, trans);

    float model_height = (bbox.max[1] - bbox.min[1]);
    float model_scale  = (model_height > 0.0f) ? 1.0f / model_height : 1.0f;
    float cam_z        = max(bbox.min[2], bbox.max[2]) + 1.0f;

    LOG("CAM Z : %f\n", cam_z);
    LOG("SCALE : %f\n", model_scale);

    // create scale matrix
    mat4 scale = {0};
    m4_make_scale(model_scale, model_scale, model_scale, scale);

    struct tymer frame_timer = TIMER_INIT;
    double       avg_time    = 0.0;

    bool  view_update = true;
    float rot_angle_x = 0.0f, rot_angle_y = 0.0f;
    float scroll_amount = cam_z;

    mat4 view;
    vec3 target = {0.0f, 0.0f, 1.0f};
    vec3 up     = {0.0f, 1.0f, 0.0f};
    vec3_set(raster_state.cam_pos, 0.0f, 0.0f, scroll_amount);

    draw_onstart(&arena);

    for (bool running = true; running; /* blank */)
    {
        SDL_Event event;
        while (SDL_PollEvent(&event))
        {
            if (SDL_QUIT == event.type) running = false;

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

        TIMED_BLOCK_NAMED("blank_test") {}

        TIMED_BLOCK_NAMED("new_rot_mats")
        {
            if (view_update)
            {
                raster_state.cam_pos[2] = scroll_amount;
                m4_lookat(raster_state.cam_pos, target, up, view);

                mat4 rotx, roty;
                m4_make_rot(rotx, DEG2RAD(-rot_angle_x), 'x');
                m4_make_rot(roty, DEG2RAD(rot_angle_y), 'y');

                // combine rotations by multiplying matrices (SRT)
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

        TIMED_BLOCK_NAMED("draw_object")
        {
            draw_object(&arena);
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
        text_write(2, 14, "verts : %zu", raster_state.obj.num_verts);
        text_write(2, 26, "cam: (%0.2f, %0.2f, %0.2f)", raster_state.cam_pos[0], raster_state.cam_pos[1], raster_state.cam_pos[2]);

#ifdef BENCH
        debug_frame_end();
        debug_render_info();
#endif

        grafika_present();

        TIMER_UPDATE(frame_timer);
    }

    draw_onexit();

    text_shutdown();
    grafika_shutdown();

    arena_destroy(&arena);

    LOG("Exiting...\n");
    return 0;
}

//
// DEBUG
//
#define BENCH_IMPLEMENTATION
#include "bench.h"