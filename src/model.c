#include <stdio.h>
#include <stdlib.h>

// #include "graphics/raster_mthods/edging.h"
// #include "graphics/raster_mthods/stepping.h"
// #include "graphics/raster_mthods/matrix.h"
// #include "graphics/raster_mthods/avx.h"
// #include "graphics/raster_mthods/phong.h"
// #include "graphics/raster_mthods/normal_map.h"
#include "graphics/raster_mthods/parallax.h"

int main(int argc, char *argv[])
{
    UNUSED(argc), UNUSED(argv);

    grafika_startup();
    text_startup(GRAFIKA_SCREEN_WIDTH, GRAFIKA_SCREEN_HEIGHT, 12);

    state.obj = obj_load("res/Square/square.obj"); // parallax mapping
    // state.obj = obj_load("res/Cube/cube.obj");
    // state.obj = obj_load("bunny.obj");
    // state.obj = obj_load("res/Wooden Box/wooden crate.obj"); // normal mapping
    // state.obj = obj_load("res/Dog House/Doghouse.obj");      // normal mapping
    // state.obj = obj_load("res/Lorry/lorry.obj");
    // state.obj = obj_load("res/Camera/Camera.obj");
    // state.obj = obj_load("res/Plane/Plane.obj");

    state.tex = tex_load(state.obj.mats[0].map_Kd);

    mat4 proj; // projection matrix
    m4_proj(DEG2RAD(60.0f), (float)GRAFIKA_SCREEN_WIDTH / (float)GRAFIKA_SCREEN_HEIGHT, 0.1f, 100.0f, proj);

    // Calculate model height
    BoundingBox_t bbox    = state.obj.bbox;
    float         centerx = -(bbox.min[0] + bbox.max[0]) * 0.5f;
    float         centery = -(bbox.min[1] + bbox.max[1]) * 0.5f;
    float         centerz = -(bbox.min[2] + bbox.max[2]) * 0.5f;

    mat4 trans; // translation matrix
    m4_make_trans(centerx, centery, centerz, trans);

    // Scale model height to 1
    float model_scale = 1.0f;
#ifdef LH_COORDINATE_SYSTEM
    const float model_scale = -1.0f / (bbox.max[1] - bbox.min[1]);
#else
    if ((centerx + centery + centerz) != 0.0f)
        model_scale = 1.0f / (bbox.max[1] - bbox.min[1]);
#endif

    mat4 scale; // scale matrix
    m4_make_scale(model_scale, model_scale, model_scale, scale);

    LOG("bounding boxe : x(%f, %f), y(%f, %f), z(%f, %f)\n", bbox.min[0], bbox.max[0], bbox.min[1], bbox.max[1], bbox.min[2], bbox.max[2]);
    LOG("model - scale : %f, center (%f, %f, %f)\n", model_scale, centerx, centery, centerz);

    float rotationAngleX = 0.0f, rotationAngleY = 0.0f;
    float scrollAmount = 4.0f;

    state.cam_pos[0] = 0.0f;
    state.cam_pos[1] = 0.0f;
    state.cam_pos[2] = scrollAmount;

    timer_t frame_timer = {0};
    TIMER_START(frame_timer);

    int    frame_counter          = 0;
    double frame_time_accumulated = 0.0, avg_time = 0.0;

    while (!rend.quit)
    {
        SDL_Event event;
        while (SDL_PollEvent(&event))
        {
            if (SDL_QUIT == event.type)
            {
                rend.quit = true;
            }

            if (event.type == SDL_MOUSEMOTION && (event.motion.state & SDL_BUTTON(SDL_BUTTON_LEFT)))
            {
                rotationAngleY += (float)event.motion.xrel * 0.4f; // Cap these angles?
                rotationAngleX += (float)event.motion.yrel * 0.4f;
            }

            if (event.type == SDL_MOUSEWHEEL)
            {
                scrollAmount -= (float)event.wheel.y * 0.1f;
                if (scrollAmount < 1.5f)
                    scrollAmount = 1.5f;
            }
        }

        mat4 view;
        state.cam_pos[2] = scrollAmount;
        vec3 target      = {0.0f, 0.0f, 0.0f};
        vec3 up          = {0.0f, 1.0f, 0.0f};
        m4_lookat(state.cam_pos, target, up, view);

        // Create rotation matrices
        mat4 rotx, roty;
        m4_make_rot(rotx, DEG2RAD(-rotationAngleX), 'x');
        m4_make_rot(roty, DEG2RAD(-rotationAngleY), 'y');

        // Combine rotations by multiplying matrices (SRT)
        mat4 rot;
        m4_mul_m4(rotx, roty, rot);
        m4_mul_m4(rot, trans, state.model);
        m4_mul_m4(scale, state.model, state.model);

        m4_mul_m4(view, state.model, state.MVP);
        m4_mul_m4(proj, state.MVP, state.MVP);

        grafika_clear();
        {
            draw_object();
        }
        grafika_present();

        TIMER_UPDATE(frame_timer);

        text_write(2, 2, "%0.2fms", avg_time);
        text_write(2, 12, "verts : %zu", state.obj.num_verts);
        text_write(2, 22, "cam: (%0.2f, %0.2f, %0.2f)", state.cam_pos[0], state.cam_pos[1], state.cam_pos[2]);

        if (frame_counter == 32)
        {
            avg_time               = frame_time_accumulated / frame_counter;
            frame_counter          = 0;
            frame_time_accumulated = 0.0;
        }
        else
        {
            frame_time_accumulated += TIMER_ELAPSED_MS(frame_timer);
            frame_counter++;
        }
    }

    obj_destroy(&state.obj);
    draw_onexit();
    grafika_shutdown();
    text_shutdown();

    return 0;
}