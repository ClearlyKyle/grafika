#include <stdio.h>
#include <stdlib.h>

#include "raster_mthods/edging.h"
// #include "raster_mthods/stepping.h"
// #include "raster_mthods/matrix.h"
// #include "raster_mthods/tiled.h"
// #include "raster_mthods/avx.h"

int main(int argc, char *argv[])
{
    (void)argc, (void)argv;

    grafika_startup();

    state.obj = obj_load("res/Cube/cube.obj");
    // state.obj = obj_load("plane.obj");
    // state.obj = obj_load("bunny.obj");
    // state.obj = obj_load("res/Dog House/Doghouse.obj");
    // state.obj = obj_load("res/Lorry/lorry.obj");
    // state.obj = obj_load("res/Camera/Camera.obj");
    // state.obj = obj_load("res/Plane/Plane.obj");

    state.tex = tex_load(state.obj.mats[0].map_Kd, true);

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
#ifdef LH_COORDINATE_SYSTEM
    const float model_scale = -1.0f / (bbox.max[1] - bbox.min[1]);
#else
    const float model_scale = 1.0f / (bbox.max[1] - bbox.min[1]);
#endif

    mat4 scale; // scale matrix
    m4_make_scale(model_scale, model_scale, model_scale, scale);

    float rotationAngleX = 0.0f, rotationAngleY = 0.0f;
    float scrollAmount = 2.0f;

    timer_t frame_timer = {0};
    TIMER_START(frame_timer);

    int    frame_counter          = 0;
    double frame_time_accumulated = 0.0;

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
                rotationAngleY += event.motion.xrel * 0.4f;
                rotationAngleX += event.motion.yrel * 0.4f;
            }

            if (event.type == SDL_MOUSEWHEEL)
            {
                scrollAmount -= event.wheel.y * 0.1f;
                if (scrollAmount < 1.5f)
                    scrollAmount = 1.5f;
            }
        }

        mat4 view;
        vec3 eye    = {0.0f, 0.0f, scrollAmount};
        vec3 target = {0.0f, 0.0f, 0.0f};
        vec3 up     = {0.0f, 1.0f, 0.0f};
        m4_lookat(eye, target, up, view);

        // Create rotation matrices
        mat4 rotx, roty;
        m4_make_rot(rotx, DEG2RAD(-rotationAngleX), 'x');
        m4_make_rot(roty, DEG2RAD(-rotationAngleY), 'y');

        // Combine rotations by multiplying matrices (SRT)
        mat4 model, rot;
        m4_mul_m4(rotx, roty, rot);
        m4_mul_m4(rot, trans, model);
        m4_mul_m4(scale, model, model);

        m4_mul_m4(view, model, state.MVP);
        m4_mul_m4(proj, state.MVP, state.MVP);

        grafika_clear();
        {
            draw_object();
            // draw_depthbuffer();
        }
        grafika_present();

        TIMER_UPDATE(frame_timer);

        if (frame_counter == 64)
        {
            const double avg_time = frame_time_accumulated / frame_counter;

            char buff[16] = {0};
            sprintf_s(buff, sizeof(buff), "%0.2fms", avg_time);
            SDL_SetWindowTitle(rend.window, buff);

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
    grafika_destroy();

    return 0;
}