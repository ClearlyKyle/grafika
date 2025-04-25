#include <stdio.h>
#include <stdlib.h>

//#include "graphics/raster_mthods/edging.h"
//#include "graphics/raster_mthods/stepping.h"
//#include "graphics/raster_mthods/matrix.h"
//#include "graphics/raster_mthods/avx.h"
//#include "graphics/raster_mthods/phong.h"
//#include "graphics/raster_mthods/normal_map.h"
#include "graphics/raster_mthods/parallax.h"

int main(int argc, char *argv[])
{
    UNUSED(argc), UNUSED(argv);

    grafika_startup();
    text_startup(GRAFIKA_SCREEN_WIDTH, GRAFIKA_SCREEN_HEIGHT, 12);

    // edging, stepping, matrix, avx, phong
    // state.obj = obj_load("res/Cube/cube.obj");
    // state.obj = obj_load("res/Lorry/lorry.obj");
    // state.obj = obj_load("res/Camera/Camera.obj");
    // state.obj = obj_load("res/Plane/Plane.obj");

    // edging, stepping, matrix, avx, phong, normal mapping
    // state.obj = obj_load("res/Wooden Box/wooden crate.obj");
    // state.obj = obj_load("res/Dog House/Doghouse.obj");

    // edging, stepping, matrix, avx, phong, normal, parallax
    state.obj = obj_load("res/Square/square.obj");

    ASSERT(state.obj.mats, "Object must have atleast a diffuse texture\n");
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

    mat4 scale = {0}; // scale matrix
    m4_make_scale(model_scale, model_scale, model_scale, scale);

    LOG("bounding boxe : x(%f, %f), y(%f, %f), z(%f, %f)\n", bbox.min[0], bbox.max[0], bbox.min[1], bbox.max[1], bbox.min[2], bbox.max[2]);
    LOG("model - scale : %f, center (%f, %f, %f)\n", model_scale, centerx, centery, centerz);

    float rotationAngleX = 0.0f, rotationAngleY = 0.0f;
    float scrollAmount = 4.0f;

    state.cam_pos[0] = 0.0f;
    state.cam_pos[1] = 0.0f;
    state.cam_pos[2] = scrollAmount;

    tymer_t frame_timer = {0};
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
//
// DEBUG
//

static struct debug_state _global_debug_state = {0};
struct debug_state       *global_debug_state  = &_global_debug_state;

static void _collate_debug_records();

static void debug_frame_end(void)
{
    // swap what the current active table we are using
    // const uint64_t event_array_index = global_debug_state->array_and_event_index >> 32;
    // uint64_t       next_array_index  = event_array_index + 1;

    ++global_debug_state->current_table_index;
    if (global_debug_state->current_table_index >= DEBUG_FRAME_COLLECTION_COUNT)
    {
        global_debug_state->current_table_index = 0;
    }

    const uint64_t table_and_event_index = Atomic_Exchange_int64(&global_debug_state->table_and_event_index,
                                                                 (uint64_t)global_debug_state->current_table_index << 32);

    const uint32_t table_index = (const uint32_t)(table_and_event_index >> 32);
    const uint32_t event_index = (const uint32_t)(table_and_event_index & 0xFFFFFFFF); // effectivly one past the last event written, a "event_count"

    // LOG("Total events : %u\n", event_index);

    global_debug_state->table[table_index].event_count = event_index; // storing the current count of events

    _collate_debug_records(table_index);

    // for (uint32_t record_index = 0; record_index < DEBUG_MAX_RECORDS; record_index++)
    //{
    //     struct debug_record *record = global_debug_state->records + record_index;

    //    record->hit_count = 0;
    //    // Atomic_Exchange_uint32(&record->hit_count, (uint64_t)0);
    //}
}

static void _collate_debug_records(const uint32_t frame_index)
{
    // for each TIMED.. block, we are resetting the hit and cycle count
    // remember the snapshots are what were rendered at the right of our counters,
    // its a history of previous values stored, a saved "snapshot" of previous values
    struct debug_table *table = &global_debug_state->table[frame_index];

    // clear the old collected frame data that was in this slot
    for (uint32_t record_index = 0; record_index < DEBUG_MAX_RECORDS; record_index++)
    {
        struct debug_snapshot *dst = global_debug_state->snapshots[record_index] + frame_index;

        dst->hit_count = 0;
        dst->cycles    = 0;
    }

    for (uint32_t event_index = 0; event_index < table->event_count; event_index++)
    {
        struct debug_event *event = table->events + event_index;

        struct debug_snapshot *dst_snapshot = global_debug_state->snapshots[event->record_index] + frame_index;
        struct debug_record   *src_record   = global_debug_state->records + event->record_index;

        dst_snapshot->file_name     = src_record->file_name;
        dst_snapshot->function_name = src_record->function_name;
        dst_snapshot->line_number   = src_record->line_number;

        if (event->type == debug_event_type_begin_block)
        {
            // LOG("Begin : %s\n", src_record->function_name);

            dst_snapshot->hit_count += 1;
            dst_snapshot->cycles -= event->cycles;
        }
        else if (event->type == debug_event_type_end_block)
        {
            dst_snapshot->cycles += event->cycles;
        }
        else
        {
            assert(false);
        }
    }
}

struct debug_statistic
{
    double   min, max, avg;
    uint32_t count;
};

#define DEBUG_STAT_INIT \
    {                   \
        .min = DBL_MAX, .max = DBL_MIN, .avg = 0.0, .count = 0}

static inline void update_debug_statistic(struct debug_statistic *stat, double val)
{
    stat->count++;

    if (stat->min > val)
        stat->min = val;

    if (stat->max < val)
        stat->max = val;

    stat->avg += val;
}

static inline void end_debug_statistic(struct debug_statistic *stat)
{
    if (stat->count)
    {
        stat->avg = stat->avg / (double)stat->count;
    }
    else
    {
        stat->min = 0.0f;
        stat->max = 0.0f;
    }
}

static void debug_render_info(void)
{
    for (uint32_t record_index = 0; record_index < DEBUG_MAX_RECORDS; record_index++)
    {
        struct debug_statistic hit_count_stat       = DEBUG_STAT_INIT;
        struct debug_statistic cycle_count_stat     = DEBUG_STAT_INIT;
        struct debug_statistic hits_over_count_stat = DEBUG_STAT_INIT;

        struct debug_snapshot *snapshot = &global_debug_state->snapshots[record_index][0];

        if (snapshot->function_name)
        {
            for (uint32_t frame_index = 0; frame_index < DEBUG_FRAME_COLLECTION_COUNT; frame_index++)
            {
                struct debug_snapshot *snap = snapshot + frame_index;

                update_debug_statistic(&hit_count_stat, (double)snap->hit_count);
                update_debug_statistic(&cycle_count_stat, (double)snap->cycles);

                double hits_over_count = 0.0;
                if (snap->hit_count > 0)
                {
                    hits_over_count = (double)snap->cycles / (double)snap->hit_count;
                }
                update_debug_statistic(&hits_over_count_stat, hits_over_count);
            }
            end_debug_statistic(&hit_count_stat);
            end_debug_statistic(&cycle_count_stat);
            end_debug_statistic(&hits_over_count_stat);

            text_write(2, 38 + (12 * (int)record_index),
                       "%-24s %4ucalls %10ucy %10uc/call",
                       snapshot->function_name,
                       (uint32_t)hit_count_stat.avg,        // how many times the function was called
                       (uint32_t)cycle_count_stat.avg,      // number of cycles total for all calls of the function
                       (uint32_t)(hits_over_count_stat.avg) // average number of cycles per function call
            );
        }
    }
}