#include <stdio.h>
#include <stdlib.h>

// TODO : move utils to common.h
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
    // LOG("Table index  : %u\n", table_index);

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