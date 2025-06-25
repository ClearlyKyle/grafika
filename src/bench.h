#ifndef __BENCH_H__
#define __BENCH_H__

#include <stdint.h>
#include "utils.h"

#ifdef _MSC_VER
#include <intrin.h> // For MSVC intrinsic functions
#else
#include <windows.h>
#endif

//
// TIMER
//

struct tymer
{
    Uint64 start, elapsed, perf_frequency;
};

#define TIMER_INIT                                       \
    (struct tymer)                                       \
    {                                                    \
        .perf_frequency = SDL_GetPerformanceFrequency(), \
        .start = 0, .elapsed = 0,                        \
    }

#define TIMER_START(TIMER)                                      \
    {                                                           \
        (TIMER).perf_frequency = SDL_GetPerformanceFrequency(); \
        (TIMER).start          = SDL_GetPerformanceCounter();   \
        (TIMER).elapsed        = 0;                             \
    }

#define TIMER_UPDATE(TIMER)                                                 \
    {                                                                       \
        const Uint64 _new_time_##TIMER = SDL_GetPerformanceCounter();       \
        (TIMER).elapsed                = _new_time_##TIMER - (TIMER).start; \
        (TIMER).start                  = _new_time_##TIMER;                 \
    }

#define TIMER_ELAPSED_S(TIMER) ((double)(TIMER).elapsed / (double)(TIMER).perf_frequency)

#define TIMER_ELAPSED_MS(TIMER) (TIMER_ELAPSED_S(TIMER) * 1000.0)

//
// PERFORNAMCE BENCHMARKING
//

struct debug_record
{
    char    *file_name;
    char    *function_name;
    uint32_t line_number;
};

struct debug_event
{
    uint64_t cycles;
    uint64_t time;
    uint32_t record_index;
    uint16_t thread_index;
    uint16_t type; // enum debug_event_type
};

struct debug_snapshot
{
    char    *file_name;
    char    *function_name;
    uint32_t line_number;
    uint32_t hit_count;
    uint64_t cycles;
    uint64_t timing;
};

#define DEBUG_FRAME_COLLECTION_COUNT (16)    // how many frames to collect for
#define DEBUG_MAX_RECORDS            (16)    // each one is a different "TIMED_BLOCK_BEGIN"
#define DEBUG_MAX_EVENTS             (65535) // each one is when a "BEGIN" or "END" block is called

struct debug_table
{
    uint32_t           event_count;
    struct debug_event events[DEBUG_MAX_EVENTS];
};

struct debug_state
{
    // "table_and_event_index"
    //  > "table"   in the upper 32 bits, has the index of which column of the events array we are currently
    //              filling. Can think of this as switching to another array at the end of a frame
    //              so we can have some history too
    //  > "event"   index of which event we adding to. Will either store the start or end of a timed section
    //              stores the actual timed period

    // table_and_event_index - holds what current "table" we are loading timing data into
    // previous_table_index - is the last table index in which data was stored into, so will
    //                          be the index of the table with the previous frame events collected
    volatile uint64_t table_and_event_index;
    uint32_t          previous_table_index;

    struct debug_table table[DEBUG_FRAME_COLLECTION_COUNT]; // frame table?

    uint32_t            record_count;
    struct debug_record records[DEBUG_MAX_RECORDS];

    // this is our final data collected and ready to be used for something
    // int arr[2][3] = { {1, 2, 3}, {4, 5, 6} };
    struct debug_snapshot snapshots[DEBUG_MAX_RECORDS][DEBUG_FRAME_COLLECTION_COUNT];
};

extern struct debug_state *global_debug_state;

int debug_timed_block_begin(const int counter, char *file_name, int line_number, char *function_name);
int debug_timed_block_end(const int counter, const uint64_t current_cycles);

void debug_frame_end(void);
void debug_render_info(void);

#define TOKEN_PASTE(x, y)  x##y
#define TOKEN_PASTE2(x, y) TOKEN_PASTE(x, y)

#ifdef BENCH
// __COUNTER__ + 1, the first timed block will be __COUNTER__ == 0, so if we return 0 from debug_timed_block_begin,
// the for loop will not be ran
#define TIMED_BLOCK_NAMED(name)                                                                                    \
    for (int TOKEN_PASTE2(counter, __LINE__) = debug_timed_block_begin(__COUNTER__ + 1, __FILE__, __LINE__, name); \
         TOKEN_PASTE2(counter, __LINE__);                                                                          \
         TOKEN_PASTE2(counter, __LINE__) = debug_timed_block_end(TOKEN_PASTE2(counter, __LINE__), __rdtsc()))

#define TIMED_BLOCK_BEGIN(name) \
    const int _tbb_counter_##name = debug_timed_block_begin(__COUNTER__ + 1, __FILE__, __LINE__, #name)

#define TIMED_BLOCK_END(name) \
    debug_timed_block_end(_tbb_counter_##name, __rdtsc())

// #define TIMED_FUNCTION(func, ...)
//     TIMED_BLOCK_BEGIN(func);
//     func(__VA_ARGS__);
//     TIMED_BLOCK_END(func);
#else
#define TIMED_BLOCK_NAMED(name)
#define TIMED_BLOCK_BEGIN(name)
#define TIMED_BLOCK_END(name)
#define TIMED_FUNCTION(func, ...)
#endif

#endif // __BENCH_H__

#ifdef BENCH_IMPLEMENTATION

enum debug_event_type
{
    debug_event_type_begin_block,
    debug_event_type_end_block,
};

static struct debug_state _global_debug_state = {0};
struct debug_state       *global_debug_state  = &_global_debug_state;

static inline uint32_t Get_Thread_ID(void)
{
    return (uint32_t)GetCurrentThreadId();
}

static inline uint64_t Atomic_Exchange_int64(uint64_t volatile *where_to_put_new_value, uint64_t new_value)
{
#ifdef _MSC_VER
    // Using _InterlockedExchange64 intrinsic for MSVC
    return _InterlockedExchange64((volatile LONG64 *)where_to_put_new_value, (LONG64)new_value);
#else
    // Using atomic_exchange for MinGW (C11)
    // __int64 _InterlockedExchange64(__int64 volatile *Target, __int64 Value)
    return (uint64_t)_InterlockedExchange64((__int64 volatile *)where_to_put_new_value, (__int64)new_value);
#endif
}

// Atomic Add for 32-bit integers
static inline uint32_t Atomic_Add32(uint32_t volatile *original_value, uint32_t value_to_be_added)
{
#ifdef _MSC_VER
    // MSVC and MinGW use the same intrinsic _InterlockedExchangeAdd for 32-bit
    return _InterlockedExchangeAdd((long volatile *)original_value, value_to_be_added);
#else
    // For MinGW and MSVC, they both use _InterlockedExchangeAdd
    // __LONG32 _InterlockedExchangeAdd(__LONG32 volatile *Addend, __LONG32 Value)
    return (uint32_t)InterlockedExchangeAdd((__LONG32 volatile *)original_value, (__LONG32)value_to_be_added);
#endif
}

// Atomic Add for 64-bit integers
static inline uint64_t Atomic_Add64(uint64_t volatile *original_value, uint64_t value_to_be_added)
{
#ifdef _MSC_VER
    // MSVC uses _InterlockedExchangeAdd64 for 64-bit
    return _InterlockedExchangeAdd64((volatile LONG64 *)original_value, value_to_be_added);
#else
    // MinGW uses InterlockedExchangeAdd64 for 64-bit
    // __int64 _InterlockedExchangeAdd64(__int64 volatile *Addend, __int64 Value)
    return (uint64_t)InterlockedExchangeAdd64((__int64 volatile *)original_value, (__int64)value_to_be_added);
#endif
}

int debug_timed_block_begin(const int counter, char *file_name, int line_number, char *function_name)
{
    struct debug_record *record = global_debug_state->records + counter;
    record->file_name           = file_name;
    record->function_name       = function_name;
    record->line_number         = cast_i32_to_u32(line_number);

    const uint64_t table_and_event_index = Atomic_Add64(&global_debug_state->table_and_event_index, 1);
    const uint32_t table_index           = (const uint32_t)(table_and_event_index >> 32);
    const uint32_t event_index           = (const uint32_t)(table_and_event_index & 0xFFFFFFFF);

    ASSERT(table_index < DEBUG_FRAME_COLLECTION_COUNT, "table_index error");
    ASSERT(event_index < DEBUG_MAX_EVENTS, "event_index error %s\n", function_name);

    struct debug_event *event = global_debug_state->table[table_index].events + event_index;
    // event->thread_index       = (uint16_t)Get_Thread_ID(); // NOTE : casting
    event->record_index = cast_i32_to_u32(counter);
    event->type         = (uint16_t)(debug_event_type_begin_block); // NOTE : casting
    event->cycles       = __rdtsc();
    event->time         = SDL_GetPerformanceCounter();

    return counter;
}

int debug_timed_block_end(const int counter, const uint64_t current_cycles)
{
    const uint64_t table_and_event_index = Atomic_Add64(&global_debug_state->table_and_event_index, 1);
    const uint32_t table_index           = (const uint32_t)(table_and_event_index >> 32);
    const uint32_t event_index           = (const uint32_t)(table_and_event_index & 0xFFFFFFFF);
    ASSERT(event_index < DEBUG_MAX_EVENTS, "event index too big!\n");

    struct debug_event *event = global_debug_state->table[table_index].events + event_index;
    // event->thread_index       = (uint16_t)Get_Thread_ID(); // NOTE : casting
    event->record_index = cast_i32_to_u32(counter);
    event->type         = (uint16_t)(debug_event_type_end_block); // NOTE : casting
    event->cycles       = current_cycles;
    event->time         = SDL_GetPerformanceCounter();

    return 0;
}

static void _collate_debug_records(const uint32_t frame_index)
{
    // for each TIMED.. block, we are resetting the hit and cycle count
    // remember the snapshots are what were rendered at the right of our counters,
    // its a history of previous values stored, a saved "snapshot" of previous values
    struct debug_table *table = &global_debug_state->table[frame_index];

    // clear the old collected frame data that was in this slot
    for (uint32_t index = 0; index < DEBUG_MAX_RECORDS; index++)
    {
        struct debug_snapshot *dst = &global_debug_state->snapshots[index][frame_index];
        *dst                       = (struct debug_snapshot){0};

        // dst->cycles    = 0;
        // dst->hit_count = 0;
        // dst->timing    = 0;
    }

    for (uint32_t event_index = 0; event_index < table->event_count; event_index++)
    {
        struct debug_event *event = table->events + event_index;

        struct debug_snapshot *dst_snapshot = &global_debug_state->snapshots[event->record_index][frame_index];
        struct debug_record   *src_record   = global_debug_state->records + event->record_index;

        dst_snapshot->file_name     = src_record->file_name;
        dst_snapshot->function_name = src_record->function_name;
        dst_snapshot->line_number   = src_record->line_number;

        if (event->type == debug_event_type_begin_block)
        {
            // LOG("Begin : %s\n", src_record->function_name);

            dst_snapshot->hit_count += 1;
            dst_snapshot->cycles -= event->cycles;
            dst_snapshot->timing -= event->time;
        }
        else if (event->type == debug_event_type_end_block)
        {
            dst_snapshot->cycles += event->cycles;
            dst_snapshot->timing += event->time;
        }
        else
        {
            assert(false);
        }
    }
}

void debug_frame_end(void)
{
    // swap what the current active table we are using
    // const uint64_t event_array_index = global_debug_state->array_and_event_index >> 32;
    // uint64_t       next_array_index  = event_array_index + 1;

    ++global_debug_state->previous_table_index;
    if (global_debug_state->previous_table_index >= DEBUG_FRAME_COLLECTION_COUNT)
    {
        global_debug_state->previous_table_index = 0;
    }

    const uint64_t table_and_event_index = Atomic_Exchange_int64(&global_debug_state->table_and_event_index,
                                                                 (uint64_t)global_debug_state->previous_table_index << 32);

    const uint32_t table_index = (const uint32_t)(table_and_event_index >> 32);
    const uint32_t event_index = (const uint32_t)(table_and_event_index & 0xFFFFFFFF); // effectivly one past the last event written, a "event_count"

    // LOG("Total events : %u\n", event_index);
    // LOG("Table index  : %u\n", table_index);

    global_debug_state->table[table_index].event_count = event_index; // storing the current count of events

    // table_index is the last table we collected data for
    _collate_debug_records(table_index);
}

struct debug_statistic
{
    double   min, max, avg;
    uint32_t count;
};

#define DEBUG_STAT_INIT {.min = DBL_MAX, .max = -DBL_MAX, .avg = 0.0, .count = 0}

static inline void update_debug_statistic(struct debug_statistic *stat, double val)
{
    stat->count++;

    if (stat->min > val) stat->min = val;
    if (stat->max < val) stat->max = val;

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

void debug_render_info(void)
{
    // uint32_t frame_index = global_debug_state->previous_table_index;

    const uint64_t perf_freq           = SDL_GetPerformanceFrequency();
    const int      text_starting_y_pos = 38;

    // for each timed block we have collected data for
    //      go through the frame history for it and get stats

    // TODO : we should store somewhere the actualy record count we have
    for (uint32_t record_index = 0; record_index < DEBUG_MAX_RECORDS; record_index++)
    {
        struct debug_snapshot *record = &global_debug_state->snapshots[record_index][0];

        if (!record->function_name) continue;

        struct debug_statistic hit_count_stat       = DEBUG_STAT_INIT;
        struct debug_statistic cycle_count_stat     = DEBUG_STAT_INIT;
        struct debug_statistic hits_over_count_stat = DEBUG_STAT_INIT;
        struct debug_statistic timing_stat          = DEBUG_STAT_INIT;

        for (uint32_t frame_index = 0; frame_index < DEBUG_FRAME_COLLECTION_COUNT; frame_index++)
        {
            struct debug_snapshot *snap = record + frame_index;

            update_debug_statistic(&hit_count_stat, (double)snap->hit_count);
            update_debug_statistic(&cycle_count_stat, (double)snap->cycles);

            double ms_time = ((double)snap->timing / (double)perf_freq) * 1000.0;
            update_debug_statistic(&timing_stat, ms_time);

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
        end_debug_statistic(&timing_stat);

        text_write(2, text_starting_y_pos + (12 * (int)record_index),
                   "%-24s %4ucalls %10ucy %09.6fms",
                   record->function_name,
                   (uint32_t)hit_count_stat.avg,   // how many times the function was called
                   (uint32_t)cycle_count_stat.avg, // number of cycles total for all calls of the function
                   (timing_stat.avg)               // average time in ms
        );
    }
}

#endif // BENCH_IMPLEMENTATION