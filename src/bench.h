#ifndef __BENCH_H__
#define __BENCH_H__

#include <stdint.h>
#include "utils.h"

enum debug_event_type
{
    debug_event_type_begin_block,
    debug_event_type_end_block,
};

struct debug_record
{
    char    *file_name;
    char    *function_name;
    uint32_t line_number;
};

struct debug_event
{
    uint64_t cycles;
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
};

#define DEBUG_FRAME_COLLECTION_COUNT (16)   // how many frames to collect for
#define DEBUG_MAX_RECORDS            (16)   // each one is a different "TIMED_BLOCK_BEGIN"
#define DEBUG_MAX_EVENTS             (4096) // each one is when a "BEGIN" or "END" block is called

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
    volatile uint64_t table_and_event_index;
    uint32_t          current_table_index;

    struct debug_table table[DEBUG_FRAME_COLLECTION_COUNT]; // frame table?

    uint32_t            record_count;
    struct debug_record records[DEBUG_MAX_RECORDS];

    // this is our final data collected and ready to be used for something
    struct debug_snapshot snapshots[DEBUG_MAX_RECORDS][DEBUG_FRAME_COLLECTION_COUNT];
};

extern struct debug_state *global_debug_state;

#ifdef _MSC_VER
    #include <intrin.h> // For MSVC intrinsic functions
#else
    #include <windows.h>
#endif

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

static inline int debug_timed_block_begin(const int counter, char *file_name, int line_number, char *function_name)
{
    struct debug_record *record = global_debug_state->records + counter;
    record->file_name           = file_name;
    record->function_name       = function_name;
    record->line_number         = (uint32_t)line_number; // NOTE : casting, do we even need this?

    const uint64_t table_and_event_index = Atomic_Add64(&global_debug_state->table_and_event_index, 1);
    const uint32_t table_index           = (const uint32_t)(table_and_event_index >> 32);
    const uint32_t event_index           = (const uint32_t)(table_and_event_index & 0xFFFFFFFF);

    ASSERT(table_index < DEBUG_FRAME_COLLECTION_COUNT, "table_index error");
    ASSERT(event_index < DEBUG_MAX_EVENTS, "event_index error");

    struct debug_event *event = global_debug_state->table[table_index].events + event_index;
    event->record_index       = (uint32_t)counter; // NOTE : casting
    // event->thread_index       = (uint16_t)Get_Thread_ID();                // NOTE : casting
    event->type   = (uint16_t)(debug_event_type_begin_block); // NOTE : casting
    event->cycles = (uint64_t)__rdtsc();                      // NOTE : casting

    return counter; // to avoid returning 0 for the "for" loop
}

static inline int debug_timed_block_end(const int counter, const uint64_t current_cycles)
{
    const uint64_t table_and_event_index = Atomic_Add64(&global_debug_state->table_and_event_index, 1);
    const uint32_t table_index           = (const uint32_t)(table_and_event_index >> 32);
    const uint32_t event_index           = (const uint32_t)(table_and_event_index & 0xFFFFFFFF);
    ASSERT(event_index < DEBUG_MAX_EVENTS, "event index too big!");

    struct debug_event *event = global_debug_state->table[table_index].events + event_index;
    event->record_index       = (uint32_t)counter; // NOTE : casting
    // event->thread_index       = (uint16_t)Get_Thread_ID();              // NOTE : casting
    event->type   = (uint16_t)(debug_event_type_end_block); // NOTE : casting
    event->cycles = current_cycles;

    return 0;
}

#define TOKEN_PASTE(x, y)  x##y
#define TOKEN_PASTE2(x, y) TOKEN_PASTE(x, y)

#define TIMED_BLOCK_NAMED(name)                                                                                \
    for (int TOKEN_PASTE2(counter, __LINE__) = debug_timed_block_begin(__COUNTER__, __FILE__, __LINE__, name); \
         TOKEN_PASTE2(counter, __LINE__);                                                                      \
         TOKEN_PASTE2(counter, __LINE__) = debug_timed_block_end(TOKEN_PASTE2(counter, __LINE__), __rdtsc()))

#define TIMED_BLOCK_BEGIN(name)                  \
    const int _tbb_counter_##name = __COUNTER__; \
    debug_timed_block_begin(_tbb_counter_##name, __FILE__, __LINE__, #name)

#define TIMED_BLOCK_END(name) \
    debug_timed_block_end(_tbb_counter_##name, __rdtsc())

#define TIMED_FUNCTION(func, ...) \
    TIMED_BLOCK_BEGIN(func);      \
    func(__VA_ARGS__);            \
    TIMED_BLOCK_END(func);

#endif // __BENCH_H__