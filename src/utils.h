#ifndef __UTILS_H__
#define __UTILS_H__

#include <assert.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdint.h>

#define UNUSED(VAR)      ((void)(VAR))
#define IS_POWER_OF_2(x) ((x) > 0 && ((x) & ((x) - 1)) == 0)

#define KILOBYTES(number) ((number) * 1024ull)
#define MEGABYRES(number) (KILOBYTES(number) * 1024ull)

//
// LOGGING
//
#define ATTRIBATE_FORMAT_PRINTF(VAL1, VAL2) __attribute__((format(gnu_printf, (VAL1), (VAL2))))

static inline void _log_printf(const char *level, const char *fmt, ...) ATTRIBATE_FORMAT_PRINTF(2, 3);
static inline void _log_printf(const char *level, const char *fmt, ...)
{
    va_list va_args;
    va_start(va_args, fmt);
    {
        fprintf(stdout, "%s ", level);
        vfprintf(stdout, fmt, va_args);
    }
    fflush(stdout);
    va_end(va_args);
}

#ifndef NDEBUG
    #define LOG(...)  _log_printf("[LOG]", __VA_ARGS__)
    #define LOGW(...) _log_printf("[WAR]", __VA_ARGS__)
    #define LOGE(...) _log_printf("[ERR]", __VA_ARGS__)
#else
    #define LOG(...)
    #define LOGW(...)
    #define LOGE(...)
#endif

//
// ASSERT
//
#ifndef NDEBUG
    #define ASSERT(EXP, ...)  \
        if (!(EXP))           \
        {                     \
            LOG(__VA_ARGS__); \
            assert(EXP);      \
        }
#else
    #define ASSERT(EXP, ...)
#endif

//
// ARENA
//
struct arena
{
    uint8_t *base;
    size_t   capacity;
    size_t   used;
};

#define MEGABYTE(n) ((n) * (1024 * 1024))

static inline struct arena arena_create(size_t size)
{
    struct arena arena = {0};

    arena.capacity = size;
    arena.used     = 0;
    arena.base     = malloc(arena.capacity);
    ASSERT(arena.base, "Unable to allocate %ubytes for arena\n", size);

    return arena;
}

static inline void arena_destroy(struct arena *arena)
{
    ASSERT(arena, "Invalid arena to free\n");
    if (arena->base) free(arena->base);
}

static void *arena_alloc_aligned(struct arena *arena, size_t size, size_t alignment)
{
    ASSERT(arena != NULL, "arena is null\n");
    ASSERT(IS_POWER_OF_2(alignment), "alignment must be a power of 2\n");
    ASSERT(arena->used + size + alignment - 1 <= arena->capacity, "oom\n");

    uintptr_t ptr         = (uintptr_t)arena->base + arena->used;
    uintptr_t aligned_ptr = (ptr + alignment - 1) & ~(alignment - 1);
    size_t    padding     = aligned_ptr - ptr;

    arena->used += size + padding;

    LOG("allocated: %zu (req) + %zu (padding), remaining: %zu\n",
        size, padding, arena->capacity - arena->used);

    return (void *)aligned_ptr;
}

static void *ma_push_size(struct arena *arena, const size_t size)
{
    ASSERT(arena != NULL, "arena is null");
    ASSERT(arena->used + size <= arena->capacity, "oom");

    uint8_t *result = arena->base + arena->used;
    assert(result);

    arena->used += size;
    return (void *)result;
}

// main arena will reset after scope exit
static inline struct arena arena_create_tmp(struct arena *a)
{
    return (struct arena){
        .base     = a->base + a->used,
        .capacity = a->capacity - a->used,
        .used     = 0,
    };
}

#define ma_push_struct(arena, type) (type *)ma_push_size((arena), sizeof(type))

#endif // __UTILS_H__