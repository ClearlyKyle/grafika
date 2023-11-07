#ifndef __UTILS_H__
#define __UTILS_H__

#include <assert.h>
#include <stdarg.h>

#define UNUSED(VAR) ((void)(VAR))

#define SAFE_FREE(POINTER_TO_DATA)    \
    {                                 \
        if ((POINTER_TO_DATA))        \
        {                             \
            free(POINTER_TO_DATA);    \
            (POINTER_TO_DATA) = NULL; \
        }                             \
    }
#define ATTRIBATE_FORMAT_PRINTF(VAL1, VAL2) __attribute__((format(gnu_printf, (VAL1), (VAL2))))

static inline void _log_printf(const char *fmt, ...) ATTRIBATE_FORMAT_PRINTF(1, 2);
static inline void _log_printf(const char *fmt, ...)
{
    va_list va_args;
    va_start(va_args, fmt);
    {
        fprintf(stdout, "[LOG] ");
        vfprintf(stdout, fmt, va_args);
    }
    fflush(stdout);
    va_end(va_args);
}

#ifndef NDEBUG
#define LOG(...) _log_printf(__VA_ARGS__)
// #define LOG(format, ...) fprintf(stderr, format, ##__VA_ARGS__)
#define ASSERT(EXP, ...)                \
    if (!(EXP))                         \
    {                                   \
        LOG(__VA_ARGS__); \
        assert(EXP);                    \
    }
#else
#define ASSERT(EXP, MSG, ...)
#define STATIC_ASSERT(CONDITION, MSG)
#define LOG(MSG, ...)
#endif

#endif // __UTILS_H__