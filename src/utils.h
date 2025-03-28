#ifndef __UTILS_H__
#define __UTILS_H__

#include <assert.h>
#include <stdarg.h>

#define UNUSED(VAR)      ((void)(VAR))
#define IS_POWER_OF_2(x) ((x) > 0 && ((x) & ((x)-1)) == 0)

#define KILOBYTES(number) ((number)*1024ull)
#define MEGABYRES(number) (KILOBYTES(number) * 1024ull)

#define SAFE_FREE(POINTER_TO_DATA)    \
    {                                 \
        if ((POINTER_TO_DATA))        \
        {                             \
            free(POINTER_TO_DATA);    \
            (POINTER_TO_DATA) = NULL; \
        }                             \
    }

// Log levels
typedef enum
{
    LOG_LEVEL_INFO,
    LOG_LEVEL_WARN,
    LOG_LEVEL_ERROR
} log_level_t;

// Helper macro to get the log level as a string
#define LOG_LEVEL_STR(level) \
    (level == LOG_LEVEL_INFO ? "INFO" : (level == LOG_LEVEL_WARN ? "WARN" : (level == LOG_LEVEL_ERROR ? "ERROR" : "UNKNOWN")))

#define ATTRIBATE_FORMAT_PRINTF(VAL1, VAL2) __attribute__((format(gnu_printf, (VAL1), (VAL2))))

#ifndef NDEBUG
static inline void _log_printf(log_level_t level, const char *file, int line, const char *fmt, ...) ATTRIBATE_FORMAT_PRINTF(4, 5);
static inline void _log_printf(log_level_t level, const char *file, int line, const char *fmt, ...)
{
    va_list args;
    va_start(args, fmt);

    // fprintf(stdout, "[%s] [%s] (%s:%d) ", time_str, LOG_LEVEL_STR(level), file, line);
    fprintf(stdout, "[%s](%s:%d) ", LOG_LEVEL_STR(level), file, line);

    // Print the actual log message
    vfprintf(stdout, fmt, args);

    fflush(stdout); // Flush to ensure the output is immediate
    va_end(args);
}

    #define LOG(level, fmt, ...) _log_printf(level, __FILE__, __LINE__, fmt, ##__VA_ARGS__)
    #define LOG_INFO(fmt, ...)   LOG(LOG_LEVEL_INFO, fmt, ##__VA_ARGS__)
    #define LOG_WARN(fmt, ...)   LOG(LOG_LEVEL_WARN, fmt, ##__VA_ARGS__)
    #define LOG_ERROR(fmt, ...)  LOG(LOG_LEVEL_ERROR, fmt, ##__VA_ARGS__)

    #define ASSERT(expr)                                          \
        do                                                        \
        {                                                         \
            if (!(expr))                                          \
            {                                                     \
                fprintf(stderr, "Assertion failed: %s\n", #expr); \
                abort();                                          \
            }                                                     \
        } while (0)

    #define ASSERTM(expr, fmt, ...)                               \
        do                                                        \
        {                                                         \
            if (!(expr))                                          \
            {                                                     \
                fprintf(stderr, "Assertion failed: %s\n", #expr); \
                LOG(LOG_LEVEL_ERROR, fmt, ##__VA_ARGS__);         \
                abort();                                          \
            }                                                     \
        } while (0)

#else
    #define LOG(level, fmt, ...)
    #define LOG_INFO(fmt, ...)
    #define LOG_WARN(fmt, ...)
    #define LOG_ERROR(fmt, ...)

    #define ASSERT(expr)
    #define ASSERTM(expr, fmt, ...)
#endif

#endif // __UTILS_H__