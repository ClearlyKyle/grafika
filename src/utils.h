#ifndef __UTILS_H__
#define __UTILS_H__

#include <assert.h>

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

#ifndef NDEBUG

#define ASSERT(EXP, MSG, ...)                \
    if (!(EXP))                              \
    {                                        \
        LOG("[ASSERT] " MSG, ##__VA_ARGS__); \
        assert(EXP);                         \
    }

#define LOG(FORMAT, ...) fprintf(stderr, "[LOG] " FORMAT, ##__VA_ARGS__)
// #define ASSERT(EXPR, ...)
//     ((EXPR) ? (void)0 : (void)LOG(#EXPR, ##__VA_ARGS__), assert(EXPR))
#else
#define ASSERT(EXP, MSG, ...)
#define STATIC_ASSERT(CONDITION, MSG)
#define LOG(MSG, ...)
#endif

#endif // __UTILS_H__