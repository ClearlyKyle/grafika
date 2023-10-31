#ifndef __UTILS_H__
#define __UTILS_H__

#include "assert.h"

#define ALIGN_ME(VAL) __attribute__((aligned((VAL))))

#define SAFE_FREE(POINTER_TO_DATA)    \
    {                                 \
        if ((POINTER_TO_DATA))        \
        {                             \
            free(POINTER_TO_DATA);    \
            (POINTER_TO_DATA) = NULL; \
        }                             \
    }

#ifdef NDEBUG
#define ASSERT(EXPR, ...)             \
    if (!(EXPR))                      \
    {                                 \
        fprintf(stderr, __VA_ARGS__); \
        assert(EXPR);                 \
    }

#define STATIC_ASSERT(CONDITION, MSG)            \
    {                                            \
        typedef char test[(CONDITION) ? 1 : -1]; \
        (void)(test *)#MSG;                      \
    }                                            \
    (void)0

#define LOG(MSG, ...) fprintf(stderr, MSG, __VA_ARGS__)
#else
#define ASSERT(EXPR, ...)
#define STATIC_ASSERT(CONDITION, MSG)
#define LOG(MSG, ...)
#endif

#endif // __UTILS_H__