#ifndef __UTILS_H__
#define __UTILS_H__

#ifndef NDEBUG
    #define ASSERT(EXPR, ...)                                                   \
        if (!(EXPR))                                                            \
        {                                                                       \
            fprintf(stderr, "[ASSERT] %s %s:%d - ", #EXPR, __FILE__, __LINE__); \
            fprintf(stderr, __VA_ARGS__);                                       \
            abort();                                                            \
        }
    #define SAFE_FREE(POINTER_TO_DATA)    \
        {                                 \
            if ((POINTER_TO_DATA))        \
            {                             \
                free(POINTER_TO_DATA);    \
                (POINTER_TO_DATA) = NULL; \
            }                             \
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
    #define SAFE_FREE(POINTER_TO_DATA)
    #define STATIC_ASSERT(CONDITION, MSG)
    #define LOG(MSG, ...)
#endif

#endif // __UTILS_H__