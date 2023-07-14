#ifndef __UTILS_H__
#define __UTILS_H__

#ifndef NDEBUG
#    define ASSERT(EXPR, ...)                                                   \
        if (!(EXPR))                                                            \
        {                                                                       \
            fprintf(stderr, "[ASSERT] %s %s:%d - ", #EXPR, __FILE__, __LINE__); \
            fprintf(stderr, __VA_ARGS__);                                       \
            abort();                                                            \
        }
#else
#    define ASSERT(EXPR, ...)
#endif

#endif // __UTILS_H__