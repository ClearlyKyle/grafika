#ifndef __UTILS_H__
#define __UTILS_H__

#define ASSERT(EXPR, ...)                                                   \
    if (!(EXPR))                                                            \
    {                                                                       \
        fprintf(stderr, "[ASSERT] %s %s:%d - ", #EXPR, __FILE__, __LINE__); \
        fprintf(stderr, __VA_ARGS__);                                       \
        abort();                                                            \
    }

#endif // __UTILS_H__