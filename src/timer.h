#ifndef __TIMER_H__
#define __TIMER_H__

#include "SDL2/SDL.h"

typedef struct timer
{
    Uint64 start, elapsed, perf_frequency;
} timer_t;

#define TIMER_START(TIMER)                                      \
    {                                                           \
        (TIMER).perf_frequency = SDL_GetPerformanceFrequency(); \
        (TIMER).start          = SDL_GetPerformanceCounter();   \
        (TIMER).elapsed        = 0;                             \
    }

#define TIMER_UPDATE(TIMER)                                  \
    {                                                        \
        const Uint64 new_time = SDL_GetPerformanceCounter(); \
        (TIMER).elapsed       = new_time - (TIMER).start;    \
        (TIMER).start         = new_time;                    \
    }

#define TIMER_ELAPSED_S(TIMER) ((double)(TIMER).elapsed / (double)(TIMER).perf_frequency)

#define TIMER_ELAPSED_MS(TIMER) (TIMER_ELAPSED_S(TIMER) * 1000.0)

#endif // __TIMER_H__