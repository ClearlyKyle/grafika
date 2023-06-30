#ifndef __TIMER_H__
#define __TIMER_H__

#include "SDL2/SDL.h"

typedef struct timer
{
    Uint64 start, elapsed, perf_frequency;
} timer_t;

static struct timer Timer_Init_Start(void)
{
    struct timer t   = {0};
    t.perf_frequency = SDL_GetPerformanceFrequency();
    t.start          = SDL_GetPerformanceCounter();
    t.elapsed        = 0;
    return t;
}

static void Timer_Update(struct timer *const timer)
{
    const Uint64 new_time = SDL_GetPerformanceCounter();
    timer->elapsed        = new_time - timer->start;
    timer->start          = new_time;
}

static double Timer_Get_Elapsed_Seconds(const struct timer *const timer)
{
    return (double)timer->elapsed / (double)timer->perf_frequency;
}

static double Timer_Get_Elapsed_MS(const struct timer *const timer)
{
    return Timer_Get_Elapsed_Seconds(timer) * 1000.0;
}

#endif // __TIMER_H__