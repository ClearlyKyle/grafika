#ifndef __GRAFIKA_H__
#define __GRAFIKA_H__

#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>

#include "utils.h"

#include "SDL2/SDL.h"
#include "SDL2/SDL_image.h"

#ifndef GRAFIKA_SCREEN_WIDTH
#define GRAFIKA_SCREEN_WIDTH (512)
#endif
#ifndef GRAFIKA_SCREEN_HEIGHT
#define GRAFIKA_SCREEN_HEIGHT (512)
#endif
#ifndef GRAFIKA_TITLE
#define GRAFIKA_TITLE ("grafika")
#endif

struct grafika
{
    SDL_Window  *window;
    SDL_Surface *surface;
    uint32_t    *pixels;
    float       *depth_buffer;
};

static struct grafika rend = {0};

static inline void grafika_present(void)
{
    ASSERT(rend.surface, "Surface was not set in the renderer\n");

    SDL_UpdateWindowSurface(rend.window);
}

static inline void grafika_setpixel(uint32_t x, uint32_t y, uint32_t colour)
{
    ASSERT(y < GRAFIKA_SCREEN_HEIGHT, "y - %u out of bounds\n", y);
    ASSERT(x < GRAFIKA_SCREEN_WIDTH, "x - %u out of bounds\n", x);

#if IS_POWER_OF_2(GRAFIKA_SCREEN_WIDTH)
    rend.pixels[(y << 9) + x] = colour; // why 9? 1 << 9 == 512, y << 9 == y * 512
#else
    rend.pixels[(y * GRAFIKA_SCREEN_WIDTH) + x] = colour;
#endif
}

static inline void grafika_clear(void)
{
    memset(rend.pixels, 0, GRAFIKA_SCREEN_WIDTH * GRAFIKA_SCREEN_HEIGHT * rend.surface->format->BytesPerPixel);
    memset(rend.depth_buffer, 0x48, GRAFIKA_SCREEN_WIDTH * GRAFIKA_SCREEN_HEIGHT * sizeof(float));
    // why 0x48? floats are 4 bytes, 0x48 per byte, thus we get 0x4848 for 4 bytes, which is +8.57f
}

static void grafika_startup(struct arena *arena)
{
    if (0 != SDL_Init(SDL_INIT_VIDEO))
        LOGE("SDL failed to initialize: %s\n", SDL_GetError()), abort();

    int img_flags = IMG_INIT_JPG | IMG_INIT_PNG;
    if (!(IMG_Init(img_flags) & img_flags))
        LOGE("IMG_Init failed to initialize: %s\n", IMG_GetError()), abort();

    rend.window = SDL_CreateWindow(GRAFIKA_TITLE,
                                   SDL_WINDOWPOS_CENTERED,
                                   SDL_WINDOWPOS_CENTERED,
                                   GRAFIKA_SCREEN_WIDTH,
                                   GRAFIKA_SCREEN_HEIGHT,
                                   0);
    ASSERT(rend.window, "Error - SDL_CreateWindow: %s\n", SDL_GetError());

    rend.surface = SDL_GetWindowSurface(rend.window);
    ASSERT(rend.surface, "Error - SDL_GetWindowSurface: %s\n", SDL_GetError());

    rend.pixels = rend.surface->pixels;
    ASSERT(rend.pixels, "Error getting surface pixels\n");

    SDL_PixelFormat *pixel_format = rend.surface->format;
    const char      *format_name  = SDL_GetPixelFormatName(pixel_format->format);
    LOG("Window pixel format: %s, bpp: %d\n", format_name, pixel_format->BitsPerPixel);
    LOG("Shift: R: %d, G:%d, B:%d, A:%d\n", pixel_format->Rshift, pixel_format->Gshift, pixel_format->Bshift, pixel_format->Ashift);

    rend.depth_buffer = arena_alloc_aligned(arena, sizeof(float) * GRAFIKA_SCREEN_WIDTH * GRAFIKA_SCREEN_HEIGHT, 16);
    ASSERT(rend.depth_buffer, "Error allocating depth buffer\n");
}

static void grafika_shutdown(void)
{
    IMG_Quit();

    if (rend.surface) SDL_FreeSurface(rend.surface);
    if (rend.window) SDL_DestroyWindow(rend.window);

    // depth cleaned up by the arena

    SDL_Quit();

    LOG("grafika_shutdown\n");
}

#endif // __GRAFIKA_H__