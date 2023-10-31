#ifndef __GRAFIKA_H__
#define __GRAFIKA_H__

#include <stdbool.h>
#include <assert.h>

#include "shrifty.h"
#include "../utils.h"

#include "SDL2/SDL.h"

#define GRAFIKA_SCREEN_WIDTH  512
#define GRAFIKA_SCREEN_HEIGHT 512
#define GRAFIKA_BPP           4

typedef struct renderer
{
    bool          quit;
    SDL_Window   *window;
    SDL_Renderer *renderer;
    SDL_Texture  *texture;
    uint32_t     *pixels;
    // ALIGN_ME(64)
    // float  depth_buffer[GRAFIKA_SCREEN_WIDTH * GRAFIKA_SCREEN_HEIGHT];
    float *depth_buffer;
} Renderer_t;

static Renderer_t rend = {0};

static void grafika_present(void)
{
    void *px    = NULL;
    int   pitch = 0;
    SDL_LockTexture(rend.texture, NULL, &px, &pitch);
    {
        memcpy_s(px, GRAFIKA_SCREEN_HEIGHT * pitch,
                 rend.pixels, GRAFIKA_SCREEN_WIDTH * GRAFIKA_SCREEN_HEIGHT * sizeof(uint32_t));
    }
    SDL_UnlockTexture(rend.texture);
    SDL_RenderClear(rend.renderer);

    // Convert the SDL_Surface to an SDL_Texture
    SDL_Texture *surfaceTexture = SDL_CreateTextureFromSurface(rend.renderer, text_state.surface);

    SDL_RenderCopyEx(rend.renderer,
                     rend.texture,
                     NULL,
                     NULL,
                     0.0,
                     NULL,
                     SDL_FLIP_VERTICAL);

    SDL_RenderCopy(rend.renderer, surfaceTexture, NULL, NULL);

    // Clear the entire surface with black color
    SDL_FillRect(text_state.surface, NULL, SDL_MapRGBA(text_state.surface->format, 0, 0, 0, 0));

    SDL_RenderPresent(rend.renderer);

    SDL_DestroyTexture(surfaceTexture);
}

static inline void grafika_setpixel(int x, int y, uint32_t colour)
{
    ASSERT(y <= GRAFIKA_SCREEN_WIDTH, "y - %d out of bounds\n", y);
    ASSERT(x <= GRAFIKA_SCREEN_HEIGHT, "x - %d out of bounds\n", x);
    rend.pixels[(y * GRAFIKA_SCREEN_WIDTH) + x] = colour;
}

static inline void grafika_clear(void)
{
    memset(rend.pixels, 0, GRAFIKA_SCREEN_WIDTH * GRAFIKA_SCREEN_HEIGHT * GRAFIKA_BPP); // clear pixels

    // float *end = &rend.depth_buffer[GRAFIKA_SCREEN_WIDTH * GRAFIKA_SCREEN_HEIGHT];
    // for (float *p = &rend.depth_buffer[0]; p != end; p++) // clear depth buffer
    //     *p = 10.0f;

    for (size_t i = 0; i < GRAFIKA_SCREEN_WIDTH * GRAFIKA_SCREEN_HEIGHT; i++)
        rend.depth_buffer[i] = 10.0f;
}

static void grafika_startup(void)
{
    if (0 != SDL_Init(SDL_INIT_VIDEO))
        fprintf(stderr, "SDL failed to initialize: %s\n", SDL_GetError()), abort();

    rend.window =
        SDL_CreateWindow(
            "Title",
            SDL_WINDOWPOS_CENTERED_DISPLAY(0),
            SDL_WINDOWPOS_CENTERED_DISPLAY(0),
            512,
            512,
            0);
    ASSERT(rend.window, "Error - SDL_CreateWindow: %s\n", SDL_GetError());

    rend.renderer =
        SDL_CreateRenderer(
            rend.window,
            -1,
            SDL_RENDERER_ACCELERATED | SDL_RENDERER_PRESENTVSYNC);
    ASSERT(rend.renderer, "Error - SDL_CreateRenderer: %s\n", SDL_GetError());

    rend.texture =
        SDL_CreateTexture(
            rend.renderer,
            SDL_PIXELFORMAT_ABGR8888,
            SDL_TEXTUREACCESS_STREAMING,
            GRAFIKA_SCREEN_WIDTH,
            GRAFIKA_SCREEN_HEIGHT);
    ASSERT(rend.texture, "Error - SDL_CreateTexture: %s\n", SDL_GetError());

    // rend.pixels = malloc(sizeof(uint32_t) * GRAFIKA_SCREEN_WIDTH * GRAFIKA_SCREEN_HEIGHT);
    rend.pixels = _aligned_malloc(sizeof(uint32_t) * GRAFIKA_SCREEN_WIDTH * GRAFIKA_SCREEN_HEIGHT, 32);
    ASSERT(rend.pixels, "Error allocating pixel buffer\n");

    rend.depth_buffer = _aligned_malloc(sizeof(float) * GRAFIKA_SCREEN_WIDTH * GRAFIKA_SCREEN_HEIGHT, 64);
    ASSERT(rend.depth_buffer, "Error allocating depth buffer\n");

    SDL_SetRenderTarget(rend.renderer, NULL);
    SDL_SetRenderDrawColor(rend.renderer, 0, 0, 0, 0xFF);
    SDL_SetRenderDrawBlendMode(rend.renderer, SDL_BLENDMODE_BLEND); // enable alpha blending
}

static void grafika_shutdown(void)
{
    if (rend.texture)
        SDL_DestroyTexture(rend.texture), rend.texture = NULL;
    if (rend.renderer)
        SDL_DestroyRenderer(rend.renderer), rend.renderer = NULL;
    if (rend.window)
        SDL_DestroyWindow(rend.window), rend.window = NULL;
    SDL_Quit();

    if (rend.pixels)
        _aligned_free(rend.pixels), rend.pixels = NULL;
    if (rend.depth_buffer)
        _aligned_free(rend.depth_buffer), rend.depth_buffer = NULL;
}

#endif // __GRAFIKA_H__