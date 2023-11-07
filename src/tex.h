#ifndef __TEX_H__
#define __TEX_H__

#include <stdbool.h>

#include "utils.h"
#include "SDL2/SDL_image.h"

typedef struct tex
{
    int            w, h, bpp;
    unsigned char *data;
    SDL_Surface   *surface;
} tex_t;

static void _flip_surface(SDL_Surface *surface)
{
    SDL_LockSurface(surface);

#if 0
    const int pitch = surface->pitch; // row size

    ASSERT(pitch > 0, "Error with pitch of the surface (pitch:%d)\n", pitch);
    char *temp   = malloc(sizeof(char) * (size_t)pitch); // intermediate buffer
    char *pixels = (char *)surface->pixels;

    for (int i = 0; i < surface->h / 2; ++i)
    {
        // get pointers to the two rows to swap
        char *row1 = pixels + i * pitch;
        char *row2 = pixels + (surface->h - i - 1) * pitch;

        // swap rows
        memcpy(temp, row1, (const size_t)pitch);
        memcpy(row1, row2, (const size_t)pitch);
        memcpy(row2, temp, (const size_t)pitch);
    }

    free(temp);
#else
    const int pitch = surface->pitch; // row size

    ASSERT(pitch > 0, "Error with pitch of the surface (pitch:%d)\n", pitch);
    char *pixels = (char *)surface->pixels;

    for (int i = 0; i < surface->h / 2; ++i)
    {
        // Get pointers to the two rows to swap
        char *row1 = pixels + i * pitch;
        char *row2 = pixels + (surface->h - i - 1) * pitch;

        // Swap rows using SSE2 intrinsics
        for (int j = 0; j < pitch; j += 16) // Assuming 16 bytes (128 bits) can be loaded and stored at once
        {
            __m128i xmm_row1 = _mm_load_si128((__m128i *)(row1 + j));
            __m128i xmm_row2 = _mm_load_si128((__m128i *)(row2 + j));
            _mm_store_si128((__m128i *)(row1 + j), xmm_row2);
            _mm_store_si128((__m128i *)(row2 + j), xmm_row1);
        }
    }
#endif

    SDL_UnlockSurface(surface);
}

tex_t tex_load(const char *file_path)
{
    tex_t t = {0};

    SDL_Surface *surface = IMG_Load(file_path);
    ASSERT(surface, "Failed to load image\n\t|%s|\n", file_path);

    _flip_surface(surface);

    t.surface = surface;
    t.h       = surface->h;
    t.w       = surface->w;
    t.bpp     = surface->format->BytesPerPixel;
    t.data    = surface->pixels;

    LOG("IMAGE bpp(%d) w(%d) h(%d)\t%s\n", t.bpp, t.w, t.h, file_path);

    return t;
}

void tex_destroy(tex_t *t)
{
    ASSERT(t, "tex_destroy - texture is not valid\n");

    SDL_FreeSurface(t->surface), t->surface = NULL;
    *t = (tex_t){0};
}

#endif // __TEX_H__