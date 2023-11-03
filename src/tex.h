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

tex_t tex_load(const char *file_path)
{
    tex_t t = {0};

    SDL_Surface *surface = IMG_Load(file_path);
    ASSERT(surface, "Failed to load image : %s\n", file_path);

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