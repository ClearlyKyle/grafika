#ifndef __SHRIFTY_H__ // Fonts
#define __SHRIFTY_H__

#include "SDL2/SDL_ttf.h"
#include "SDL2/SDL.h"

#include "utils.h"

typedef struct text
{
    TTF_Font    *font;
    SDL_Surface *surface;
} text_t;

static char   shrifty_text_buffer[128] = {0}; // should this be static?
static text_t text_state               = {0};

static inline void text_startup(int w, int h, int font_size)
{
    if (0 != TTF_Init())
        fprintf(stderr, "TTF failed to initialize: %s\n", TTF_GetError()), abort();

    // surface to hold all written text
    text_state.surface = SDL_CreateRGBSurface(0, w, h, 32, 0x000000FF, 0x0000FF00, 0x00FF0000, 0xFF000000);
    ASSERT(text_state.surface, "Failed to create surface SDL_CreateRGBSurface : %s\n", SDL_GetError());

    text_state.font = TTF_OpenFont("res/font/Arial.ttf", font_size);
    ASSERT(text_state.font, "Failed to load the font : %s\n", TTF_GetError());
}

static inline void text_shutdown(void)
{
    if (text_state.surface)
        SDL_FreeSurface(text_state.surface);
    if (text_state.font)
        TTF_CloseFont(text_state.font);
    if (TTF_WasInit())
        TTF_Quit();
}

static void text_write(int x, int y, const char *formatted_text, ...) ATTRIBATE_FORMAT_PRINTF(3, 4);
static void text_write(int x, int y, const char *formatted_text, ...)
{
    ASSERT(formatted_text, "Text is NULL\n");

    va_list lst;
    va_start(lst, formatted_text);
    vsnprintf(shrifty_text_buffer, sizeof(shrifty_text_buffer), formatted_text, lst);
    va_end(lst);

    SDL_Color    text_colour  = {255, 255, 255, 255}; // White color for text
    SDL_Surface *text_surface = TTF_RenderText_Blended(text_state.font, shrifty_text_buffer, text_colour);
    ASSERT(text_surface, "TTF_RenderText_Blended FAILED:\nError - %s\ntext  - %s\n", TTF_GetError(), formatted_text);

    SDL_Rect text_dest = {x, y, text_surface->w, text_surface->h};

    SDL_BlitSurface(text_surface, NULL, text_state.surface, &text_dest);

    if (text_surface)
        SDL_FreeSurface(text_surface);
}

#endif // __SHRIFTY_H__