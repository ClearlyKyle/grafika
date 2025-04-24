#ifndef __SHRIFTY_H__
#define __SHRIFTY_H__

#include "SDL2/SDL_ttf.h"
#include "SDL2/SDL.h"

#include "utils.h"

#define SHRIFTY_GLYPH_CACHE_COUNT (128)
struct glyph
{
    SDL_Surface *surface;
    int          w, h;
};

struct shrifty
{
    TTF_Font    *font;
    SDL_Surface *surface;

    // struct glyph glyph_cache[SHRIFTY_GLYPH_CACHE_COUNT];
};

static struct glyph   shrifty_glyph_cache[SHRIFTY_GLYPH_CACHE_COUNT];
static char           shrifty_text_buffer[128];
static struct shrifty shrifty_state;

static void text_startup(SDL_Surface *surface, int font_size)
{
    if (0 != TTF_Init()) LOGE("TTF failed to initialize: %s\n", TTF_GetError()), abort();

    shrifty_state.font = TTF_OpenFont("c:\\WINDOWS\\Fonts\\COUR.TTF", font_size);
    ASSERT(shrifty_state.font, "Failed to load the font : %s\n", TTF_GetError());

    TTF_SetFontHinting(shrifty_state.font, TTF_HINTING_LIGHT);
    TTF_SetFontKerning(shrifty_state.font, 0);

    shrifty_state.surface = surface;

    {
        int minx, maxx, miny, maxy, advance;
        TTF_GlyphMetrics(shrifty_state.font, 'a', &minx, &maxx, &miny, &maxy, &advance);
        LOG("Metrics : %c - %d, %d, %d, %d, %d\n", 'a', minx, maxx, miny, maxy, advance);
    }
    {
        int minx, maxx, miny, maxy, advance;
        TTF_GlyphMetrics(shrifty_state.font, 'm', &minx, &maxx, &miny, &maxy, &advance);
        LOG("Metrics : %c - %d, %d, %d, %d, %d\n", 'm', minx, maxx, miny, maxy, advance);
    }
}

static void text_write(int x, int y, const char *formatted_text, ...) ATTRIBATE_FORMAT_PRINTF(3, 4);
static void text_write(int x, int y, const char *formatted_text, ...)
{
    ASSERT(shrifty_state.font, "Font state is NULL\n");
    ASSERT(formatted_text, "Text is NULL\n");

    va_list lst;
    va_start(lst, formatted_text);
    int needed = vsnprintf(shrifty_text_buffer, sizeof(shrifty_text_buffer), formatted_text, lst);
    if (needed >= 128) LOGE("Writing outside of shrifty_text_buffer\n");
    va_end(lst);

    for (char *c = shrifty_text_buffer; *c; c++)
    {
        ASSERT((*c) > 0, "c is out of bounds for a uint16_t");
        // ASSERT((*c) < SHRIFTY_GLYPH_CACHE_COUNT, "character cache is too small");

        Uint16 character = (Uint16)(*c);

        struct glyph glyph = {0};

        if (!(shrifty_glyph_cache[character].surface))
        {
            SDL_Surface *surf = TTF_RenderGlyph_Blended(shrifty_state.font, character, (SDL_Color){255, 255, 255, 255});
            if (surf)
            {
                ASSERT(surf, "Failed to create surface for character '%c'", *c);

                int advance; // TODO : set this for all? so can only use monospaced fonts
                TTF_GlyphMetrics(shrifty_state.font, character, NULL, NULL, NULL, NULL, &advance);

                glyph = (struct glyph){.surface = surf,
                                       .w       = advance,
                                       .h       = surf->h};

                shrifty_glyph_cache[character] = glyph;

                LOG("Chaching : %c - width : %d\n", *c, advance);
            }
        }
        else
        {
            glyph = shrifty_glyph_cache[character];
        }

        SDL_Rect dest = {x, y, glyph.w, glyph.h};
        SDL_UpperBlit(glyph.surface, NULL, shrifty_state.surface, &dest);

        x += glyph.w; // Advance cursor
    }
}

void text_shutdown(void)
{
    if (shrifty_state.font) TTF_CloseFont(shrifty_state.font);
    if (TTF_WasInit()) TTF_Quit();

    // for (size_t i = 0; i < SHRIFTY_GLYPH_CACHE_COUNT; i++)
    //{
    //     if (shrifty_state.glyph_cache[i].surface)
    //     {
    //         SDL_FreeSurface(shrifty_state.glyph_cache[i].surface);
    //         shrifty_state.glyph_cache[i].surface = NULL;
    //     }
    // }

    LOG("text_shutdown\n");
}

#endif // __SHRIFTY_H__
