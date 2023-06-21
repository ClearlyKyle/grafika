#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <stdbool.h>

static inline void verline(int x, int y0, int y1, uint32_t colour)
{
    for (int y = y0; y <= y1; y++)
    {
        grafika_setpixel(x, y, colour);
    }
}

void draw_square(int x, int y, int size, uint32_t colour)
{
    int i = 0;

    for (i = x; i < x + size; i++) // Draw the top side
        grafika_setpixel(i, y, colour);

    for (i = y; i < y + size; i++) // Draw the right side
        grafika_setpixel(x + size - 1, i, colour);

    for (i = x; i < x + size; i++) // Draw the bottom side
        grafika_setpixel(i, y + size - 1, colour);

    for (i = y; i < y + size; i++) // Draw the left side
        grafika_setpixel(x, i, colour);
}

static void update(void)
{
    draw_square(1, 1, 300, 0xFF223311);
}

int main(int argc, char *argv[])
{
    grafika_startup();

    while (!rend.quit)
    {
        SDL_Event ev;
        while (SDL_PollEvent(&ev))
        {
            switch (ev.type)
            {
            case SDL_QUIT:
                rend.quit = true;
                break;
            default:
                break;
            }
        }

        grafika_clear();
        update();
        grafika_present();
    }

    grafika_destroy();

    return 0;
}