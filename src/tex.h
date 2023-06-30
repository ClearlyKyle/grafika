#ifndef __TEX_H__
#define __TEX_H__

#include <stdbool.h>

#include "utils.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

typedef struct tex
{
    int            w, h, bpp;
    unsigned char *data;
} tex_t;

tex_t tex_load(const char *file_path, bool flip)
{
    tex_t t = {0};

    stbi_set_flip_vertically_on_load(flip ? 1 : 0);

    t.data = stbi_load(file_path, &t.w, &t.h, &t.bpp, 0);
    ASSERT(t.data, "Failed to open file : '%s'\n", file_path);

    return t;
}

void tex_destroy(tex_t *t)
{
    if (t->data)
    {
        free(t->data);
        t->data = NULL;
    }
    *t = (tex_t){0};
}
#endif // __TEX_H__