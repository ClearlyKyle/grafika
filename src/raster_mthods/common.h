#ifndef __COMMON_H__
#define __COMMON_H__

#include "../grafika.h"
#include "../shrifty.h"
#include "../utils.h"
#include "../tex.h"
#include "../obj.h"
#include "../matematika.h"
#include "../bench.h"

#include <xmmintrin.h>
#include <immintrin.h>

struct triangle
{
    vec3 pos[3];
    vec2 tex[3];
};

struct rasterstate
{
    mat4       model;
    mat4       MVP;
    vec3       cam_pos;
    tex_t      tex;
    struct obj obj;
};

void draw_object(struct arena *arena);
void draw_onexit(void);

#endif // __COMMON_H__