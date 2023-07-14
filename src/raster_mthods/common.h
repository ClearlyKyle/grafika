#ifndef __COMMON_H__
#define __COMMON_H__

#include "../grafika.h"
#include "../tex.h"
#include "../obj.h"
#include "../timer.h"
#include "../matematika.h"

typedef struct triangle
{
    vec3 pos[3];
    vec2 tex[3];
} triangle_t;

typedef struct rasterstate
{
    obj_t obj;
    tex_t tex;
    mat4  MVP;
} rasterstate_t;

static rasterstate_t state = {0};

#endif // __COMMON_H__