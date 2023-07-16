#ifndef __COMMON_H__
#define __COMMON_H__

#include "../grafika.h"
#include "../../tex.h"
#include "../../obj.h"
#include "../../timer.h"
#include "../../matematika.h"

typedef struct triangle
{
    vec3 pos[3];
    vec2 tex[3];
} triangle_t;

typedef struct rasterstate
{
    mat4  model;
    mat4  MVP;
    tex_t tex;
    obj_t obj;
} rasterstate_t;

static rasterstate_t state = {0};

#endif // __COMMON_H__