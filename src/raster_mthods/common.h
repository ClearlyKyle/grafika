#ifndef __COMMON_H__
#define __COMMON_H__

#include "../grafika.h"
#include "../tex.h"
#include "../obj.h"
#include "../timer.h"
#include "../matematika.h"

typedef vec3 triangle[3];

typedef struct rasterstate
{
    mat4 proj;
    mat4 view;
    mat4 model;

    obj_t obj;
    tex_t tex;

    // NOTE : Should this be with the graphics?
} rasterstate_t;

static float zbuffer[GRAFIKA_SCREEN_HEIGHT * GRAFIKA_SCREEN_WIDTH];

static rasterstate_t state = {0};

#endif // __COMMON_H__