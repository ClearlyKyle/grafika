#ifndef __MATEMATIKA_H__
#define __MATEMATIKA_H__

#include <math.h>

#if 0
#    define LH_COORDINATE_SYSTEM
#else
#    define RH_COORDINATE_SYSTEM /* Blender uses RH */
#endif

// GCC have funny inline rules, this will help
#if defined(_MSC_VER)
#    define _INLINE __forceinline
#else
#    define _INLINE static inline
#endif

#define max(a, b) (((a) > (b)) ? (a) : (b))
#define min(a, b) (((a) < (b)) ? (a) : (b))

#define _PI          3.14159265358979323846264338327950288 /* pi */
#define _PIf         ((float)_PI)
#define DEG2RAD(DEG) ((DEG)*_PIf / 180.0f)

typedef float vec3[3];
typedef float vec4[4];
typedef vec4  mat4[4];

_INLINE
float v3dot(const vec3 v0, const vec3 v1)
{
    return (v0[0] * v1[0] + v0[1] * v1[1] + v0[2] * v1[2]);
}

_INLINE
void v3cross(const vec3 v1, const vec3 v2, vec3 res)
{
    res[0] = v1[1] * v2[2] - v1[2] * v2[1];
    res[1] = v1[2] * v2[0] - v1[0] * v2[2];
    res[2] = v1[0] * v2[1] - v1[1] * v2[0];
}

_INLINE
void v3add(const vec3 v1, const vec3 v2, vec3 res)
{
    res[0] = v1[0] + v2[0];
    res[1] = v1[1] + v2[1];
    res[2] = v1[2] + v2[2];
}

_INLINE
void v3sub(const vec3 v0, const vec3 v1, vec3 res)
{
    res[0] = v0[0] - v1[0];
    res[1] = v0[1] - v1[1];
    res[2] = v0[2] - v1[2];
}

_INLINE
void v3div(const vec3 v, const float val, vec3 res)
{
    res[0] = v[0] / val;
    res[1] = v[1] / val;
    res[2] = v[2] / val;
}

_INLINE
void m4mulm4(const mat4 m1, const mat4 m2, mat4 dest)
{
    float a00 = m1[0][0], a01 = m1[0][1], a02 = m1[0][2], a03 = m1[0][3],
          a10 = m1[1][0], a11 = m1[1][1], a12 = m1[1][2], a13 = m1[1][3],
          a20 = m1[2][0], a21 = m1[2][1], a22 = m1[2][2], a23 = m1[2][3],
          a30 = m1[3][0], a31 = m1[3][1], a32 = m1[3][2], a33 = m1[3][3],

          b00 = m2[0][0], b01 = m2[0][1], b02 = m2[0][2], b03 = m2[0][3],
          b10 = m2[1][0], b11 = m2[1][1], b12 = m2[1][2], b13 = m2[1][3],
          b20 = m2[2][0], b21 = m2[2][1], b22 = m2[2][2], b23 = m2[2][3],
          b30 = m2[3][0], b31 = m2[3][1], b32 = m2[3][2], b33 = m2[3][3];

    dest[0][0] = a00 * b00 + a10 * b01 + a20 * b02 + a30 * b03;
    dest[0][1] = a01 * b00 + a11 * b01 + a21 * b02 + a31 * b03;
    dest[0][2] = a02 * b00 + a12 * b01 + a22 * b02 + a32 * b03;
    dest[0][3] = a03 * b00 + a13 * b01 + a23 * b02 + a33 * b03;
    dest[1][0] = a00 * b10 + a10 * b11 + a20 * b12 + a30 * b13;
    dest[1][1] = a01 * b10 + a11 * b11 + a21 * b12 + a31 * b13;
    dest[1][2] = a02 * b10 + a12 * b11 + a22 * b12 + a32 * b13;
    dest[1][3] = a03 * b10 + a13 * b11 + a23 * b12 + a33 * b13;
    dest[2][0] = a00 * b20 + a10 * b21 + a20 * b22 + a30 * b23;
    dest[2][1] = a01 * b20 + a11 * b21 + a21 * b22 + a31 * b23;
    dest[2][2] = a02 * b20 + a12 * b21 + a22 * b22 + a32 * b23;
    dest[2][3] = a03 * b20 + a13 * b21 + a23 * b22 + a33 * b23;
    dest[3][0] = a00 * b30 + a10 * b31 + a20 * b32 + a30 * b33;
    dest[3][1] = a01 * b30 + a11 * b31 + a21 * b32 + a31 * b33;
    dest[3][2] = a02 * b30 + a12 * b31 + a22 * b32 + a32 * b33;
    dest[3][3] = a03 * b30 + a13 * b31 + a23 * b32 + a33 * b33;
}

_INLINE
void m4mulv4(const mat4 m, const vec4 v, vec4 res)
{
    res[0] = m[0][0] * v[0] + m[1][0] * v[1] + m[2][0] * v[2] + m[3][0] * v[3];
    res[1] = m[0][1] * v[0] + m[1][1] * v[1] + m[2][1] * v[2] + m[3][1] * v[3];
    res[2] = m[0][2] * v[0] + m[1][2] * v[1] + m[2][2] * v[2] + m[3][2] * v[3];
    res[3] = m[0][3] * v[0] + m[1][3] * v[1] + m[2][3] * v[2] + m[3][3] * v[3];
}

_INLINE
void v3norm(vec3 v)
{
    float length = sqrtf(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);

    if (length == 0.0f)
    {
        v[0] = v[1] = v[2] = 0.0f;
        return;
    }

    v[0] /= length;
    v[1] /= length;
    v[2] /= length;
}

// mat4 functions ---------------------------------------------------------------------------------

_INLINE
void m4identity(mat4 m)
{
    memset(m, 0, sizeof(mat4));
    m[0][0] = m[1][1] = m[2][2] = m[3][3] = 1.0f;
}

_INLINE
void m4_lookat(const vec3 eye, const vec3 center, const vec3 up, mat4 res)
{
    vec3 forward;
    v3_sub(center, eye, forward);
    v3_norm(forward);

#ifdef RH_COORDINATE_SYSTEM // Right handed
    vec3 right;
    v3_cross(forward, up, right);
    v3_norm(right);

    vec3 upp;
    v3_cross(right, forward, upp);

    res[0][0] = right[0];
    res[0][1] = upp[0];
    res[0][2] = -forward[0];
    res[0][3] = 0.0f;

    res[1][0] = right[1];
    res[1][1] = upp[1];
    res[1][2] = -forward[1];
    res[1][3] = 0.0f;

    res[2][0] = right[2];
    res[2][1] = upp[2];
    res[2][2] = -forward[2];
    res[2][3] = 0.0f;

    res[3][0] = -v3_dot(right, eye);
    res[3][1] = -v3_dot(upp, eye);
    res[3][2] = v3_dot(forward, eye);
    res[3][3] = 1.0f;
#else // Left handed
    vec3 left;
    v3_cross(up, forward, left);
    v3_norm(left);

    vec3 upp;
    v3_cross(forward, left, upp);

    res[0][0] = left[0];
    res[0][1] = upp[0];
    res[0][2] = forward[0];
    res[0][3] = 0.0f;

    res[1][0] = left[1];
    res[1][1] = upp[1];
    res[1][2] = forward[1];
    res[1][3] = 0.0f;

    res[2][0] = left[2];
    res[2][1] = upp[2];
    res[2][2] = forward[2];
    res[2][3] = 0.0f;

    res[3][0] = -v3_dot(left, eye);
    res[3][1] = -v3_dot(upp, eye);
    res[3][2] = -v3_dot(forward, eye);
    res[3][3] = 1.0f;
#endif
}

_INLINE
void m4_proj(float fovy, float aspect, float nearZ, float farZ, mat4 res)
{
    memset(res, 0, sizeof(mat4));

    const float f  = 1.0f / tanf(fovy * 0.5f);
    const float fn = 1.0f / (nearZ - farZ);

#ifdef RH_COORDINATE_SYSTEM // RH [-1, 1]
    res[0][0] = f / aspect;
    res[1][1] = f;
    res[2][2] = (nearZ + farZ) * fn;
    res[2][3] = -1.0f;
    res[3][2] = 2.0f * nearZ * farZ * fn;
#else // LH [-1, 1]
    res[0][0] = f / aspect;
    res[1][1] = f;
    res[2][2] = -(nearZ + farZ) * fn;
    res[2][3] = 1.0f;
    res[3][2] = 2.0f * nearZ * farZ * fn;
#endif
}

_INLINE
void m4transmake(float x, float y, float z, mat4 res)
{
    m4identity(res);
    res[3][0] = x;
    res[3][1] = y;
    res[3][2] = z;
}

_INLINE
void makeAABB(vec3 pos[3], int AABB[4])
{
    /* Get the bounding box of the triangle,
        setting pos[0] as starting values */
    float fminX = pos[0][0];
    float fminY = pos[0][1];
    float fmaxX = pos[0][0];
    float fmaxY = pos[0][1];

    for (int i = 1; i < 3; ++i)
    {
        /* Update minimum and maximum values for x and y */
        const float x = pos[i][0];
        const float y = pos[i][1];

        if (x < fminX)
            fminX = x;
        if (y < fminY)
            fminY = y;

        if (x > fmaxX)
            fmaxX = x;
        if (y > fmaxY)
            fmaxY = y;
    }
    /* Precompute constants */
    const int screenWidthMinusOne  = GRAFIKA_SCREEN_WIDTH - 1;
    const int screenHeightMinusOne = GRAFIKA_SCREEN_HEIGHT - 1;

    /* Clamp values to valid range */
    AABB[0] = max(0, min((int)fminX, screenWidthMinusOne));  // minX
    AABB[1] = max(0, min((int)fminY, screenHeightMinusOne)); // minY
    AABB[2] = max(0, min((int)fmaxX, screenWidthMinusOne));  // maxX
    AABB[3] = max(0, min((int)fmaxY, screenHeightMinusOne)); // maxY
}

#endif // __MATEMATIKA_H__