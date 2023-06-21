#ifndef __OBJ_H__
#define __OBJ_H__

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

typedef struct vertindices
{
    int v_idx, vn_idx, vt_idx;
} vertindices_t;

typedef struct obj
{
    float *pos;
    float *norms;
    float *texs;

    vertindices_t *indices;
    size_t         num_verts;

    size_t num_pos;
    size_t num_norms;
    size_t num_texs;
    size_t num_f_rows;
} obj_t;

obj_t obj_load(const char *filename)
{
    FILE *fp = NULL;
    fopen_s(&fp, filename, "r");
    assert(fp);

    char linebuffer[256]; // Adjust the size

    int posCount    = 0;
    int texCount    = 0;
    int normalCount = 0;
    int frowCount   = 0;
    while (fgets(linebuffer, sizeof(linebuffer), fp) != NULL)
    {
        if (linebuffer[0] == 'v' && linebuffer[1] == 'n')
            normalCount++;
        else if (linebuffer[0] == 'v' && linebuffer[1] == 't')
            texCount++;
        else if (linebuffer[0] == 'v')
            posCount++;
        else if (linebuffer[0] == 'f')
            frowCount++;
    }

    obj_t obj      = {0};
    obj.num_norms  = normalCount;
    obj.num_texs   = texCount;
    obj.num_pos    = posCount;
    obj.num_f_rows = frowCount;
    obj.num_verts  = frowCount * 3; // Assume triangulation

    obj.pos   = malloc(sizeof(float) * posCount * 3);
    obj.norms = malloc(sizeof(float) * normalCount * 3);
    obj.texs  = malloc(sizeof(float) * texCount * 2);

    obj.indices = malloc(sizeof(vertindices_t) * frowCount * 3);

    // Reset file pointer to the start of the file
    fseek(fp, 0, SEEK_SET);

    int nrm_idx = 0, pos_idx = 0, tex_idx = 0, indi_idx = 0;
    while (fgets(linebuffer, sizeof(linebuffer), fp) != NULL)
    {
        if (linebuffer[0] == 'v' && linebuffer[1] == 'n')
        {
            sscanf(linebuffer, "vn %f %f %f",
                   &obj.norms[nrm_idx + 0],
                   &obj.norms[nrm_idx + 1],
                   &obj.norms[nrm_idx + 2]);
            nrm_idx += 3;
        }
        else if (linebuffer[0] == 'v' && linebuffer[1] == 't')
        {
            sscanf(linebuffer, "vt %f %f",
                   &obj.texs[tex_idx + 0],
                   &obj.texs[tex_idx + 1]);
            tex_idx += 2;
        }
        else if (linebuffer[0] == 'v')
        {
            sscanf(linebuffer, "v %f %f %f",
                   &obj.pos[pos_idx + 0],
                   &obj.pos[pos_idx + 1],
                   &obj.pos[pos_idx + 2]);
            pos_idx += 3;
        }
        else if (linebuffer[0] == 'f')
        {
            int fv[3], fvt[3], fvn[3];
            int result = sscanf(linebuffer, "f %d/%d/%d %d/%d/%d %d/%d/%d",
                                &fv[0], &fvt[0], &fvn[0],
                                &fv[1], &fvt[1], &fvn[1],
                                &fv[2], &fvt[2], &fvn[2]);

            if (result != 9)
                printf("OBJ ERROR! the face values are of a new format, update this code!\n");

            for (size_t f = 0; f < 3; f++)
            {
                vertindices_t *fptr = &obj.indices[indi_idx + f];
                fptr->v_idx         = fv[f] - 1;
                fptr->vt_idx        = fvt[f] - 1;
                fptr->vn_idx        = fvn[f] - 1;
            }
            indi_idx += 3;
        }
    }

    // printf("poss : %d\n", posCount);
    // printf("texs : %d\n", texCount);
    // printf("norms: %d\n", normalCount);
    // printf("faces: %d\n", frowCount);
    // printf("verts: %d\n", frowCount * 3);

    // for (size_t i = 0; i < obj.num_pos; i++)
    //{
    //     const int index = i * 3;
    //     printf("v %f, %f, %f\n",
    //            obj.pos[index + 0],
    //            obj.pos[index + 1],
    //            obj.pos[index + 2]);
    // }

    // for (size_t i = 0; i < obj.num_norms; i++)
    //{
    //     const int index = i * 3;
    //     printf("vn %f, %f, %f\n",
    //            obj.norms[index + 0],
    //            obj.norms[index + 1],
    //            obj.norms[index + 2]);
    // }

    // for (size_t i = 0; i < obj.num_texs; i++)
    //{
    //     const int index = i * 2;
    //     printf("vt %f, %f\n",
    //            obj.texs[index + 0],
    //            obj.texs[index + 1]);
    // }

    // for (size_t i = 0; i < obj.num_verts; i++)
    //{
    //     printf("%d, %d, %d\n",
    //            obj.indices[i].v_idx,
    //            obj.indices[i].vt_idx,
    //            obj.indices[i].vn_idx);
    // }

    fclose(fp);

    return obj;
}

void obj_destroy(obj_t *obj)
{
    if (!obj)
        return;

    if (obj->pos)
        free(obj->pos);
    if (obj->norms)
        free(obj->norms);
    if (obj->texs)
        free(obj->texs);
    if (obj->indices)
        free(obj->indices);

    *obj = (obj_t){0};

    printf("obj destroyed\n");
}

#endif // __OBJ_H__