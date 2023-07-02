#ifndef __OBJ_H__
#define __OBJ_H__

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "utils.h"

#define MAX_CHAR_LINE_BUFFER 256

typedef struct material
{
    char *newmtl;
    float Ns;
    float Ka[3];
    float Kd[3];
    float Ks[3];
    float Ke[3];
    float Ni;
    float d;
    int   illum;

    char *map_Bump;
    char *map_Kd;
    char *map_Ns;
    char *refl;
} material_t;

typedef struct vertindices
{
    int v_idx, vn_idx, vt_idx;
} vertindices_t;

typedef struct BoundingBox
{
    float min[3], max[3]; // x, y, z
} BoundingBox_t;

typedef struct obj
{
    size_t num_pos; // X, Y, Z of vertex
    size_t num_norms;
    size_t num_texs;
    size_t num_f_rows; // How many rows starting with f
    size_t num_verts;  // How many individual vertices

    BoundingBox_t bbox;

    float *pos;
    float *norms;
    float *texs;

    vertindices_t *indices; // Index values for 1 vertex
} obj_t;

static char *_read_string_until(char *inputString, const char stopChar)
{
    // TODO : this should take in some kind of length
    assert(inputString);

    while (*inputString != '\0' && *inputString != stopChar)
        inputString++;

    return inputString;
}

static void removeNewline(char *str)
{
    assert(str);

    const size_t length = strcspn(str, "\n"); // TODO : remove this function call

    if (str[length] == '\n')
        str[length] = '\0'; // Replace newline with null-termination
}

static void _material_file(char *line)
{
    char *material_file_path = _read_string_until(line, ' ');

    removeNewline(++material_file_path); // NOTE : this might be unsafe

    printf("material file path : |%s|\n", material_file_path);

    FILE *fp = NULL;
    fopen_s(&fp, material_file_path, "r");

    char linebuffer[256]; // Adjust the size
    while (fgets(linebuffer, sizeof(linebuffer), fp) != NULL)
    {
        if (strncmp(linebuffer, "map_Bump", 8) == 0)
        {
            char *bump_path = _read_string_until(linebuffer, ' ');
            removeNewline(++bump_path); // NOTE : this might be unsafe
            printf("Bump map path: |%s|\n", bump_path);
        }
        else if (strncmp(linebuffer, "map_Kd", 6) == 0)
        {
            char *bump_path = _read_string_until(linebuffer, ' ');
            removeNewline(++bump_path); // NOTE : this might be unsafe
            printf("Bump map path: |%s|\n", bump_path);
        }
    }
    fclose(fp);
}

static void parse_f_line(char *fline, const int num_vertex_values, vertindices_t index_data[6])
{
    int   vertex_index, texture_index, normal_index;
    char *current_char = fline + 2; // Skip the leading 'f' character and whitespace

    int collected_indices = 0;

    while (*current_char != '\0')
    {
        vertex_index = atoi(current_char); // Extract vertex index

        // Find the position of the next '/' character
        while (*current_char != '/' && *current_char != ' ' && *current_char != '\0')
            current_char++;

        if (*current_char == '\0')
            break; // Reached the end of the line

        current_char++; // Move to the next character after '/'

        // Extract texture index
        texture_index = atoi(current_char);

        // Find the position of the next '/' character
        while (*current_char != '/' && *current_char != ' ' && *current_char != '\0')
            current_char++;

        if (*current_char == '\0')
            break; // Reached the end of the line

        current_char++; // Move to the next character after '/'

        // Extract normal index
        normal_index = atoi(current_char);

        // Find the position of the next space or null character
        while (*current_char != ' ' && *current_char != '\0')
            current_char++;

        if (*current_char == '\0')
            break; // Reached the end of the line

        current_char++; // Move to the next character after the space

        // Print the extracted values
        index_data[collected_indices].v_idx  = vertex_index - 1;
        index_data[collected_indices].vt_idx = texture_index - 1;
        index_data[collected_indices].vn_idx = normal_index - 1;
        collected_indices++;
    }
}

obj_t obj_load(const char *filename)
{
    FILE *fp = NULL;
    fopen_s(&fp, filename, "r");
    ASSERT(fp, "Error with opening object file : '%s'\n", filename);

    char linebuffer[MAX_CHAR_LINE_BUFFER]; // Adjust the size

    int posCount    = 0;
    int texCount    = 0;
    int normalCount = 0;
    int frowCount   = 0;
    int vertCount   = 0;
    while (fgets(linebuffer, sizeof(linebuffer), fp) != NULL)
    {
        if (linebuffer[0] == '#')
            continue;

        if (linebuffer[0] == 'v' && linebuffer[1] == 'n')
            normalCount++;
        else if (linebuffer[0] == 'v' && linebuffer[1] == 't')
            texCount++;
        else if (linebuffer[0] == 'v')
            posCount++;
        else if (linebuffer[0] == 'f')
        {
            // Count the number of spaces in the line
            int spaceCount = 0;
            for (int i = 0; i < strlen(linebuffer); i++)
                if (linebuffer[i] == ' ')
                    spaceCount++;

            if (spaceCount == 4)
            {
                vertCount += 3;
                frowCount++;
            }
            vertCount += 3;
            frowCount++;
        }
    }

    obj_t obj      = {0};
    obj.num_norms  = normalCount;
    obj.num_texs   = texCount;
    obj.num_pos    = posCount;
    obj.num_f_rows = frowCount;
    obj.num_verts  = vertCount; // Assume triangulation

    obj.pos   = malloc(sizeof(float) * posCount * 3);
    obj.norms = malloc(sizeof(float) * normalCount * 3);
    obj.texs  = malloc(sizeof(float) * texCount * 2);

    obj.indices = malloc(sizeof(vertindices_t) * vertCount);

    printf("poss : %d\n", posCount);
    printf("norms: %d\n", normalCount);
    printf("texs : %d\n", texCount);
    printf("faces: %d\n", frowCount);
    printf("verts: %d\n", frowCount * 3);

    // Reset file pointer to the start of the file
    fseek(fp, 0, SEEK_SET);

    // bounding box
    BoundingBox_t bbox;
    bbox.min[0] = bbox.min[1] = bbox.min[2] = 1e9;  // Initialize with large values
    bbox.max[0] = bbox.max[1] = bbox.max[2] = -1e9; // Initialize with small values

    int nrm_idx = 0, pos_idx = 0, tex_idx = 0, indi_idx = 0;
    while (fgets(linebuffer, sizeof(linebuffer), fp) != NULL)
    {
        if (linebuffer[0] == '#')
            continue;

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
            float pos[3] = {0};
            sscanf(linebuffer, "v %f %f %f", &pos[0], &pos[1], &pos[2]);

            bbox.min[0] = (pos[0] < bbox.min[0]) ? pos[0] : bbox.min[0];
            bbox.min[1] = (pos[1] < bbox.min[1]) ? pos[1] : bbox.min[1];
            bbox.min[2] = (pos[2] < bbox.min[2]) ? pos[2] : bbox.min[2];

            bbox.max[0] = (pos[0] > bbox.max[0]) ? pos[0] : bbox.max[0];
            bbox.max[1] = (pos[1] > bbox.max[1]) ? pos[1] : bbox.max[1];
            bbox.max[2] = (pos[2] > bbox.max[2]) ? pos[2] : bbox.max[2];

            obj.pos[pos_idx + 0] = pos[0];
            obj.pos[pos_idx + 1] = pos[1];
            obj.pos[pos_idx + 2] = pos[2];

            pos_idx += 3;
        }
        else if (linebuffer[0] == 'f')
        {
            int fv[4] = {0}, fvt[4] = {0}, fvn[4] = {0};

            int space_counter = 0;

            // Count the number of spaces in the line
            for (int i = 0; i < strlen(linebuffer); i++)
            {
                if (linebuffer[i] == ' ')
                    space_counter++;
            }

            // TODO : Read up to a space and save values
            if (space_counter >= 4)
            {
                assert(space_counter == 4);

                vertindices_t face_data[6] = {0};
                parse_f_line(linebuffer, space_counter, face_data);

                vertindices_t *fptr = &obj.indices[indi_idx + 0];

                fptr[0] = face_data[0];
                fptr[1] = face_data[1];
                fptr[2] = face_data[2];

                fptr[3] = face_data[0];
                fptr[4] = face_data[1];
                fptr[5] = face_data[3];

                indi_idx += 6;
            }
            else if (space_counter == 3)
            {
                int result = sscanf(linebuffer, "f %d/%d/%d %d/%d/%d %d/%d/%d",
                                    &fv[0], &fvt[0], &fvn[0],
                                    &fv[1], &fvt[1], &fvn[1],
                                    &fv[2], &fvt[2], &fvn[2]);
                assert(result == 9);
            }

            // NOTE : We could remove this copy here
            for (size_t f = 0; f < 3; f++) // Load in the first 3 values: f 0/0/0 1/1/1 2/2/2
            {
                vertindices_t *fptr = &obj.indices[indi_idx + f];
                fptr->v_idx         = fv[f] - 1;
                fptr->vt_idx        = fvt[f] - 1;
                fptr->vn_idx        = fvn[f] - 1;
            }
            indi_idx += 3;
        }
        if (linebuffer[0] == 'm' && linebuffer[1] == 't') // mtllib
        {
            // Parse material file
            _material_file(linebuffer);
        }
    }

    obj.bbox = bbox;

#if 0 // DEBUG
    printf("poss : %d\n", posCount);
    printf("texs : %d\n", texCount);
    printf("norms: %d\n", normalCount);
    printf("faces: %d\n", frowCount);
    printf("verts: %d\n", frowCount * 3);

    printf("Bounding box : min(%f, %f, %f), max(%f, %f, %f)\n",
           bbox.min[0], bbox.min[1], bbox.min[2],
           bbox.max[0], bbox.max[1], bbox.max[2]);

    for (size_t i = 0; i < obj.num_pos; i++)
    {
        const int index = i * 3;
        printf("v %f, %f, %f\n",
               obj.pos[index + 0],
               obj.pos[index + 1],
               obj.pos[index + 2]);
    }

    for (size_t i = 0; i < obj.num_norms; i++)
    {
        const int index = i * 3;
        printf("vn %f, %f, %f\n",
               obj.norms[index + 0],
               obj.norms[index + 1],
               obj.norms[index + 2]);
    }

    for (size_t i = 0; i < obj.num_texs; i++)
    {
        const int index = i * 2;
        printf("vt %f, %f\n",
               obj.texs[index + 0],
               obj.texs[index + 1]);
    }

    for (size_t i = 0; i < obj.num_verts; i++)
    {
        printf("%d, %d, %d\n",
               obj.indices[i].v_idx,
               obj.indices[i].vt_idx,
               obj.indices[i].vn_idx);
    }
#endif

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