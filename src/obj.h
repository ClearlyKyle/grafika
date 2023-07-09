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
    char *name;
    char *map_Kd;
} material_t;

typedef struct vertindices
{
    int v_idx, vn_idx, vt_idx;
} vertindices_t;

typedef struct shape
{
    int mat_idx, idx_offset, idx_end; // idx_offset is the offset into the 'indices' array
} shape_t;

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

    int         num_of_mats;
    material_t *mats;

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

static void _material_file(char *line, material_t **mat_data, int *num_of_mats)
{
    char *material_file_path = _read_string_until(line, ' ');

    removeNewline(++material_file_path); // NOTE : this might be unsafe

    printf("material file path : |%s|\n", material_file_path);

    FILE *material_fp = NULL;
    fopen_s(&material_fp, material_file_path, "r");
    assert(material_fp);

    material_t *curr_mat    = NULL;
    int         mat_counter = 0;

    char linebuffer[256]; // Adjust the size
    while (fgets(linebuffer, sizeof(linebuffer), material_fp) != NULL)
    {
        // Check for a new material
        if (strncmp(linebuffer, "newmtl", 6) == 0)
        {
            // Create new material struct
            mat_counter++;
            printf("Creating new mat\n");
            *mat_data = realloc(*mat_data, sizeof(material_t) * mat_counter);

            curr_mat = mat_data[mat_counter - 1];

            // Get material name
            char *material_name = _read_string_until(linebuffer, ' ');
            removeNewline(++material_name); // NOTE : this might be unsafe

            const size_t str_len = strlen(material_name);
            curr_mat->name       = malloc(sizeof(char) * str_len);

            strncpy_s(curr_mat->name, str_len, material_name, str_len);
        }
        if (strncmp(linebuffer, "map_Kd", 6) == 0)
        {
            char *diffuse_path = _read_string_until(linebuffer, ' ');
            removeNewline(++diffuse_path); // NOTE : this might be unsafe

            const size_t str_len = strlen(diffuse_path) + 1;
            curr_mat->map_Kd     = calloc(str_len, sizeof(char));

            strncpy_s(curr_mat->map_Kd, str_len, diffuse_path, str_len);
        }
    }

    *num_of_mats = mat_counter;
    fclose(material_fp);
}

static void parse_f_line(char *fline, const int num_vertex_values, vertindices_t *index_data)
{
    (void)num_vertex_values;

    int   vertex_index, texture_index, normal_index;
    char *current_char = fline + 2; // Skip the leading 'f' character and whitespace

    int collected_indices = 0;

    while (*current_char != '\0')
    {
        vertex_index = atoi(current_char); // Extract vertex index

        // Find the position of the next /
        while (*current_char != '/')
        {
            current_char++;
            if (*current_char == '\0' || *current_char == '\n')
                return; // Reached the end of the line
        }
        current_char++; // Move to the next character after the space

        // Extract texture index
        texture_index = atoi(current_char);

        // Find the position of the next /
        while (*current_char != '/')
        {
            current_char++;
            if (*current_char == '\0' || *current_char == '\n')
                return; // Reached the end of the line
        }
        current_char++; // Move to the next character after the space

        // Extract normal index
        normal_index = atoi(current_char);

        // Print the extracted values
        index_data[collected_indices].v_idx  = vertex_index - 1;
        index_data[collected_indices].vt_idx = texture_index - 1;
        index_data[collected_indices].vn_idx = normal_index - 1;
        collected_indices++;

        // Find the position of the next space or null character
        while (*current_char != ' ')
        {
            current_char++;
            if (*current_char == '\0' || *current_char == '\n')
                return; // Reached the end of the line
        }
        current_char++; // Move to the next character after the space
    }
    // printf("%s", fline);
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
            for (int i = 0; i < (int)strlen(linebuffer); i++)
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
            sscanf_s(linebuffer, "vn %f %f %f",
                     &obj.norms[nrm_idx + 0],
                     &obj.norms[nrm_idx + 1],
                     &obj.norms[nrm_idx + 2]);
            nrm_idx += 3;
        }
        else if (linebuffer[0] == 'v' && linebuffer[1] == 't')
        {
            sscanf_s(linebuffer, "vt %f %f",
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
            int space_counter = 0;

            // Count the number of spaces in the line
            for (int i = 0; i < (int)strlen(linebuffer); i++)
            {
                if (linebuffer[i] == ' ')
                    space_counter++;
            }

            // TODO : Read up to a space and save values
            if (space_counter >= 4)
            {
                assert(space_counter == 4);

                vertindices_t face_data[4] = {0};
                parse_f_line(linebuffer, space_counter, face_data);

                vertindices_t *fptr = &obj.indices[indi_idx];

                fptr[0] = face_data[0];
                fptr[1] = face_data[1];
                fptr[2] = face_data[2];

                fptr[3] = face_data[0];
                fptr[4] = face_data[2];
                fptr[5] = face_data[3];

                indi_idx += 6;
            }
            else if (space_counter == 3)
            {
                assert(space_counter == 3);

                vertindices_t face_data[3] = {0};
                parse_f_line(linebuffer, space_counter, face_data);

                vertindices_t *fptr = &obj.indices[indi_idx];

                fptr[0] = face_data[0];
                fptr[1] = face_data[1];
                fptr[2] = face_data[2];

                indi_idx += 3;
            }
        }
        if (linebuffer[0] == 'm' && linebuffer[1] == 't') // mtllib
        {
            // Parse material file
            _material_file(linebuffer, &obj.mats, &obj.num_of_mats);
        }
    }

    obj.bbox = bbox;

    printf("poss : %d\n", posCount);
    printf("norms: %d\n", normalCount);
    printf("texs : %d\n", texCount);
    printf("faces: %d\n", frowCount);
    printf("verts: %d\n", frowCount * 3);

    printf("num mats: %d\n", obj.num_of_mats);
    // for (int i = 0; i < obj.num_of_mats; i++)
    //{
    //     printf("%s - diffuse %s\n", obj.mats[i].name, obj.mats[i].map_Kd);
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
    if (obj->mats)
        free(obj->mats);

    *obj = (obj_t){0};

    printf("obj destroyed!\n");
}

#endif // __OBJ_H__