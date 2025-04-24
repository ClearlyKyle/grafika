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
    char *map_bump;
    char *disp;
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

    size_t      num_of_mats;
    material_t *mats;

    float *pos;
    float *norms;
    float *texs;

    vertindices_t *indices; // Index values for 1 vertex
} obj_t;

/*
dest    - needs to be big enough to store src + null
destsz  - strlen(src) + 1;
src     - source string to copy
srcsz   - how much of the src string to copy, no including the nul
*/
static void _str_copy(char *dest, const size_t destsz, const char *src, const size_t srcsz)
{
    if (destsz >= (srcsz + 1))
    {
        memmove(dest, src, srcsz);
        dest[srcsz] = '\0';
    }
    else
    {
        fprintf(stderr, "destination for string |%s| is too small\n", src);
        memmove(dest, src, destsz);
        dest[destsz - 1] = '\0';
    }
}

static char *_read_string_until(char *inputString, const char stopChar)
{
    // TODO : this should take in some kind of length
    ASSERT(inputString, "string should not be NULL");

    while (*inputString != '\0' && *inputString != stopChar)
        inputString++;

    return inputString;
}

static void _remove_new_line(char *str)
{
    ASSERT(str, "string should not be NULL");

    const size_t length = strcspn(str, "\n"); // TODO : remove this function call

    // TODO : refactor this
    if (str[length] == '\n' || str[length] == '\r')
        str[length] = '\0'; // Replace newline with null-termination

    if (str[length - 1] == '\n' || str[length - 1] == '\r')
        str[length - 1] = '\0'; // Replace newline with null-termination
}

static void _material_file(char *line, material_t **mat_data, size_t *num_of_mats)
{
    char *material_file_path = _read_string_until(line, ' ');

    _remove_new_line(++material_file_path); // NOTE : this might be unsafe

    FILE *material_fp = NULL;
    material_fp       = fopen(material_file_path, "r");
    ASSERT(material_fp, "Error opening matrial file: %s\n", material_file_path);

    material_t *curr_mat    = NULL;
    size_t      mat_counter = 0;

    char linebuffer[256]; // Adjust the size
    while (fgets(linebuffer, sizeof(linebuffer), material_fp) != NULL)
    {
        // Check for a new material
        if (strncmp(linebuffer, "newmtl", 6) == 0)
        {
            // Create new material struct
            mat_counter++;
            LOG("Creating new mat\n");
            *mat_data = realloc(*mat_data, sizeof(material_t) * mat_counter);

            curr_mat = mat_data[mat_counter - 1];
            memset(curr_mat, 0, sizeof(material_t));

            // Get material name
            char *material_name = _read_string_until(linebuffer, ' ');
            _remove_new_line(++material_name); // NOTE : this might be unsafe

            const size_t str_len  = strlen(material_name); // not including the null terminating character
            const size_t dst_size = str_len + 1;
            curr_mat->name        = malloc(sizeof(char) * dst_size);

            _str_copy(curr_mat->name, dst_size, material_name, str_len);
            // curr_mat->name[str_len - 1] = '\0';
        }
        if (strncmp(linebuffer, "map_Kd", 6) == 0) // Diffuse Map
        {
            char *diffuse_path = _read_string_until(linebuffer, ' ');
            _remove_new_line(++diffuse_path); // NOTE : this might be unsafe

            const size_t str_len  = strlen(diffuse_path);
            const size_t dst_size = str_len + 1;
            curr_mat->map_Kd      = calloc(dst_size, sizeof(char));

            _str_copy(curr_mat->map_Kd, dst_size, diffuse_path, str_len);
            // curr_mat->map_Kd[str_len - 1] = '\0';

            LOG("curr_mat->map_Kd = |%s|\n", curr_mat->map_Kd);
        }
        if (strncmp(linebuffer, "map_bump", 8) == 0 || strncmp(linebuffer, "map_Bump", 8) == 0) // Bump Map
        {
            char *bump_path = _read_string_until(linebuffer, ' ');
            _remove_new_line(++bump_path); // NOTE : this might be unsafe

            const size_t str_len  = strlen(bump_path);
            const size_t dst_size = str_len + 1;

            curr_mat->map_bump = calloc(dst_size, sizeof(char));

            _str_copy(curr_mat->map_bump, dst_size, bump_path, str_len);
            // curr_mat->map_bump[str_len - 1] = '\0';

            LOG("curr_mat->map_bump = |%s|\n", curr_mat->map_bump);
        }
        if (strncmp(linebuffer, "disp", 4) == 0) // Displacement Map
        {
            char *disp_path = _read_string_until(linebuffer, ' ');
            _remove_new_line(++disp_path); // NOTE : this might be unsafe

            const size_t str_len  = strlen(disp_path);
            const size_t dst_size = str_len + 1;
            curr_mat->disp        = calloc(dst_size, sizeof(char));

            _str_copy(curr_mat->disp, dst_size, disp_path, str_len);
            // curr_mat->disp[str_len - 1] = '\0';

            LOG("curr_mat->disp = |%s|\n", curr_mat->disp);
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
    // LOG("%s", fline);
}

obj_t obj_load(const char *filename)
{
    FILE *fp = NULL;
    fp       = fopen(filename, "r");
    ASSERT(fp, "Error with opening object file : '%s'\n", filename);

    char linebuffer[MAX_CHAR_LINE_BUFFER]; // Adjust the size

    size_t posCount    = 0;
    size_t texCount    = 0;
    size_t normalCount = 0;
    size_t frowCount   = 0;
    size_t vertCount   = 0;
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
            int space_counter = 0;

            // Count the number of spaces in the line
            for (int i = 0; i < (int)strlen(linebuffer); i++)
            {
                if (linebuffer[i] == ' ')
                    space_counter++;
            }

            if (space_counter >= 4)
            {
                ASSERT(space_counter == 4, "space_counter = %d", space_counter);

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
                ASSERT(space_counter == 3, "space_counter = %d", space_counter);

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

    LOG("poss : %zu\n", posCount);
    LOG("norms: %zu\n", normalCount);
    LOG("texs : %zu\n", texCount);
    LOG("faces: %zu\n", frowCount);
    LOG("verts: %zu\n", frowCount * 3);

    LOG("num mats: %zu\n", obj.num_of_mats);
    for (size_t i = 0; i < obj.num_of_mats; i++)
    {
        LOG("%s - diffuse %s\n", obj.mats[i].name, obj.mats[i].map_Kd);
        LOG("%s - normal  %s\n", obj.mats[i].name, obj.mats[i].map_bump);
        LOG("%s - disp    %s\n", obj.mats[i].name, obj.mats[i].disp);
    }

    fclose(fp);
    return obj;
}

void obj_destroy(struct obj *obj)
{
    if (!obj) return;

    if (obj->pos) free(obj->pos);
    if (obj->norms) free(obj->norms);
    if (obj->texs) free(obj->texs);
    if (obj->indices) free(obj->indices);

    if (obj->mats)
    {
        for (size_t i = 0; i < obj->num_of_mats; i++)
        {
            if (obj->mats[i].name) free(obj->mats[i].name);
            if (obj->mats[i].map_Kd) free(obj->mats[i].map_Kd);
            if (obj->mats[i].map_bump) free(obj->mats[i].map_bump);
            if (obj->mats[i].disp) free(obj->mats[i].disp);
        }

        free(obj->mats);
    }

    LOG("obj destroyed!\n");
}

#endif // __OBJ_H__