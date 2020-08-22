#include <algorithm>
#include <map>


const int MAX_COORD = 1 << 16;

struct Cell;

struct Cell_Point
{
    int index;
    int e[3];
};

struct Cell_Face
{
    int normal_axis;
    bool normal_flip;
    bool split;
    Cell_Point* points[4];
    Cell_Face* children[4];
    int index;
    Cell* owner;
    Cell* neighbor;
};

struct Cell
{
    bool split;
    Cell* children[8]; // ---, +--, -+-, ++-, --+, +-+, -++, +++
    Cell_Face* faces[6];  // -x, +x, -y, +y, -z, +z
    int index;
};

#define Alloc(Type) ((Type*) calloc(sizeof(Type), 1))


Cell_Point* make_point(int u, int v, int w)
{
    static std::map<std::pair<int, std::pair<int, int>>, Cell_Point*> m;
    Cell_Point* r = m[{ u, { v, w } }];
    if (r) return r;

    r = Alloc(Cell_Point);
    r->index = -1;
    r->e[0] = u;
    r->e[1] = v;
    r->e[2] = w;

    m[{ u, { v, w } }] = r;
    return r;
}

Cell_Point* midpoint(Cell_Point* a, Cell_Point* b)
{
    return make_point((a->e[0] + b->e[0]) / 2,
                      (a->e[1] + b->e[1]) / 2,
                      (a->e[2] + b->e[2]) / 2);
}

Cell_Face* make_face(int normal_axis, Cell_Point* a, Cell_Point* b, Cell_Point* c, Cell_Point* d)
{
    Cell_Face* face = Alloc(Cell_Face);
    face->index = -1;
    face->normal_axis = normal_axis;
    face->points[0] = a;
    face->points[1] = b;
    face->points[2] = c;
    face->points[3] = d;
    return face;
}

Cell_Face* midface(Cell_Face* a, Cell_Face* b)
{
    assert(a->normal_axis == b->normal_axis);
    Cell_Point* p0 = midpoint(a->points[0], b->points[0]);
    Cell_Point* p1 = midpoint(a->points[1], b->points[1]);
    Cell_Point* p2 = midpoint(a->points[2], b->points[2]);
    Cell_Point* p3 = midpoint(a->points[3], b->points[3]);
    return make_face(a->normal_axis, p0, p1, p2, p3);
}

void split_face(Cell_Face* face)
{
    if (face->split) return;
    face->split = true;

    Cell_Point* p0  = face->points[0];
    Cell_Point* p1  = face->points[1];
    Cell_Point* p2  = face->points[2];
    Cell_Point* p3  = face->points[3];
    Cell_Point* p01 = midpoint(p0, p1);
    Cell_Point* p12 = midpoint(p1, p2);
    Cell_Point* p23 = midpoint(p2, p3);
    Cell_Point* p30 = midpoint(p3, p0);
    Cell_Point* pin = midpoint(p01, p23);

    face->split = true;
    face->children[0] = make_face(face->normal_axis, p0, p01, pin, p30);
    face->children[1] = make_face(face->normal_axis, p01, p1, p12, pin);
    face->children[2] = make_face(face->normal_axis, p30, pin, p23, p3);
    face->children[3] = make_face(face->normal_axis, pin, p12, p2, p23);
}

Cell* make_cell(Cell_Face* x0, Cell_Face* x1, Cell_Face* y0, Cell_Face* y1, Cell_Face* z0, Cell_Face* z1)
{
    Cell* cell = Alloc(Cell);
    cell->index = -1;
    cell->faces[0] = x0;
    cell->faces[1] = x1;
    cell->faces[2] = y0;
    cell->faces[3] = y1;
    cell->faces[4] = z0;
    cell->faces[5] = z1;
    return cell;
}

void split_cell(Cell* cell)
{
    if (cell->split) return;
    cell->split = true;

    for (int i = 0; i < 6; i++)
        split_face(cell->faces[i]);

    Cell_Face* fx[3];
    fx[0] = cell->faces[0];
    fx[1] = midface(cell->faces[0], cell->faces[1]);
    fx[2] = cell->faces[1];
    split_face(fx[1]);

    Cell_Face* fy[3];
    fy[0] = cell->faces[2];
    fy[1] = midface(cell->faces[2], cell->faces[3]);
    fy[2] = cell->faces[3];
    split_face(fy[1]);

    Cell_Face* fz[3];
    fz[0] = cell->faces[4];
    fz[1] = midface(cell->faces[4], cell->faces[5]);
    fz[2] = cell->faces[5];
    split_face(fz[1]);

    for (int z = 0; z < 2; z++)
        for (int y = 0; y < 2; y++)
            for (int x = 0; x < 2; x++)
            {
                int xy = y * 2 + x;
                int xz = z * 2 + x;
                int yz = z * 2 + y;

                Cell_Face* x0 = fx[x + 0]->children[yz];
                Cell_Face* x1 = fx[x + 1]->children[yz];
                Cell_Face* y0 = fy[y + 0]->children[xz];
                Cell_Face* y1 = fy[y + 1]->children[xz];
                Cell_Face* z0 = fz[z + 0]->children[xy];
                Cell_Face* z1 = fz[z + 1]->children[xy];

                int i = z * 4 + y * 2 + x;
                cell->children[i] = make_cell(x0, x1, y0, y1, z0, z1);
            }

}

Cell_Point* face_center(Cell_Face* face)
{
    Cell_Point* p0 = midpoint(face->points[0], face->points[1]);
    Cell_Point* p1 = midpoint(face->points[2], face->points[3]);
    return midpoint(p0, p1);
}

Cell_Point* cell_center(Cell* cell)
{
    Cell_Point* pl = face_center(cell->faces[0]);
    Cell_Point* ph = face_center(cell->faces[1]);
    return midpoint(pl, ph);
}


void assign_cell_to_face(Cell* cell, Cell_Face* face, int* next_face_index)
{
    if (face->split)
    {
        for (int i = 0; i < 4; i++)
            assign_cell_to_face(cell, face->children[i], next_face_index);
        return;
    }

    assert(face->index == -1);
    assert(!face->neighbor);
    if (face->owner)
    {
        assert(cell->index > face->owner->index);
        face->neighbor = cell;
        face->index = (*next_face_index)++;
    }
    else
    {
        face->owner = cell;
    }
}

void assign_indices(Cell* cell, int* next_cell_index, int* next_face_index)
{
    if (!cell) return;
    assert(cell->index == -1);

    if (cell->split)
    {
        for (int i = 0; i < 8; i++)
            assign_indices(cell->children[i], next_cell_index, next_face_index);
        return;
    }

    cell->index = (*next_cell_index)++;
    for (int i = 0; i < 6; i++)
        assign_cell_to_face(cell, cell->faces[i], next_face_index);
}


void assign_face_boundary_index(Cell* cell, Cell_Face* face, int* next_face_index, bool normal_flip)
{
    if (face->split)
    {
        for (int i = 0; i < 4; i++)
            assign_face_boundary_index(cell, face->children[i], next_face_index, normal_flip);
        return;
    }

    if (face->neighbor)
        return;

    assert(face->owner == cell);
    assert(face->index == -1);
    face->index = (*next_face_index)++;
    face->normal_flip = normal_flip;
}

void assign_boundary_indices(Cell* cell, int* next_face_index)
{
    if (!cell) return;
    if (cell->split)
    {
        for (int i = 0; i < 8; i++)
            assign_boundary_indices(cell->children[i], next_face_index);
        return;
    }

    for (int i = 0; i < 6; i++)
    {
        bool normal_flip = (i % 2) == 0;
        assign_face_boundary_index(cell, cell->faces[i], next_face_index, normal_flip);
    }
}


void assign_face_point_indices(Cell_Face* face, Cell_Face** faces, int* next_point_index)
{
    if (face->split)
    {
        for (int i = 0; i < 4; i++)
            assign_face_point_indices(face->children[i], faces, next_point_index);
        return;
    }

    faces[face->index] = face;
    for (int i = 0; i < 4; i++)
    {
        Cell_Point* p = face->points[i];
        if (p->index == -1)
            p->index = (*next_point_index)++;
    }
}

void assign_point_indices(Cell* cell, Cell_Face** faces, int* next_point_index)
{
    if (!cell) return;
    if (cell->split)
    {
        for (int i = 0; i < 8; i++)
            assign_point_indices(cell->children[i], faces, next_point_index);
        return;
    }

    for (int i = 0; i < 6; i++)
        assign_face_point_indices(cell->faces[i], faces, next_point_index);
}


void assign_face_points(Cell_Face* face, Cell_Point** points)
{
    if (face->split)
    {
        for (int i = 0; i < 4; i++)
            assign_face_points(face->children[i], points);
        return;
    }

    for (int i = 0; i < 4; i++)
    {
        Cell_Point* p = face->points[i];
        points[p->index] = p;
    }
}

void assign_points(Cell* cell, Cell_Point** points)
{
    if (!cell) return;
    if (cell->split)
    {
        for (int i = 0; i < 8; i++)
            assign_points(cell->children[i], points);
        return;
    }

    for (int i = 0; i < 6; i++)
        assign_face_points(cell->faces[i], points);
}



struct Cell_Tree
{
    Cell* root;
    Point(*transform)(double u, double v, double w);
};


static FILE* open_foam_file(char* dir, char* object, char* clazz, bool is_ascii)
{
    char path[256];
    sprintf(path, "%s/%s", dir, object);

    FILE* out = fopen(path, "wb");
    if (!out)
    {
        printf("ERROR: can't open %s\n", path);
        return NULL;
    }

    printf(" .. %s\n", path);

    fprintf(out, "/*--------------------------------*- C++ -*----------------------------------*\\\n");
    fprintf(out, "| =========                 |                                                 |\n");
    fprintf(out, "| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |\n");
    fprintf(out, "|  \\\\    /   O peration     | Version:  1906                                  |\n");
    fprintf(out, "|   \\\\  /    A nd           | Web:      www.OpenFOAM.com                      |\n");
    fprintf(out, "|    \\\\/     M anipulation  |                                                 |\n");
    fprintf(out, "\\*---------------------------------------------------------------------------*/\n");
    fprintf(out, "FoamFile\n");
    fprintf(out, "{\n");
    fprintf(out, "    version     2.0;\n");
    fprintf(out, "    format      %s;\n", is_ascii ? "ascii" : "binary");
    fprintf(out, "    class       %s;\n", clazz);
    if (!is_ascii)
        fprintf(out, "    arch        \"LSB;label=32;scalar=64\";\n");
    fprintf(out, "    location    \"constant/polyMesh\";\n");
    fprintf(out, "    object      %s;\n", object);
    fprintf(out, "}\n");
    fprintf(out, "\n");

    return out;
}


void recursive_split(Cell* cell, int s, int ox, int oy, int oz, int cx, int cy, int cz)
{
    if (!s) return;

    split_cell(cell);
    for (int z = 0; z < 2; z++)
        for (int y = 0; y < 2; y++)
            for (int x = 0; x < 2; x++)
            {
                int i = z * 4 + y * 2 + x;
                int cox = ox + x * (1 << (s - 1));
                int coy = oy + y * (1 << (s - 1));
                int coz = oz + z * (1 << (s - 1));
                if (cox >= cx || coy >= cy || coz >= cz)
                    cell->children[i] = NULL;
                else
                    recursive_split(cell->children[i], s - 1, cox, coy, coz, cx, cy, cz);
            }
}

Cell* make_unit_cell(int cx, int cy, int cz)
{
    int s = cx;
    if (cy > s) s = cy;
    if (cz > s) s = cz;

    Cell_Point* p[8] =
    {
        make_point(0         << (s - cx), 0         << (s - cy), 0         << (s - cz)),
        make_point(MAX_COORD << (s - cx), 0         << (s - cy), 0         << (s - cz)),
        make_point(MAX_COORD << (s - cx), MAX_COORD << (s - cy), 0         << (s - cz)),
        make_point(0         << (s - cx), MAX_COORD << (s - cy), 0         << (s - cz)),
        make_point(0         << (s - cx), 0         << (s - cy), MAX_COORD << (s - cz)),
        make_point(MAX_COORD << (s - cx), 0         << (s - cy), MAX_COORD << (s - cz)),
        make_point(MAX_COORD << (s - cx), MAX_COORD << (s - cy), MAX_COORD << (s - cz)),
        make_point(0         << (s - cx), MAX_COORD << (s - cy), MAX_COORD << (s - cz)),
    };

    Cell_Face* x0 = make_face(0, p[0], p[3], p[7], p[4]);
    Cell_Face* x1 = make_face(0, p[1], p[2], p[6], p[5]);
    Cell_Face* y0 = make_face(1, p[0], p[1], p[5], p[4]);
    Cell_Face* y1 = make_face(1, p[3], p[2], p[6], p[7]);
    Cell_Face* z0 = make_face(2, p[0], p[1], p[2], p[3]);
    Cell_Face* z1 = make_face(2, p[4], p[5], p[6], p[7]);

    Cell* root = make_cell(x0, x1, y0, y1, z0, z1);
    recursive_split(root, s, 0, 0, 0, 1 << cx, 1 << cy, 1 << cz);

    return root;
}


Point transform(Cell_Tree* tree, Cell_Point* p)
{
    double u = p->e[0] / (double) MAX_COORD;
    double v = p->e[1] / (double) MAX_COORD;
    double w = p->e[2] / (double) MAX_COORD;
    return tree->transform(u, v, w);
}

void refine(Cell_Tree* tree, Cell** cell, int limit, bool(*inside)(Point p))
{
    Cell* c = *cell;
    if (!c) return;
    if (limit <= 0) return;

    if (!c->split)
    {
        int in = 0;
        if (inside(transform(tree, c->faces[0]->points[0]))) in++;
        if (inside(transform(tree, c->faces[0]->points[1]))) in++;
        if (inside(transform(tree, c->faces[0]->points[2]))) in++;
        if (inside(transform(tree, c->faces[0]->points[3]))) in++;
        if (inside(transform(tree, c->faces[1]->points[0]))) in++;
        if (inside(transform(tree, c->faces[1]->points[1]))) in++;
        if (inside(transform(tree, c->faces[1]->points[2]))) in++;
        if (inside(transform(tree, c->faces[1]->points[3]))) in++;
        if (!in) return;

        split_cell(c);
    }

    for (int i = 0; i < 8; i++)
        refine(tree, &c->children[i], limit - 1, inside);
}

void carve_out(Cell_Tree* tree, Cell** cell, int limit, bool(*inside)(Point p))
{
    Cell* c = *cell;
    if (!c) return;

    if (!c->split)
    {
        int in = 0;
        if (inside(transform(tree, c->faces[0]->points[0]))) in++;
        if (inside(transform(tree, c->faces[0]->points[1]))) in++;
        if (inside(transform(tree, c->faces[0]->points[2]))) in++;
        if (inside(transform(tree, c->faces[0]->points[3]))) in++;
        if (inside(transform(tree, c->faces[1]->points[0]))) in++;
        if (inside(transform(tree, c->faces[1]->points[1]))) in++;
        if (inside(transform(tree, c->faces[1]->points[2]))) in++;
        if (inside(transform(tree, c->faces[1]->points[3]))) in++;
        if (!in) return;

        if (limit > 0 && in < 8)
        {
            split_cell(c);
        }
        else
        {
            *cell = NULL;
            return;
        }
    }

    bool ok = false;
    for (int i = 0; i < 8; i++)
    {
        carve_out(tree, &c->children[i], limit - 1, inside);
        if (c->children[i])
            ok = true;
    }

    if (!ok)
    {
        c->split = false;
        *cell = NULL;
    }
}


void output_cell_tree(Cell_Tree* tree, char* dir, bool ascii)
{
    Cell* root = tree->root;

    int next_cell_index = 0;
    int next_face_index = 0;
    int next_point_index = 0;
    assign_indices(root, &next_cell_index, &next_face_index);

    int boundary_start = next_face_index;
    assign_boundary_indices(root, &next_face_index);
    int boundary_count = next_face_index - boundary_start;

    Cell_Face** faces = (Cell_Face**) malloc(sizeof(Cell_Face*) * next_face_index);
    assign_point_indices(root, faces, &next_point_index);

    Cell_Point** points = (Cell_Point**) malloc(sizeof(Cell_Face*) * next_point_index);
    assign_points(root, points);

    for (int i = 0; i < next_face_index; i++)
    {
        Cell_Face* face = faces[i];
        if (face->normal_axis == 1)
            face->normal_flip = !face->normal_flip;
    }

    FILE* out;
    printf("Writing the OpenFOAM structure:\n");

    if (out = open_foam_file(dir, "boundary", "polyBoundaryMesh", ascii))
    {
        fprintf(out, "1\n");
        fprintf(out, "(\n");
        fprintf(out, "    boundaries\n");
        fprintf(out, "    {\n");
        fprintf(out, "        type patch;\n");
        fprintf(out, "        nFaces %d;\n", boundary_count);
        fprintf(out, "        startFace %d;\n", boundary_start);
        fprintf(out, "    }\n");
        fprintf(out, ")\n");
        fclose(out);
    }

    if (out = open_foam_file(dir, "faces", ascii ? "faceList" : "faceCompactList", ascii))
    {
        if (ascii)
        {
            fprintf(out, "%d\n(\n", next_face_index);
        }
        else
        {
            fprintf(out, "%d\n(", next_face_index + 1);
            for (int i = 0; i <= next_face_index; i++)
            {
                uint32_t start = i * 4;
                fwrite(&start, 4, 1, out);
            }
            fprintf(out, ")\n%d\n(", next_face_index * 4);
        }
        for (int i = 0; i < next_face_index; i++)
        {
            uint32_t index[4];
            for (int j = 0; j < 4; j++)
                index[j] = faces[i]->points[faces[i]->normal_flip ? 3 - j : j]->index;

            if (ascii)
                fprintf(out, "4(%d %d %d %d)\n", index[0], index[1], index[2], index[3]);
            else
                fwrite(index, 4, 4, out);
        }
        fprintf(out, ")\n");
        fclose(out);
    }

    if (out = open_foam_file(dir, "neighbour", "labelList", ascii))
    {
        fprintf(out, "%d\n", next_face_index);
        fprintf(out, ascii ? "(\n" : "(");
        for (int i = 0; i < next_face_index; i++)
        {
            uint32_t index = faces[i]->neighbor ? faces[i]->neighbor->index : -1;
            if (ascii)
                fprintf(out, "%d\n", index);
            else
                fwrite(&index, 4, 1, out);
        }
        fprintf(out, ")\n");
        fclose(out);
    }

    if (out = open_foam_file(dir, "owner", "labelList", ascii))
    {
        fprintf(out, "%d\n", next_face_index);
        fprintf(out, ascii ? "(\n" : "(");
        for (int i = 0; i < next_face_index; i++)
        {
            uint32_t index = faces[i]->owner->index;
            if (ascii)
                fprintf(out, "%d\n", index);
            else
                fwrite(&index, 4, 1, out);
        }
        fprintf(out, ")\n");
        fclose(out);
    }

    if (out = open_foam_file(dir, "points", "vectorField", ascii))
    {
        fprintf(out, "%d\n", next_point_index);
        fprintf(out, ascii ? "(\n" : "(");
        for (int i = 0; i < next_point_index; i++)
        {
            Point p = transform(tree, points[i]);
            if (ascii)
                fprintf(out, "(%f %f %f)\n", p.x, p.y, p.z);
            else
                fwrite(p.e, 8, 3, out);
        }
        fprintf(out, ")\n");
        fclose(out);
    }
}
