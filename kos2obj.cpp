#define _USE_MATH_DEFINES
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#if defined(_WIN32)
#include <direct.h>
#elif defined(__linux__)
#include <sys/stat.h>
#endif


const double DEG = 3.14159265359 / 180.0;

#include "geometry.inl"
#include "model.inl"
#include "triangle_mesh.inl"
#include "cell_tree.inl"


void contours_to_mesh(Mesh* mesh, Contour* c, int n)
{
    int stride = c->n_points;
    int first_base = -1;
    int last_base = -1;
    for (int i = 0; i < 2 * n; i++)
    {
        double sx = 1.0;
        double sy = 1.0;
        double sz = 1.0;

        int j = n - i - 1;
        if (i >= n)
        {
            j = i - n;
            sy *= -1;
        }

        int base = reserve_points(mesh, c[j].n_points);
        for (int i = 0; i < c[j].n_points; i++)
        {
            Point p = c[j].points[i];
            mesh->points[base + i] = { p.x * sx, p.y * sy, p.z * sz };
        }

        if (first_base == -1)
        {
            first_base = base;
        }
        else
        {
            int* indices = reserve_triangles(mesh, (stride - 1) * 2);
            for (int j = 0; j < stride - 1; j++)
            {
                *(indices++) = last_base + j;
                *(indices++) = last_base + j + 1;
                *(indices++) = base + j;

                *(indices++) = last_base + j + 1;
                *(indices++) = base + j + 1;
                *(indices++) = base + j;
            }
        }

        last_base = base;
    }

    int* indices = reserve_triangles(mesh, (n - 1) * 4);
    for (int i = 0; i < n - 1; i++)
    {
        int a1 = first_base + stride * (i);
        int b1 = first_base + stride * (i + 1);
        int c1 = first_base + stride * (2 * n - i - 2);
        int d1 = first_base + stride * (2 * n - i - 1);

        *(indices++) = b1;
        *(indices++) = d1;
        *(indices++) = a1;

        *(indices++) = b1;
        *(indices++) = c1;
        *(indices++) = d1;

        int a2 = a1 + stride - 1;
        int b2 = b1 + stride - 1;
        int c2 = c1 + stride - 1;
        int d2 = d1 + stride - 1;

        *(indices++) = b2;
        *(indices++) = a2;
        *(indices++) = d2;

        *(indices++) = b2;
        *(indices++) = d2;
        *(indices++) = c2;
    }
}

void generate_meshes(Model* model, Plane cut_planes[], int n_cut_planes, Mesh* out_hull, Mesh out_skeg[])
{
    printf("Generating mesh for model %s\n", model->name);
    printf(" .. hull\n");

    // model contours
    int n_contours = model->n_cross + 2;
    Contour* contours = (Contour*) calloc(sizeof(Contour), n_contours);
    contours[0] = sub_contour(&model->back, model->back.i_sharp[0], model->back.n_points);
    for (int i = 0; i < model->n_cross; i++)
        contours[i + 1] = model->cross[i];
    contours[n_contours - 1].n_points = 2;
    contours[n_contours - 1].points = (Point*) calloc(sizeof(Point), 2);
    contours[n_contours - 1].points[0] = model->front.points[model->front.n_points - 1];
    contours[n_contours - 1].points[1] = model->front.points[model->front.n_points - 1];

    for (int i = 0; i < n_contours; i++)
        contours[i] = interpolate_contour(&contours[i]);

    // hull
    Contour ship_contours[INTERPOLATED];
    transpose_interpolated(contours, n_contours, ship_contours);
    for (int i = 0; i < INTERPOLATED; i++)
    {
        Contour c = interpolate_contour(&ship_contours[i]);
        ship_contours[i] = c;
    }

    Mesh hull_mesh = { 0 };
    contours_to_mesh(&hull_mesh, ship_contours, INTERPOLATED);
    *out_hull = hull_mesh;

    // skeg
    printf(" .. skeg\n");
    Contour skeg_contours[INTERPOLATED];
    transpose_interpolated(ship_contours, INTERPOLATED, skeg_contours);
    for (int i = 0; i < INTERPOLATED; i++)
    {
        double skeg_thickness = 0.01;
        skeg_contours[i] = skeg_contour(&skeg_contours[i], 0.5 * skeg_thickness);
    }

    Mesh skeg_mesh = { 0 };
    contours_to_mesh(&skeg_mesh, skeg_contours, INTERPOLATED);

    // cut skeg
    printf(" .. cutting skeg\n");
    for (int i = 0; i < n_cut_planes - 1; i++)
    {
        Mesh cut_mesh = copy_mesh(skeg_mesh);
        cut_plane(&cut_mesh, cut_planes[i], true);
        cut_plane(&cut_mesh, { cut_planes[i + 1].p, mul(cut_planes[i + 1].n, -1) }, true);
        out_skeg[i] = cut_mesh;
    }
}



char root_path[256];
char constant_path[256];
char trisurface_path[256];
char polymesh_path[256];

char kos_path[256];
char stl_path[16][256];

int ends_with(const char* str, const char* suffix)
{
    if (!str || !suffix) return 0;
    size_t lenstr    = strlen(str);
    size_t lensuffix = strlen(suffix);
    if (lensuffix > lenstr) return 0;
    return strncmp(str + lenstr - lensuffix, suffix, lensuffix) == 0;
}


int main(int argc, char** argv)
{
#define SKIP_INPUT 1

#if SKIP_INPUT
    sprintf(kos_path, "m808.kos");
#else
    do
    {
        printf("KOS path: ");
        scanf("%s", kos_path);
    }
    while (!ends_with(kos_path, ".kos") && !ends_with(kos_path, ".KOS"));
#endif

    FILE* kos = fopen(kos_path, "rt");
    if (!kos)
    {
        printf("Error: Can't open %s\n", kos_path);
        return EXIT_FAILURE;
    }


    sprintf(root_path,       "%.*s", (int) strlen(kos_path) - 4, kos_path);
    sprintf(constant_path,   "%s/constant", root_path);
    sprintf(trisurface_path, "%s/triSurface", constant_path);
    sprintf(polymesh_path,   "%s/polyMesh", constant_path);

#if defined(_WIN32)
    _mkdir(root_path);
    _mkdir(constant_path);
    _mkdir(trisurface_path);
    _mkdir(polymesh_path);
#elif defined(__linux__)
    mkdir(root_path, 0777);
    mkdir(constant_path, 0777);
    mkdir(trisurface_path, 0777);
    mkdir(polymesh_path, 0777);
#endif

    sprintf(stl_path[0], "%s/A.stl", trisurface_path);
    sprintf(stl_path[1], "%s/B.stl", trisurface_path);
    sprintf(stl_path[2], "%s/C.stl", trisurface_path);

    FILE* stlA = fopen(stl_path[0], "wb");
    FILE* stlB = fopen(stl_path[1], "wb");
    FILE* stlC = fopen(stl_path[2], "wb");
    if (!stlA || !stlB || !stlC)
    {
        printf("Error: Can't open %s or %s or %s\n", stl_path[0], stl_path[1], stl_path[2]);
        return EXIT_FAILURE;
    }


    double T, Fx, Fy, Fz;
#if SKIP_INPUT
    T = 0.2;
    Fx = 0.0;
    Fy = 0.0;
    Fz = 0.0;
#else
    printf("T: ");   scanf("%lf", &T);
    printf("FIx: "); scanf("%lf", &Fx);
    printf("FIy: "); scanf("%lf", &Fy);
    printf("FIz: "); scanf("%lf", &Fz);
#endif


    Model model;
    load_model(kos, &model);


    static Mesh meshes[4];

    Plane cut_planes[4];
    {
        cut_planes[0] = { { -1.131, 0, 0.066 }, { cos( 30 * DEG), 0, sin( 30 * DEG) } };
        cut_planes[1] = { { -0.792, 0, 0     }, { cos( 30 * DEG), 0, sin( 30 * DEG) } };
        cut_planes[2] = { { -0.419, 0, 0     }, { cos( 30 * DEG), 0, sin( 30 * DEG) } };
        cut_planes[3] = { { -0.061, 0, 0     }, { cos(-45 * DEG), 0, sin(-45 * DEG) } };
    }

    generate_meshes(&model, cut_planes, 4, &meshes[0], meshes + 1);

    printf("Applying the model transform\n");
    for (int i = 0; i < 4; i++)
    {
        rotate   (&meshes[i], -Fy * DEG, 1);
        rotate   (&meshes[i],  Fx * DEG, 0);
        rotate   (&meshes[i], -Fz * DEG, 2);
        translate(&meshes[i], 0, 0, -T);
        cut_plane(&meshes[i], { { 0, 0, 0 }, { 0, 0, -1 } }, false);
        rotate   (&meshes[i], 180 * DEG, 0);
    }

    //
    // STL output
    //

    printf("Writing to STL\n");

    add_mesh(&meshes[0], meshes[3]);
    output_mesh_to_stl(stlA, &meshes[0]);
    output_mesh_to_stl(stlB, &meshes[2]);
    output_mesh_to_stl(stlC, &meshes[1]);

    //
    // OpenFOAM mesh output
    //

    for (int i = 0; i < 4; i++)
        rotate(&meshes[i], 180 * DEG, 0);

    printf("Computing the Z buffer\n");

    static const int X_SIZE = 2048;
    static const int Y_SIZE = 2048;
    static double z_buffer_data[Y_SIZE + 1][X_SIZE + 1];
    static Z_Buffer z_buffer;
    z_buffer.size_x = X_SIZE;
    z_buffer.size_y = Y_SIZE;
    z_buffer.stride = Y_SIZE + 1;
    z_buffer.data = &z_buffer_data[0][0];

    set_z_buffer_bounds(&z_buffer, &meshes[0]);
    for (int i = 0; i < 3; i++)
        mesh_render_z_buffer(&z_buffer, &meshes[i]);

    printf("Carving an OpenFOAM structure\n");

    Cell_Tree tree = {};
    tree.transform = [](double u, double v, double w) -> Point
    {
        double r =  1 + 3 * v;
        double x = -r * cos(u * M_PI);
        double y =  r * sin(u * M_PI) - 2.5;
        double z = -2 + 2 * w;
        return { x, y, z };
    };

    tree.root = make_unit_cell(9, 7, 7);
    carve_out(&tree, &tree.root, 12, [](Point p) -> bool
    // tree.root = make_unit_cell(7, 5, 5);
    // carve_out(&tree, &tree.root, 10, [](Point p) -> bool
    {
        if (!in_bounds(&z_buffer, p.x, p.y)) return false;

        int x, y;
        to_pixel(&z_buffer, p.x, p.y, &x, &y);
        double z = z_buffer.data[y * z_buffer.stride + x];
        return z < p.z;
    });

    /*tree.root = make_unit_cell(1, 0, 0);
    split_cell(tree.root->children[0]);
    for (int z = 0; z < 2; z++)
        for (int y = 0; y < 2; y++)
            for (int x = 0; x < 2; x++)
                if (!x || y || z)
                    tree.root->children[0]->children[z * 4 + y * 2 + x] = NULL;*/

    /*refine(&tree, &tree.root, 7, [](Point p) -> bool
    {
        Point d = sub(p, { 0, 1.5, 0 });
        return len(d) < (0.5 * 0.5);
    });*/

    bool ascii = true;
    output_cell_tree(&tree, polymesh_path, ascii);


    fclose(kos);
    fclose(stlA);
    fclose(stlB);
    fclose(stlC);
    printf("Done\n");
    return EXIT_SUCCESS;
}
