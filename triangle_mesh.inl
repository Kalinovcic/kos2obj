struct Mesh
{
    int n_points;
    int c_points;
    Point* points;

    int n_indices;
    int c_indices;
    int* indices;
};

Mesh copy_mesh(Mesh m)
{
    Mesh r = { 0 };
    r.n_points  = r.c_points  = m.n_points;
    r.n_indices = r.c_indices = m.n_indices;
    r.points  = (Point*) malloc(sizeof(Point) * r.c_points);
    r.indices = (int*)   malloc(sizeof(int)   * r.c_indices);
    memcpy(r.points,  m.points,  sizeof(Point) * r.n_points);
    memcpy(r.indices, m.indices, sizeof(int)   * r.n_indices);
    return r;
}

void remove_unreferenced_points(Mesh* mesh)
{
    int* sum = (int*) malloc(sizeof(int) * mesh->n_points);
    memset(sum, 0xFF, sizeof(int) * mesh->n_points);

    for (int i = 0; i < mesh->n_indices; i++)
        sum[mesh->indices[i]] = 0;

    int n = mesh->n_points;
    mesh->n_points = 0;
    for (int i = 0; i < n; i++)
    {
        if (sum[i] == 0)
            mesh->points[mesh->n_points++] = mesh->points[i];
        if (i)
            sum[i] += sum[i - 1];
    }

    for (int i = 0; i < mesh->n_indices; i++)
        mesh->indices[i] += sum[mesh->indices[i]];

    free(sum);
}

int reserve_points(Mesh* mesh, int count)
{
    mesh->n_points += count;
    while (mesh->n_points > mesh->c_points)
    {
        mesh->c_points = mesh->c_points ? (mesh->c_points) * 2 : 1024;
        mesh->points = (Point*) realloc(mesh->points, sizeof(Point) * mesh->c_points);
    }
    return mesh->n_points - count;
}

int* reserve_triangles(Mesh* mesh, int count)
{
    mesh->n_indices += count * 3;
    while (mesh->n_indices > mesh->c_indices)
    {
        mesh->c_indices = mesh->c_indices ? (mesh->c_indices) * 2 : 1024;
        mesh->indices = (int*) realloc(mesh->indices, sizeof(int) * mesh->c_indices);
    }
    return mesh->indices + mesh->n_indices - count * 3;
}

void remove_triangle(Mesh* mesh, int first_index)
{
    mesh->indices[first_index + 0] = mesh->indices[mesh->n_indices - 3];
    mesh->indices[first_index + 1] = mesh->indices[mesh->n_indices - 2];
    mesh->indices[first_index + 2] = mesh->indices[mesh->n_indices - 1];
    mesh->n_indices -= 3;
}

void add_mesh(Mesh* mesh, Mesh src)
{
    int base = reserve_points(mesh, src.n_points);
    memcpy(mesh->points + base, src.points, src.n_points * sizeof(Point));

    int* indices = reserve_triangles(mesh, src.n_indices / 3);
    for (int i = 0; i < src.n_indices; i++)
        indices[i] = src.indices[i] + base;
}

void translate(Mesh* mesh, double dx, double dy, double dz)
{
    for (int i = 0; i < mesh->n_points; i++)
    {
        mesh->points[i].x += dx;
        mesh->points[i].y += dy;
        mesh->points[i].z += dz;
    }
}

void rotate(Mesh* mesh, double angle, int axis)
{
    int x = (axis + 1) % 3;
    int y = (axis + 2) % 3;

    double c = cos(angle);
    double s = sin(angle);

    for (int i = 0; i < mesh->n_points; i++)
    {
        Point p = mesh->points[i];
        mesh->points[i].e[x] = c * p.e[x] - s * p.e[y];
        mesh->points[i].e[y] = s * p.e[x] + c * p.e[y];
    }
}


Mesh* comp_mesh;
Point comp_middle;

int cut_comparator(const void* pa, const void* pb)
{
    Point a = comp_middle;
    Point b = comp_mesh->points[*(int*) pa];
    Point c = comp_mesh->points[*(int*) pb];
    double ang_b = atan2(b.z - a.z, b.y - a.y);
    double ang_c = atan2(c.z - a.z, c.y - a.y);
    if (ang_c > ang_b) return  1;
    if (ang_c < ang_b) return -1;
    return 0;
}

void cut_plane(Mesh* mesh, Plane plane, bool fill)
{
    double epsilon = 1e-6;

    int added_base = mesh->n_points;
    for (int i = 0; i < mesh->n_indices; i += 3)
    {
        int index[3];
        index[0] = mesh->indices[i + 0];
        index[1] = mesh->indices[i + 1];
        index[2] = mesh->indices[i + 2];

        Point p[3];
        p[0] = mesh->points[index[0]];
        p[1] = mesh->points[index[1]];
        p[2] = mesh->points[index[2]];

        double d[3];
        d[0] = dot_plane(plane, p[0]) + epsilon;
        d[1] = dot_plane(plane, p[1]) + epsilon;
        d[2] = dot_plane(plane, p[2]) + epsilon;

        if (d[0] >= 0 && d[1] >= 0 && d[2] >= 0)
            continue;

        remove_triangle(mesh, i);
        i -= 3;

        if (d[0] < 0 && d[1] < 0 && d[2] < 0)
            continue;

        for (int j = 0; j < 3; j++)
        {
            bool tri  = (d[j] >= 0 && d[(j + 1) % 3] <  0 && d[(j + 2) % 3] <  0);
            bool quad = (d[j] <  0 && d[(j + 1) % 3] >= 0 && d[(j + 2) % 3] >= 0);
            if (tri || quad)
            {
                int base = reserve_points(mesh, 2);
                mesh->points[base + 0] = intersect_line_plane(plane, p[j], p[(j + 1) % 3]);
                mesh->points[base + 1] = intersect_line_plane(plane, p[j], p[(j + 2) % 3]);

                if (tri)
                {
                    int* indices = reserve_triangles(mesh, 1);
                    *(indices++) = index[j];
                    *(indices++) = base + 0;
                    *(indices++) = base + 1;
                }
                else
                {
                    int* indices = reserve_triangles(mesh, 2);
                    *(indices++) = index[(j + 1) % 3];
                    *(indices++) = index[(j + 2) % 3];
                    *(indices++) = base + 1;
                    *(indices++) = index[(j + 1) % 3];
                    *(indices++) = base + 1;
                    *(indices++) = base + 0;
                }
            }
        }
    }

    if (fill)
    {
        int n = 0;
        for (int i = 0; i < mesh->n_points; i++)
            if (fabs(dot_plane(plane, mesh->points[i])) < epsilon)
                n++;

        int* indices = (int*) malloc(sizeof(int) * n);
        n = 0;
        for (int i = 0; i < mesh->n_points; i++)
            if (fabs(dot_plane(plane, mesh->points[i])) < epsilon)
                indices[n++] = i;

        Point minp = {  1e10,  1e10,  1e10 };
        Point maxp = { -1e10, -1e10, -1e10 };
        for (int i = 0; i < n; i++)
        {
            Point p = mesh->points[indices[i]];
            for (int j = 0; j < 3; j++)
                minp.e[j] = (p.e[j] < minp.e[j] ? p.e[j] : minp.e[j]),
                maxp.e[j] = (p.e[j] > maxp.e[j] ? p.e[j] : maxp.e[j]);
        }

        comp_mesh = mesh;
        comp_middle = mul(add(minp, maxp), 0.5);
        qsort(indices, n, sizeof(int), cut_comparator);

        int middle = reserve_points(mesh, 1);
        mesh->points[middle] = comp_middle;

        int* tri = reserve_triangles(mesh, n);
        for (int i = 0; i < n; i++)
        {
            *(tri++) = middle;
            *(tri++) = indices[i];
            *(tri++) = indices[(i + 1) % n];
        }

        free(indices);
    }
}


void output_mesh_to_stl(FILE* stl, Mesh* mesh)
{
    char header[80] = "Moj STL file, je iza ovoga, da da";
    int32_t triangles = mesh->n_indices / 3;

    fwrite(header, 80, 1, stl);
    fwrite(&triangles, 4, 1, stl);

    for (int i = 0; i < mesh->n_indices; i += 3)
    {
        Point a = mesh->points[mesh->indices[i + 0]];
        Point b = mesh->points[mesh->indices[i + 1]];
        Point c = mesh->points[mesh->indices[i + 2]];
        Point n = norm(cross(sub(b, a), sub(c, a)));
        int16_t attrib = 0;

        float f[12] = {
            (float) n.x, (float) n.y, (float) n.z,
            (float) a.x, (float) a.y, (float) a.z,
            (float) b.x, (float) b.y, (float) b.z,
            (float) c.x, (float) c.y, (float) c.z,
        };

        fwrite(f, 4, 12, stl);
        fwrite(&attrib, 2, 1, stl);
    }
}

void output_mesh_to_obj(FILE* obj, Mesh* mesh, int base_index)
{
    remove_unreferenced_points(mesh);

    for (int i = 0; i < mesh->n_points; i++)
    {
        Point p = mesh->points[i];
        fprintf(obj, "v %f %f %f\n", p.x, p.y, p.z);
    }

    for (int i = 0; i < mesh->n_indices; i += 3)
    {
        int* tri = &mesh->indices[i];
        fprintf(obj, "f %d %d %d\n", base_index + tri[0], base_index + tri[1], base_index + tri[2]);
    }
}


double mesh_z(Mesh* mesh, double x, double y)
{
    double result = 100;
    Point p = { x, y, 0 };

    for (int i = 0; i < mesh->n_indices; i += 3)
    {
        Point a = mesh->points[mesh->indices[i + 0]];
        Point b = mesh->points[mesh->indices[i + 1]];
        Point c = mesh->points[mesh->indices[i + 2]];
        if (a.z > result && b.z > result && c.z > result) continue;

        double aa, ab, ac;
        if ((ac = orientation_xy(a, b, p)) > 0) continue;
        if ((aa = orientation_xy(b, c, p)) > 0) continue;
        if ((ab = orientation_xy(c, a, p)) > 0) continue;

        double inv_area = 1.0 / orientation_xy(a, b, c);
        aa *= inv_area;
        ab *= inv_area;
        ac *= inv_area;

        double z = a.z * aa + b.z * ab + c.z * ac;
        if (z < result)
            result = z;
    }

    return result;
}

#define Min(a, b) ((a) < (b) ? (a) : (b))
#define Max(a, b) ((a) > (b) ? (a) : (b))

struct Z_Buffer
{
    double lx, ly;
    double hx, hy;

    int size_x;
    int size_y;
    int stride;
    double* data;
};

bool in_bounds(Z_Buffer* buffer, double x, double y)
{
    if (x < buffer->lx || x > buffer->hx) return false;
    if (y < buffer->ly || y > buffer->hy) return false;
    return true;
}

void to_pixel(Z_Buffer* buffer, double x, double y, int* px, int* py)
{
    *px = (int)((x - buffer->lx) / (buffer->hx - buffer->lx) * (buffer->size_x - 1) + 0.5);
    *py = (int)((y - buffer->ly) / (buffer->hy - buffer->ly) * (buffer->size_y - 1) + 0.5);
}

void to_xy(Z_Buffer* buffer, double* px, double* py, int x, int y)
{
    *px = x / (double)(buffer->size_x - 1) * (buffer->hx - buffer->lx) + buffer->lx;
    *py = y / (double)(buffer->size_y - 1) * (buffer->hy - buffer->ly) + buffer->ly;
}

void set_z_buffer_bounds(Z_Buffer* buffer, Mesh* mesh)
{
    double lx = +1e10, ly = +1e10;
    double hx = -1e10, hy = -1e10;
    for (int i = 0; i < mesh->n_points; i++)
    {
        double x = mesh->points[i].x;
        double y = mesh->points[i].y;
        if (x < lx) lx = x;
        if (y < ly) ly = y;
        if (x > hx) hx = x;
        if (y > hy) hy = y;
    }

    buffer->lx = lx;
    buffer->ly = ly;
    buffer->hx = hx;
    buffer->hy = hy;
}

void mesh_render_z_buffer(Z_Buffer* buffer, Mesh* mesh)
{
    for (int i = 0; i < mesh->n_indices; i += 3)
    {
        Point a = mesh->points[mesh->indices[i + 0]];
        Point b = mesh->points[mesh->indices[i + 1]];
        Point c = mesh->points[mesh->indices[i + 2]];
        double inv_area = 1.0 / orientation_xy(a, b, c);

        int ax, ay, bx, by, cx, cy;
        to_pixel(buffer, a.x, a.y, &ax, &ay);
        to_pixel(buffer, b.x, b.y, &bx, &by);
        to_pixel(buffer, c.x, c.y, &cx, &cy);

        int x0 = Min(ax, Min(bx, cx));
        int x1 = Max(ax, Max(bx, cx));
        int y0 = Min(ay, Min(by, cy));
        int y1 = Max(ay, Max(by, cy));
        for (int y = y0; y <= y1; y++)
        {
            int index = y * buffer->stride + x0;
            for (int x = x0; x <= x1; x++, index++)
            {
                Point p;
                to_xy(buffer, &p.x, &p.y, x, y);

                double aa = orientation_xy(b, c, p) * inv_area;
                double ab = orientation_xy(c, a, p) * inv_area;
                double ac = orientation_xy(a, b, p) * inv_area;
                if (aa < 0 || ab < 0 || ac < 0) continue;

                double z = a.z * aa + b.z * ab + c.z * ac;
                if (z < buffer->data[index])
                    buffer->data[index] = z;
            }
        }
    }
}
