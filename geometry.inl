typedef union
{
    struct { double x, y, z; };
    double e[3];
} Point;


typedef struct
{
    double a, b, c, d, t;
} Segment;

typedef struct
{
    int n;
    Segment* seg;
    double t0;
    double t1;
} Spline;

void make_spline_from_knots(Spline* spline, int knot_count, double* values, double* times)
{
    assert(knot_count >= 2);
    int segment_count = knot_count - 1;

    Segment* segments = (Segment*) calloc(sizeof(Segment), segment_count);

    // compute natural cubic spline segments
    {
        int n = segment_count;
        double* t = times;
        double* a = values;

        double* u = (double*) calloc(sizeof(double), n);
        double* z = (double*) calloc(sizeof(double), n);
        double last_da = a[1] - a[0];
        double last_dt = t[1] - t[0];
        for (int i = 1; i < n; i++)
        {
            double dt = t[i + 1] - t[i];
            double da = a[i + 1] - a[i];
            double l = 2 * (t[i + 1] - t[i - 1]) - last_dt * u[i - 1];
            u[i] = dt / l;
            z[i] = (3 * (da / dt - last_da / last_dt) - last_dt * z[i - 1]) / l;
            last_da = da;
            last_dt = dt;
        }

        double right_c = 0.0f;
        for (int i = n - 1; i >= 0; i--)
        {
            double c = z[i] - u[i] * right_c;
            double da = a[i + 1] - a[i];
            double dc = right_c - c;
            double dt = t[i + 1] - t[i];
            segments[i].a = a[i];
            segments[i].b = da / dt - (right_c + 2 * c) * dt / 3;
            segments[i].c = c;
            segments[i].d = dc / dt / 3;
            segments[i].t = t[i];
            right_c = c;
        }

        free(u);
        free(z);
    }

    spline->n = segment_count;
    spline->seg = segments;
    spline->t0 = times[0];
    spline->t1 = times[segment_count];
}

int get_segment_at_time(Spline* spline, double t)
{
    if (t <= spline->t0) return 0;
    if (t >= spline->t1) return spline->n - 1;

    int i = 0;
    while (i + 1 < spline->n && t >= spline->seg[i + 1].t)
        i++;
    return i;
}

double spline_x(Spline* spline, double t)
{
    int i = get_segment_at_time(spline, t);
    Segment s = spline->seg[i];
    double dt = t - s.t;
    return s.a + (s.b + (s.c + s.d * dt) * dt) * dt;
}


void interpolate(int n_points, Point* points, int n_sharp, int* i_sharp, int out_n_points, Point** out_points)
{
    int capacity = n_points * 5;
    double* xs = (double*) malloc(sizeof(double) * capacity);
    double* ys = (double*) malloc(sizeof(double) * capacity);
    double* zs = (double*) malloc(sizeof(double) * capacity);
    double* ts = (double*) malloc(sizeof(double) * capacity);

    double t_epsilon = 0.01f;
    int n = 0;
    double t_cursor = 0;
    for (int i = 0; i < n_points; i++)
    {
        if (i > 0)
        {
            double dx = points[i].x - points[i - 1].x;
            double dy = points[i].y - points[i - 1].y;
            double dz = points[i].z - points[i - 1].z;
            double d = sqrt(dx * dx + dy * dy + dz * dz);
            if (d < t_epsilon)
                d = t_epsilon;
            t_cursor += d;
        }

        bool is_sharp = false;
        for (int j = 0; j < n_sharp; j++)
            if (i_sharp[j] == i + 1)
                is_sharp = true;

        if (is_sharp && i > 0)
        {
            for (int j = 2; j >= 1; j--)
            {
                xs[n] = points[i].x + 0.05 * j * (points[i - 1].x - points[i].x);
                ys[n] = points[i].y + 0.05 * j * (points[i - 1].y - points[i].y);
                zs[n] = points[i].z + 0.05 * j * (points[i - 1].z - points[i].z);
                ts[n] = t_cursor;
                t_cursor += t_epsilon;
                n++;
            }
        }

        xs[n] = points[i].x;
        ys[n] = points[i].y;
        zs[n] = points[i].z;
        ts[n] = t_cursor;
        n++;

        if (is_sharp && i + 1 < n_points)
        {
            for (int j = 1; j <= 2; j++)
            {
                t_cursor += t_epsilon;
                xs[n] = points[i].x + 0.05 * j * (points[i + 1].x - points[i].x);
                ys[n] = points[i].y + 0.05 * j * (points[i + 1].y - points[i].y);
                zs[n] = points[i].z + 0.05 * j * (points[i + 1].z - points[i].z);
                ts[n] = t_cursor;
                n++;
            }
        }
    }

    Spline x_spline = { 0 };
    Spline y_spline = { 0 };
    Spline z_spline = { 0 };
    make_spline_from_knots(&x_spline, n, xs, ts);
    make_spline_from_knots(&y_spline, n, ys, ts);
    make_spline_from_knots(&z_spline, n, zs, ts);

    *out_points = (Point*) malloc(sizeof(Point) * out_n_points);
    for (int i = 0; i < out_n_points; i++)
    {
        double t = (i / (double)(out_n_points - 1)) * t_cursor;
        double x = spline_x(&x_spline, t);
        double y = spline_x(&y_spline, t);
        double z = spline_x(&z_spline, t);
        (*out_points)[i] = { x, y, z };
    }

    free(xs);
    free(ys);
    free(zs);
    free(ts);
}


typedef struct
{
    Point p;
    Point n;
} Plane;

Point  add(Point a, Point  b) { return { a.x + b.x, a.y + b.y, a.z + b.z }; }
Point  sub(Point a, Point  b) { return { a.x - b.x, a.y - b.y, a.z - b.z }; }
Point  mul(Point a, double b) { return { a.x * b,   a.y * b,   a.z * b   }; }
double dot(Point a, Point  b) { return a.x * b.x + a.y * b.y + a.z * b.z;   }
double len(Point a) { return dot(a, a); }
Point norm(Point a) { return mul(a, 1.0 / dot(a, a)); }
double dot_plane(Plane p, Point t) { return dot(p.n, sub(t, p.p)); }

Point cross(Point a, Point b)
{
    return { (a.y * b.z) - (a.z * b.y),
             (a.z * b.x) - (a.x * b.z),
             (a.x * b.y) - (a.y * b.x) };
}

double orientation_xy(Point a, Point b, Point c)
{
    return (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x);
}

Point intersect_line_plane(Plane p, Point a, Point b)
{
    Point dir = norm(sub(b, a));
    double t = dot(p.n, sub(p.p, a)) / dot(p.n, dir);
    return add(a, mul(dir, t));
}

