typedef struct
{
    int n_points;
    Point* points;
    int n_sharp;
    int* i_sharp;
} Contour;

typedef struct
{
    char name[256];
    double length;
    double scale;

    Contour front;
    Contour back;

    int n_cross;
    int n_cross_front;
    int n_cross_back;
    Contour* cross;
} Model;

void load_contour(FILE* kos, Contour* c, int xi, int yi, int zi)
{
    fscanf(kos, "%d%d", &c->n_points, &c->n_sharp);
    c->points  = (Point*) malloc(sizeof(Point) * c->n_points);
    c->i_sharp = (int*)   malloc(sizeof(int) * c->n_sharp);

    for (int i = 0; i < c->n_sharp; i++)
        fscanf(kos, "%d", &c->i_sharp[i]);

    for (int j = 0; j < 3; j++)
        for (int i = 0; i < c->n_points; i++)
        {
            if (xi == j) fscanf(kos, "%lf", &c->points[i].x);
            if (yi == j) fscanf(kos, "%lf", &c->points[i].y);
            if (zi == j) fscanf(kos, "%lf", &c->points[i].z);
        }
}

void load_model(FILE* kos, Model* model)
{
    fscanf(kos, "%s", model->name);
    fscanf(kos, "%lf%lf", &model->length, &model->scale);
    load_contour(kos, &model->back,  0, 1, 2);
    load_contour(kos, &model->front, 0, 1, 2);

    fscanf(kos, "%d%d%d", &model->n_cross, &model->n_cross_front, &model->n_cross_back);
    model->cross = (Contour*) malloc(sizeof(Contour) * model->n_cross);
    for (int i = 0; i < model->n_cross; i++)
    {
        double x;
        fscanf(kos, "%lf", &x);
        load_contour(kos, &model->cross[i], -1, 1, 0);
        for (int j = 0; j < model->cross[i].n_points; j++)
            model->cross[i].points[j].x = x;
    }

    printf(" .. Loaded model \"%s\"\n", model->name);
}


const int INTERPOLATED = 500;  // mail: 350, veliki: 500;

Contour interpolate_contour(Contour* c)
{
    Contour result = { 0 };
    result.n_points = INTERPOLATED;
    interpolate(c->n_points, c->points, c->n_sharp, c->i_sharp, result.n_points, &result.points);
    return result;
}

Contour sub_contour(Contour* c, int from, int to)
{
    Contour r = { 0 };
    r.n_points = to - from;
    r.points = (Point*) malloc(sizeof(Point) * r.n_points);
    memcpy(r.points, c->points + from, sizeof(Point) * r.n_points);

    r.i_sharp = (int*) malloc(sizeof(int) * c->n_sharp);
    for (int i = 0; i < c->n_sharp; i++)
    {
        if (c->i_sharp[i] - 1 < from) continue;
        if (c->i_sharp[i] - 1 >= to) continue;
        r.i_sharp[r.n_sharp++] = c->i_sharp[i] - from;
    }

    return r;
}

void transpose_interpolated(Contour* in, int n_in, Contour out[INTERPOLATED])
{
    for (int i = 0; i < INTERPOLATED; i++)
    {
        Contour c = {};
        c.n_points = n_in;
        c.points = (Point*) calloc(sizeof(Point), n_in);
        for (int j = 0; j < n_in; j++)
            c.points[j] = in[j].points[i];
        out[i] = c;
    }
}

Contour skeg_contour(Contour* c, double y)
{
    int i;
    Point end = c->points[c->n_points - 1];
    for (i = 0; i < c->n_points - 1; i++)
    {
        Point a = c->points[i];
        Point b = c->points[i + 1];
        if (a.y > b.y)
        {
            Point p = a;
            a = b;
            b = p;
        }

        if (a.y <= y && y < b.y)
        {
            double t = (y - a.y) / (b.y - a.y);
            double x = a.x + t * (b.x - a.x);
            double z = a.z + t * (b.z - a.z);
            end = { x, y, z };
            break;
        }
    }

    Contour r = { 0 };
    r.n_points = INTERPOLATED + 3;
    r.points = (Point*) calloc(sizeof(Point), r.n_points);

    int cur = 0;
    while (cur < i)
        r.points[cur++] = c->points[cur];
    while (cur <= INTERPOLATED)
        r.points[cur++] = end;
    r.points[cur++] = { c->points[0].x, end.y, 0 };
    r.points[cur++] = { c->points[0].x, 0, 0 };

    return r;
}
