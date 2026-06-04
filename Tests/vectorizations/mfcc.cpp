void foo(const double *a, const double *b, double *out, int N) {
    double s = 0.0;

    for (int i = 0; i < N; ++i)
        s += a[i] * b[i];

    *out = s;
}
