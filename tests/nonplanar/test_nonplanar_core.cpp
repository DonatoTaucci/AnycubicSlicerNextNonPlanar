#include "NonPlanarCore.hpp"
#include <cassert>
#include <cstdio>
#include <random>

using namespace Slic3r::NonPlanar;

static int failures = 0;
#define CHECK(c) do { if (!(c)) { std::printf("FAIL %s:%d  %s\n", __FILE__, __LINE__, #c); ++failures; } } while (0)
static bool near(double a, double b, double e = 1e-9) { return std::abs(a - b) <= e; }

static Ring square(double x0, double y0, double x1, double y1) { return { {x0,y0},{x1,y0},{x1,y1},{x0,y1} }; }

int main()
{
    // wave
    for (double x = -20; x < 20; x += 0.37)
        for (double y = -20; y < 20; y += 0.41)
            CHECK(std::abs(wave(x, y, 12.)) <= 1. + 1e-12);
    CHECK(near(wave(3., 3., 12.), 1.));
    CHECK(near(wave(0., 5., 12.), 0.));

    CHECK(min_ramp_layers(0.4, 0.2) == 4);
    CHECK(min_ramp_layers(0.5, 0.2) == 5);
    CHECK(min_ramp_layers(0.0, 0.2) == 1);

    // Ramp: column present in layers 2..11
    {
        RampField f;
        f.init(-10, -10, 10, 10, 14, 1.0);
        for (size_t l = 2; l <= 11; ++l) f.set_layer(l, { square(-5, -5, 5, 5) });
        f.finalize(4);
        const double expect[14] = {0,0,.25,.5,.75,1,1,1,1,.75,.5,.25,0,0};
        for (size_t l = 0; l < 14; ++l) CHECK(near(f.factor(l, 0.3, -0.7), expect[l]));
        CHECK(near(f.factor(5, 8., 0.), 0.));     // outside the square
        CHECK(f.factor(99, 0., 0.) == 0.);        // unknown layer
    }
    // Ramp with hole
    {
        RampField f;
        f.init(-10, -10, 10, 10, 10, 0.5);
        for (size_t l = 0; l < 10; ++l) f.set_layer(l, { square(-8, -8, 8, 8), square(-2, -2, 2, 2) });
        f.finalize(2);
        CHECK(near(f.factor(5, 0., 0.), 0.));
        CHECK(near(f.factor(5, 5., 5.), 1.));
    }
    // Random columns: per-layer change of the factor must never exceed 1/N anywhere.
    {
        std::mt19937 rng(42);
        std::uniform_real_distribution<double> U(-15, 15);
        const int N = 5; const size_t L = 60;
        RampField f;
        f.init(-20, -20, 20, 20, L, 0.8);
        for (size_t l = 0; l < L; ++l) {
            Rings r;
            for (int k = 0; k < 3; ++k) { double x = U(rng), y = U(rng); r.push_back(square(x, y, x + 6 + U(rng)/5, y + 6)); }
            f.set_layer(l, r);
        }
        f.finalize(N);
        double worst = 0;
        for (size_t l = 1; l < L; ++l)
            for (double x = -20; x <= 20; x += 0.33)
                for (double y = -20; y <= 20; y += 0.29)
                    worst = std::max(worst, std::abs(f.factor(l, x, y) - f.factor(l - 1, x, y)));
        std::printf("max per-layer factor change: %.4f (limit %.4f)\n", worst, 1. / N);
        CHECK(worst <= 1. / N + 1e-9);

        // Thickness of every extrusion stays within [h - A/N, h + A/N] and >= h/2.
        const double A = 0.5, h = 0.2, Lw = 12.;
        const int Nmin = min_ramp_layers(A, h);
        CHECK(Nmin <= N);
        double tmin = 1e9, tmax = -1e9;
        for (size_t l = 1; l < L; ++l)
            for (double x = -20; x <= 20; x += 0.5)
                for (double y = -20; y <= 20; y += 0.5) {
                    double t = h + A * wave(x, y, Lw) * (f.factor(l, x, y) - f.factor(l - 1, x, y));
                    tmin = std::min(tmin, t); tmax = std::max(tmax, t);
                }
        std::printf("thickness range: [%.3f, %.3f] mm for h=%.2f\n", tmin, tmax, h);
        CHECK(tmin >= 0.5 * h - 1e-9);
    }
    // Horizontal fade: zero at the region edge, full inside, bounded slope, vertical rule kept.
    {
        const int N = 4; const double fade = 2.0; const size_t L = 30;
        RampField f;
        f.init(-20, -20, 20, 20, L, 0.5);
        std::mt19937 rng(7);
        std::uniform_real_distribution<double> U(-3, 3);
        for (size_t l = 1; l < L; ++l) {           // region drifting sideways layer by layer
            double s = 0.4 * double(l) + U(rng);
            f.set_layer(l, { square(-10 + s * 0.3, -8, 6 + s * 0.3, 9) });
        }
        f.finalize(N, fade);
        CHECK(near(f.factor(15, -2.0, 0.5), 1.));               // deep inside
        double edge = f.factor(15, -10 + (0.4 * 15) * 0.3 + 3.2, 9.0); // on the top edge (y = 9)
        CHECK(edge <= 0.25 + 1e-9);
        double worst_v = 0, worst_h = 0;
        for (size_t l = 1; l < L; ++l)
            for (double x = -19; x <= 19; x += 0.21)
                for (double y = -19; y <= 19; y += 0.23) {
                    worst_v = std::max(worst_v, std::abs(f.factor(l, x, y) - f.factor(l - 1, x, y)));
                    worst_h = std::max(worst_h, std::abs(f.factor(l, x + 0.21, y) - f.factor(l, x, y)) / 0.21);
                }
        std::printf("fade: max vertical step %.3f (limit %.3f), max XY slope %.3f /mm (fade limit ~%.3f)\n",
                    worst_v, 1. / N, worst_h, 1. / fade + 1. / (N * 0.5));
        CHECK(worst_v <= 1. / N + 1e-9);
        CHECK(worst_h <= 1. / fade + 1. / (N * f.cell_size()) + 1e-9);
    }
    // Memory budget enlarges the cell
    {
        RampField f;
        f.init(0, 0, 250, 250, 1250, 0.5, size_t(8) << 20);
        CHECK(f.cell_size() > 0.5);
        CHECK(double((250 / f.cell_size() + 2)) * (250 / f.cell_size() + 2) * 1250 <= double(8 << 20) * 1.1);
    }
    // Surface locator on a tilted plane z = 2 + 0.1 x (slope ~5.7 deg) and a steep wall (ignored)
    {
        std::vector<P3> v = {
            {0,0,2},{10,0,3},{10,10,3},   {0,0,2},{10,10,3},{0,10,2},     // CCW seen from above: up-facing
            {0,0,0},{0,10,0},{0,10,5},                                     // vertical, must be ignored
        };
        SurfaceLocator s; s.build(v, 30.);
        double z = 0;
        CHECK(s.query(5, 5, 2.5, 0.2, z) && near(z, 2.5, 1e-9));
        CHECK(s.query(2, 7, 2.3, 0.2, z) && near(z, 2.2, 1e-9));
        CHECK(!s.query(5, 5, 4.0, 0.2, z));     // outside the window
        CHECK(!s.query(15, 5, 2.5, 1.0, z));    // outside the mesh
        SurfaceLocator steep; steep.build(v, 3.);  // 5.7 deg > 3 deg -> nothing kept
        CHECK(steep.empty());
    }
    // warp_polyline subdivides and scales flow
    {
        const double A = 0.4, Lw = 12., h = 0.2;
        auto dz  = [&](double x, double y) { return A * wave(x, y, Lw); };
        auto dzb = [&](double x, double y) { return 0.75 * A * wave(x, y, Lw); };
        auto w = warp_polyline({ {0, 3}, {24, 3} }, Lw / 16., h, dz, dzb);
        CHECK(w.size() == 33);
        CHECK(near(w.back().x, 24.) && near(w.back().y, 3.));
        for (size_t i = 1; i < w.size(); ++i) {
            CHECK(w[i].flow >= 0.5 - 1e-9 && w[i].flow <= 1.5 + 1e-9);
            CHECK(std::hypot(w[i].x - w[i-1].x, w[i].y - w[i-1].y) <= Lw / 16. + 1e-9);
        }
        CHECK(near(w[4].dz, A * wave(3, 3, Lw)));   // x=3 at index 4 (0.75 mm steps)
        auto empty = warp_polyline({}, 1., h, dz, dzb);
        CHECK(empty.empty());
    }
    std::printf(failures ? "%d FAILURES\n" : "all tests passed\n", failures);
    return failures ? 1 : 0;
}
