// Non-planar extrusion core (std-only, unit-testable without the rest of libslic3r).
//
// Two independent effects are provided:
//
//  1. "Interlocking waves" (strength): sparse infill extrusions are displaced in Z by
//     dz = A * f(layer, x, y) * w(x, y)
//     where w is an egg-crate wave  sin(2*pi*x/L) * sin(2*pi*y/L)  shared by all layers and
//     f in [0, 1] is a ramp factor that grows by at most 1/N per layer when entering an
//     eligible column and decreases by at most 1/N per layer when leaving it. Consecutive
//     layers therefore nest into each other (shear interlock) while every extrusion keeps a
//     thickness of  h + (f_i - f_{i-1}) * A * w  >= h - A/N  (N is chosen so that this stays
//     >= h/2). Flow is scaled with the local thickness.
//
//  2. "Smooth top surfaces" (aesthetics, a.k.a. Z anti-aliasing): the exposed top surface of a
//     layer is printed following the real mesh surface instead of the flat layer plane,
//     clamped to +-(ratio * h) around the nominal layer top. This removes the stair-stepping on
//     gently sloped tops. Flow is scaled with the local thickness.
//
// All coordinates are in millimetres. XY are object coordinates (the same frame for every
// instance of the object), Z in the "mesh frame" used for slicing.

#pragma once

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace Slic3r {
namespace NonPlanar {

struct P2 { double x, y; };
struct P3 { double x, y, z; };

using Ring  = std::vector<P2>;    // closed polygon, last point != first point
using Rings = std::vector<Ring>;  // filled with the even-odd rule (contours and holes)

// Egg-crate wave in [-1, 1].
double wave(double x, double y, double wavelength);

// Per-layer ramp factor over a regular XY grid.
class RampField
{
public:
    // cell_size: desired raster resolution [mm]. max_bytes: memory budget; the cell size is
    // enlarged if needed so that cells * layers <= max_bytes.
    void init(double min_x, double min_y, double max_x, double max_y, size_t layers,
              double cell_size = 1.0, size_t max_bytes = size_t(64) << 20);
    // Mark the eligible region of a layer (even-odd fill of the rings). May be called once per layer.
    void set_layer(size_t layer, const Rings &region);
    // Convert the occupancy into ramp counts. ramp_layers >= 1.
    // fade_mm: horizontal distance from the region boundary over which the factor grows from 0
    // to 1, so that warped paths meet the neighbouring planar extrusions (walls, solid infill)
    // at dz = 0. The result is 1-Lipschitz in the layer direction (|f_i - f_{i-1}| <= 1/N)
    // and Lipschitz in XY with constant 1/fade_mm.
    void finalize(int ramp_layers, double fade_mm = 0.);

    // Ramp factor in [0, 1], bilinearly interpolated. 0 outside the grid or for unknown layers.
    double factor(size_t layer, double x, double y) const;

    bool   empty()       const { return m_layers == 0 || m_nx == 0 || m_ny == 0; }
    double cell_size()   const { return m_cell; }
    size_t layer_count() const { return m_layers; }
    int    ramp_layers() const { return m_ramp; }

private:
    double  count_at(size_t layer, long ix, long iy) const;

    double  m_min_x { 0 }, m_min_y { 0 }, m_cell { 1. };
    long    m_nx { 0 }, m_ny { 0 };
    size_t  m_layers { 0 };
    int     m_ramp { 1 };
    // occupancy (0/1) before finalize(), ramp count (0..ramp) after.
    std::vector<uint8_t> m_data;
};

// Finds the mesh surface closest to a reference height, looking straight down/up from (x, y).
class SurfaceLocator
{
public:
    // Triangles given as 3 consecutive points. Only up-facing triangles whose slope is below
    // max_slope_deg are kept (steeper parts are left planar).
    void build(const std::vector<P3> &triangle_vertices, double max_slope_deg, double cell_size = 2.0);
    // Returns true and the surface Z closest to z_ref within [z_ref - window, z_ref + window].
    bool query(double x, double y, double z_ref, double window, double &z_out) const;
    bool empty() const { return m_tris.empty(); }

private:
    struct Tri { P3 a, b, c; };
    std::vector<Tri>                  m_tris;
    std::vector<std::vector<uint32_t>> m_bins;
    double m_min_x { 0 }, m_min_y { 0 }, m_cell { 2. };
    long   m_nx { 0 }, m_ny { 0 };
};

struct WarpedPoint {
    double x, y;
    double dz;    // Z offset of this point relative to the nominal layer top
    double flow;  // flow multiplier of the segment ending at this point (1 for the first point)
};

// Subdivide a polyline so that no segment is longer than max_segment and assign to every point
// the Z offset dz(x, y). The flow multiplier of each segment is
//      (h + dz(mid) - dz_below(mid)) / h  clamped to [min_flow, max_flow].
template<class DzThis, class DzBelow>
std::vector<WarpedPoint> warp_polyline(const std::vector<P2> &pts, double max_segment, double h,
                                       DzThis &&dz, DzBelow &&dz_below,
                                       double min_flow = 0.25, double max_flow = 2.0)
{
    std::vector<WarpedPoint> out;
    if (pts.empty())
        return out;
    out.reserve(pts.size() * 2);
    out.push_back({ pts.front().x, pts.front().y, dz(pts.front().x, pts.front().y), 1. });
    for (size_t i = 1; i < pts.size(); ++i) {
        const P2 &a = pts[i - 1], &b = pts[i];
        const double dx = b.x - a.x, dy = b.y - a.y;
        const double len = std::sqrt(dx * dx + dy * dy);
        if (len <= 0.)
            continue;
        const size_t n = std::max<size_t>(1, size_t(std::ceil(len / std::max(1e-3, max_segment))));
        for (size_t k = 1; k <= n; ++k) {
            const double t  = double(k) / double(n);
            const double tm = (double(k) - 0.5) / double(n);
            const double x = a.x + dx * t, y = a.y + dy * t;
            const double xm = a.x + dx * tm, ym = a.y + dy * tm;
            double flow = (h + dz(xm, ym) - dz_below(xm, ym)) / h;
            flow = std::min(max_flow, std::max(min_flow, flow));
            out.push_back({ x, y, dz(x, y), flow });
        }
    }
    return out;
}

// Number of ramp layers needed so that the per-layer thickness change A/N stays <= max_ratio * h.
int min_ramp_layers(double amplitude, double min_layer_height, double max_ratio = 0.5);

// Largest Z offset of the nozzle tip when following the wave relative to the material under its
// flat face: slope * nozzle_radius. Used to warn about collisions.
double max_wave_slope(double amplitude, double wavelength);

} // namespace NonPlanar
} // namespace Slic3r
