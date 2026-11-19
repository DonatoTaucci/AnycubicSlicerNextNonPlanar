#include "NonPlanarCore.hpp"

#include <limits>

namespace Slic3r {
namespace NonPlanar {

static constexpr double PI_ = 3.14159265358979323846;

double wave(double x, double y, double wavelength)
{
    if (wavelength <= 0.)
        return 0.;
    const double k = 2. * PI_ / wavelength;
    return std::sin(k * x) * std::sin(k * y);
}

int min_ramp_layers(double amplitude, double min_layer_height, double max_ratio)
{
    if (amplitude <= 0. || min_layer_height <= 0. || max_ratio <= 0.)
        return 1;
    return std::max(1, int(std::ceil(amplitude / (max_ratio * min_layer_height) - 1e-9)));
}

double max_wave_slope(double amplitude, double wavelength)
{
    return wavelength > 0. ? 2. * PI_ * amplitude / wavelength : 0.;
}

// ------------------------------------------------------------------------------------------------
// RampField

void RampField::init(double min_x, double min_y, double max_x, double max_y, size_t layers,
                     double cell_size, size_t max_bytes)
{
    m_data.clear();
    m_layers = layers;
    m_min_x  = min_x;
    m_min_y  = min_y;
    m_cell   = std::max(0.05, cell_size);
    m_ramp   = 1;
    const double w = std::max(0., max_x - min_x), h = std::max(0., max_y - min_y);
    for (;;) {
        m_nx = long(std::ceil(w / m_cell)) + 1;
        m_ny = long(std::ceil(h / m_cell)) + 1;
        if (double(m_nx) * double(m_ny) * double(layers) <= double(max_bytes))
            break;
        m_cell *= 1.25;
    }
    m_data.assign(size_t(m_nx) * size_t(m_ny) * layers, 0);
}

void RampField::set_layer(size_t layer, const Rings &region)
{
    if (layer >= m_layers || m_nx == 0 || m_ny == 0)
        return;
    uint8_t *base = m_data.data() + layer * size_t(m_nx) * size_t(m_ny);
    std::vector<double> xs;
    for (long iy = 0; iy < m_ny; ++iy) {
        const double y = m_min_y + double(iy) * m_cell;
        xs.clear();
        for (const Ring &ring : region) {
            const size_t n = ring.size();
            if (n < 3)
                continue;
            for (size_t i = 0; i < n; ++i) {
                const P2 &p = ring[i], &q = ring[(i + 1) % n];
                // Half-open rule avoids double counting of vertices lying exactly on the scanline.
                if ((p.y <= y && y < q.y) || (q.y <= y && y < p.y))
                    xs.push_back(p.x + (y - p.y) * (q.x - p.x) / (q.y - p.y));
            }
        }
        if (xs.size() < 2)
            continue;
        std::sort(xs.begin(), xs.end());
        uint8_t *row = base + size_t(iy) * size_t(m_nx);
        for (size_t i = 0; i + 1 < xs.size(); i += 2) {
            long ix0 = long(std::ceil((xs[i] - m_min_x) / m_cell));
            long ix1 = long(std::floor((xs[i + 1] - m_min_x) / m_cell));
            ix0 = std::max(ix0, 0L);
            ix1 = std::min(ix1, m_nx - 1);
            for (long ix = ix0; ix <= ix1; ++ix)
                row[ix] = 1;
        }
    }
}

void RampField::finalize(int ramp_layers, double fade_mm)
{
    m_ramp = std::max(1, std::min(ramp_layers, 250));
    const size_t plane = size_t(m_nx) * size_t(m_ny);
    if (plane == 0 || m_layers == 0)
        return;
    const int N = m_ramp;

    // 1) Per layer: horizontal cap from the distance (in mm) to the nearest non-eligible node.
    if (fade_mm > 0.) {
        const double INF = std::numeric_limits<double>::max() / 4.;
        const double d1 = m_cell, d2 = m_cell * std::sqrt(2.);
        std::vector<double> dist(plane);
        for (size_t l = 0; l < m_layers; ++l) {
            uint8_t *v = m_data.data() + l * plane;
            for (size_t c = 0; c < plane; ++c)
                dist[c] = v[c] ? INF : 0.;
            auto at = [&](long ix, long iy) -> double {
                // Nodes outside the grid count as non-eligible: distance to them is one cell.
                return (ix < 0 || iy < 0 || ix >= m_nx || iy >= m_ny) ? 0. : dist[size_t(iy) * size_t(m_nx) + size_t(ix)];
            };
            // Two-pass chamfer distance transform.
            for (long iy = 0; iy < m_ny; ++iy)
                for (long ix = 0; ix < m_nx; ++ix) {
                    double &d = dist[size_t(iy) * size_t(m_nx) + size_t(ix)];
                    if (d == 0.) continue;
                    d = std::min({ d, at(ix - 1, iy) + d1, at(ix, iy - 1) + d1, at(ix - 1, iy - 1) + d2, at(ix + 1, iy - 1) + d2 });
                }
            for (long iy = m_ny - 1; iy >= 0; --iy)
                for (long ix = m_nx - 1; ix >= 0; --ix) {
                    double &d = dist[size_t(iy) * size_t(m_nx) + size_t(ix)];
                    if (d == 0.) continue;
                    d = std::min({ d, at(ix + 1, iy) + d1, at(ix, iy + 1) + d1, at(ix + 1, iy + 1) + d2, at(ix - 1, iy + 1) + d2 });
                }
            for (size_t c = 0; c < plane; ++c)
                v[c] = v[c] ? uint8_t(std::min<double>(N, std::floor(double(N) * dist[c] / fade_mm + 1e-9))) : 0;
        }
    } else {
        for (uint8_t &v : m_data)
            v = v ? uint8_t(N) : 0;
    }

    // 2) Vertical: g_l = min_k (v_k + |l - k|), i.e. the factor may change by at most one step
    //    per layer. Layers outside the object count as 0.
    for (size_t c = 0; c < plane; ++c) {
        int prev = 0;
        for (size_t l = 0; l < m_layers; ++l) {
            uint8_t &v = m_data[l * plane + c];
            prev = std::min<int>(v, prev + 1);
            v = uint8_t(prev);
        }
        prev = 0;
        for (size_t l = m_layers; l-- > 0;) {
            uint8_t &v = m_data[l * plane + c];
            prev = std::min<int>(v, prev + 1);
            v = uint8_t(prev);
        }
    }
}

double RampField::count_at(size_t layer, long ix, long iy) const
{
    if (ix < 0 || iy < 0 || ix >= m_nx || iy >= m_ny)
        return 0.;
    return double(m_data[layer * size_t(m_nx) * size_t(m_ny) + size_t(iy) * size_t(m_nx) + size_t(ix)]);
}

double RampField::factor(size_t layer, double x, double y) const
{
    if (layer >= m_layers || m_nx == 0 || m_ny == 0)
        return 0.;
    const double fx = (x - m_min_x) / m_cell, fy = (y - m_min_y) / m_cell;
    const long   ix = long(std::floor(fx)), iy = long(std::floor(fy));
    const double tx = fx - double(ix), ty = fy - double(iy);
    const double c00 = count_at(layer, ix, iy), c10 = count_at(layer, ix + 1, iy);
    const double c01 = count_at(layer, ix, iy + 1), c11 = count_at(layer, ix + 1, iy + 1);
    const double c = (c00 * (1. - tx) + c10 * tx) * (1. - ty) + (c01 * (1. - tx) + c11 * tx) * ty;
    return std::min(1., c / double(m_ramp));
}

// ------------------------------------------------------------------------------------------------
// SurfaceLocator

void SurfaceLocator::build(const std::vector<P3> &v, double max_slope_deg, double cell_size)
{
    m_tris.clear();
    m_bins.clear();
    m_cell = std::max(0.1, cell_size);
    const double min_nz = std::cos(std::max(0., std::min(89.9, max_slope_deg)) * PI_ / 180.);
    double minx = std::numeric_limits<double>::max(), miny = minx;
    double maxx = -minx, maxy = -minx;
    for (size_t i = 0; i + 2 < v.size(); i += 3) {
        const P3 &a = v[i], &b = v[i + 1], &c = v[i + 2];
        const double ux = b.x - a.x, uy = b.y - a.y, uz = b.z - a.z;
        const double wx = c.x - a.x, wy = c.y - a.y, wz = c.z - a.z;
        const double nx = uy * wz - uz * wy, ny = uz * wx - ux * wz, nz = ux * wy - uy * wx;
        const double len = std::sqrt(nx * nx + ny * ny + nz * nz);
        if (len <= 0. || nz / len < min_nz)
            continue;
        m_tris.push_back({ a, b, c });
        minx = std::min({ minx, a.x, b.x, c.x });
        miny = std::min({ miny, a.y, b.y, c.y });
        maxx = std::max({ maxx, a.x, b.x, c.x });
        maxy = std::max({ maxy, a.y, b.y, c.y });
    }
    if (m_tris.empty())
        return;
    m_min_x = minx;
    m_min_y = miny;
    m_nx = long(std::floor((maxx - minx) / m_cell)) + 1;
    m_ny = long(std::floor((maxy - miny) / m_cell)) + 1;
    m_bins.assign(size_t(m_nx) * size_t(m_ny), {});
    for (uint32_t t = 0; t < uint32_t(m_tris.size()); ++t) {
        const Tri &tr = m_tris[t];
        const long ix0 = long(std::floor((std::min({ tr.a.x, tr.b.x, tr.c.x }) - m_min_x) / m_cell));
        const long ix1 = long(std::floor((std::max({ tr.a.x, tr.b.x, tr.c.x }) - m_min_x) / m_cell));
        const long iy0 = long(std::floor((std::min({ tr.a.y, tr.b.y, tr.c.y }) - m_min_y) / m_cell));
        const long iy1 = long(std::floor((std::max({ tr.a.y, tr.b.y, tr.c.y }) - m_min_y) / m_cell));
        for (long iy = std::max(0L, iy0); iy <= std::min(m_ny - 1, iy1); ++iy)
            for (long ix = std::max(0L, ix0); ix <= std::min(m_nx - 1, ix1); ++ix)
                m_bins[size_t(iy) * size_t(m_nx) + size_t(ix)].push_back(t);
    }
}

bool SurfaceLocator::query(double x, double y, double z_ref, double window, double &z_out) const
{
    if (m_tris.empty())
        return false;
    const long ix = long(std::floor((x - m_min_x) / m_cell));
    const long iy = long(std::floor((y - m_min_y) / m_cell));
    if (ix < 0 || iy < 0 || ix >= m_nx || iy >= m_ny)
        return false;
    bool   found = false;
    double best  = std::numeric_limits<double>::max();
    for (uint32_t t : m_bins[size_t(iy) * size_t(m_nx) + size_t(ix)]) {
        const Tri &tr = m_tris[t];
        // Barycentric coordinates in XY.
        const double d = (tr.b.y - tr.c.y) * (tr.a.x - tr.c.x) + (tr.c.x - tr.b.x) * (tr.a.y - tr.c.y);
        if (std::abs(d) < 1e-12)
            continue;
        const double l1 = ((tr.b.y - tr.c.y) * (x - tr.c.x) + (tr.c.x - tr.b.x) * (y - tr.c.y)) / d;
        const double l2 = ((tr.c.y - tr.a.y) * (x - tr.c.x) + (tr.a.x - tr.c.x) * (y - tr.c.y)) / d;
        const double l3 = 1. - l1 - l2;
        const double eps = 1e-9;
        if (l1 < -eps || l2 < -eps || l3 < -eps)
            continue;
        const double z = l1 * tr.a.z + l2 * tr.b.z + l3 * tr.c.z;
        const double dist = std::abs(z - z_ref);
        if (dist <= window && dist < best) {
            best  = dist;
            z_out = z;
            found = true;
        }
    }
    return found;
}

} // namespace NonPlanar
} // namespace Slic3r
