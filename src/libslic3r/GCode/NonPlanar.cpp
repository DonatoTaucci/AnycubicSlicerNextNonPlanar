#include "NonPlanar.hpp"

#include "../ClipperUtils.hpp"
#include "../Layer.hpp"
#include "../Model.hpp"
#include "../Print.hpp"
#include "../Surface.hpp"

namespace Slic3r {

NonPlanarSettings NonPlanarSettings::from(const PrintObjectConfig &cfg)
{
    NonPlanarSettings s;
    s.infill        = cfg.nonplanar_infill.value;
    s.top_surface   = cfg.nonplanar_top_surface.value;
    s.amplitude     = std::max(0., cfg.nonplanar_amplitude.value);
    s.wavelength    = std::max(0., cfg.nonplanar_wavelength.value);
    s.ramp_layers   = std::max(1, cfg.nonplanar_ramp_layers.value);
    s.top_ratio     = std::clamp(cfg.nonplanar_top_max_offset.value / 100., 0., 1.);
    s.top_max_angle = std::clamp(cfg.nonplanar_top_max_angle.value, 0., 89.);
    return s;
}

double NonPlanarSettings::travel_clearance(double max_layer_height) const
{
    double c = 0.;
    if (this->waves())
        c = std::max(c, amplitude);
    if (top_surface)
        c = std::max(c, top_ratio * max_layer_height);
    return c > 0. ? c + 0.1 : 0.;
}

static void append_rings(const ExPolygons &expolys, NonPlanar::Rings &out)
{
    auto add = [&out](const Polygon &poly) {
        NonPlanar::Ring r;
        r.reserve(poly.points.size());
        for (const Point &p : poly.points)
            r.push_back({ unscale<double>(p.x()), unscale<double>(p.y()) });
        out.emplace_back(std::move(r));
    };
    for (const ExPolygon &ex : expolys) {
        add(ex.contour);
        for (const Polygon &h : ex.holes)
            add(h);
    }
}

NonPlanarObject::NonPlanarObject(const PrintObject &object)
    : m_settings(NonPlanarSettings::from(object.config()))
{
    if (!m_settings.any() || object.layer_count() == 0)
        return;

    const auto layers = object.layers();

    if (m_settings.waves()) {
        BoundingBox bbox;
        double      min_h = std::numeric_limits<double>::max();
        for (const Layer *layer : layers) {
            if (!layer->lslices.empty())
                bbox.merge(get_extents(layer->lslices));
            if (layer->height > 0.)
                min_h = std::min(min_h, layer->height);
        }
        if (bbox.defined && min_h < std::numeric_limits<double>::max()) {
            m_field.init(unscale<double>(bbox.min.x()), unscale<double>(bbox.min.y()),
                         unscale<double>(bbox.max.x()), unscale<double>(bbox.max.y()), layers.size(), 0.5);
            for (size_t i = 0; i < layers.size(); ++i) {
                const Layer *layer = layers[i];
                if (i == 0)
                    continue; // never on the first layer
                ExPolygons eligible;
                for (const LayerRegion *layerm : layer->regions())
                    for (const Surface &s : layerm->fill_surfaces.surfaces)
                        if (s.surface_type == stInternal)
                            eligible.push_back(s.expolygon);
                if (eligible.empty())
                    continue;
                NonPlanar::Rings rings;
                append_rings(union_ex(eligible), rings);
                m_field.set_layer(i, rings);
            }
            const int ramp = std::max(m_settings.ramp_layers, NonPlanar::min_ramp_layers(m_settings.amplitude, min_h, 0.5));
            m_field.finalize(ramp, m_settings.fade);
        }
    }

    if (m_settings.top_surface && object.model_object() != nullptr) {
        std::vector<NonPlanar::P3> tris;
        const Transform3d trafo_obj = object.trafo_centered();
        for (const ModelVolume *volume : object.model_object()->volumes) {
            if (!volume->is_model_part())
                continue;
            const Transform3d trafo = trafo_obj * volume->get_matrix();
            const indexed_triangle_set &its = volume->mesh().its;
            const bool flipped = trafo.matrix().block<3, 3>(0, 0).determinant() < 0.;
            for (const stl_triangle_vertex_indices &f : its.indices) {
                Vec3d v[3];
                for (int k = 0; k < 3; ++k)
                    v[k] = trafo * its.vertices[f[k]].cast<double>();
                if (flipped)
                    std::swap(v[1], v[2]);
                for (int k = 0; k < 3; ++k)
                    tris.push_back({ v[k].x(), v[k].y(), v[k].z() });
            }
        }
        m_surface.build(tris, m_settings.top_max_angle, 2.0);
    }
}

NonPlanarObject::Mode NonPlanarObject::mode_for(ExtrusionRole role, const Layer &layer) const
{
    if (layer.id() == 0)
        return Mode::None;
    switch (role) {
    case erInternalInfill:
        return m_settings.infill && !m_field.empty() ? Mode::Wave : Mode::None;
    case erTopSolidInfill:
    case erIroning:
        return m_settings.top_surface && !m_surface.empty() ? Mode::Top : Mode::None;
    default:
        return Mode::None;
    }
}

double NonPlanarObject::wave_dz(size_t layer_idx, double x, double y) const
{
    if (!m_settings.waves() || m_field.empty())
        return 0.;
    const double f = m_field.factor(layer_idx, x, y);
    return f == 0. ? 0. : m_settings.amplitude * f * NonPlanar::wave(x, y, m_settings.wavelength);
}

double NonPlanarObject::top_dz(const Layer &layer, double x, double y) const
{
    if (m_surface.empty())
        return 0.;
    // Layer top in the slicing (mesh) frame.
    const double top   = layer.slice_z + 0.5 * layer.height;
    const double limit = m_settings.top_ratio * layer.height;
    double z = 0.;
    if (!m_surface.query(x, y, top, layer.height, z))
        return 0.;
    return std::clamp(z - top, -limit, limit);
}

} // namespace Slic3r
