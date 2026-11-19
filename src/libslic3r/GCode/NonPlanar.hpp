// Non-planar extrusion support for the G-code generator.
// Builds, per PrintObject, the data needed by NonPlanarCore and answers Z-offset queries
// while GCode::_extrude() emits the moves. Because the offsets are applied while the G-code is
// generated (not in a post-processing script), the G-code preview shows the real 3D paths.

#pragma once

#include "NonPlanarCore.hpp"
#include "../ExtrusionEntity.hpp"

#include <memory>

namespace Slic3r {

class Layer;
class PrintObject;
class PrintObjectConfig;

struct NonPlanarSettings
{
    bool   infill        { false }; // interlocking waves in sparse infill
    bool   top_surface   { false }; // top surfaces follow the model surface
    double amplitude     { 0.4 };   // mm
    double fade          { 2.0 };   // mm, horizontal fade-in from the infill boundary
    double wavelength    { 12. };   // mm
    int    ramp_layers   { 4 };
    double top_ratio     { 0.5 };   // max |dz| on top surfaces as a fraction of the layer height
    double top_max_angle { 20. };   // deg, steeper surfaces stay planar

    static NonPlanarSettings from(const PrintObjectConfig &cfg);
    bool waves() const { return infill && amplitude > 0. && wavelength > 0.; }
    bool any()   const { return waves() || top_surface; }
    // Height above the nominal layer Z that travel moves must keep to clear warped extrusions.
    double travel_clearance(double max_layer_height) const;
};

class NonPlanarObject
{
public:
    explicit NonPlanarObject(const PrintObject &object);

    const NonPlanarSettings &settings() const { return m_settings; }
    bool                     active()   const { return m_settings.any(); }

    enum class Mode { None, Wave, Top };
    Mode mode_for(ExtrusionRole role, const Layer &layer) const;

    // Z offset of the wave field at the given object layer index and object XY [mm].
    double wave_dz(size_t layer_idx, double x, double y) const;
    // Z offset that makes a top surface follow the model, for the given layer.
    double top_dz(const Layer &layer, double x, double y) const;

    // Max segment length used to subdivide warped paths [mm].
    double max_segment() const { return std::max(0.2, m_settings.wavelength / 16.); }
    int    effective_ramp_layers() const { return m_field.ramp_layers(); }

private:
    NonPlanarSettings        m_settings;
    NonPlanar::RampField     m_field;
    NonPlanar::SurfaceLocator m_surface;
};

} // namespace Slic3r
