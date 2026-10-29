# Momentum conversion: technical reference

The k conversion of maps and cuts, and the photon-energy (k_z) scans: `kspace.py`, `cutk.py`, `kzmap.py`, `kzconv.py` in this folder.

## k conversion

For a `map`, the analyser's own slit-direction axis (`scan.k`) is already
delivered in Å⁻¹ by the MBS acquisition software. The deflector axis
(`scan.x`) is a raw angle and is converted on demand:

```
k_parallel [A^-1] = 0.5123 * sqrt(KE[eV]) * sin(angle[rad])
```

using `center_ke` from the file as the reference kinetic energy. This ignores
the inner-potential/refraction correction and treats KE as constant across
the energy window -- fine for a narrow window around a feature, worth
checking for a wide one (pass an explicit `kinetic_energy_eV` array, one
value per energy channel, to `to_kspace_cube()` for per-channel accuracy).

## kz map processing

A photon-energy scan arrives stacked as measured, because nothing in the
files says where each spectrum's Fermi level actually is. **Functions → kz
map processing...** on a kz map's viewer measures it from the spectra
themselves:

1. The slit cut opens alongside with a selection box on it. Drag the box
   over a Fermi edge -- the angle range says which channels are summed into
   the EDC, the energy range is what the fit runs over. The box starts
   around wherever the edge appears to be in the scan as a whole.
2. That same box, **by detector index**, is used on every spectrum. An index
   range is the same channels and the same analyser window everywhere; an
   energy range would mean something different in each, since these are
   exactly the spectra whose energy scales do not yet agree.
3. Each spectrum is shifted to put its own edge at zero, and the cube is
   cropped to the energy range they all still cover.
4. **Normalise each spectrum by its total intensity** (ticked by default)
   divides the beamline out: the photon flux and analyser transmission vary
   by a large factor across a wide scan, so raw counts in a kz map are
   largely a picture of the beamline. Done after cropping, so every total is
   a sum over the same range.

E_F against photon energy is plotted as it goes. A smooth curve means the
fits are sound; scatter means some of them found something other than the
edge. Spectra whose fit is not believed -- the level pinned to the end of
the window, or a step height no larger than its own error bar -- are marked
and filled in from their neighbours rather than dropped.

The fitted levels are stored on the result as `kz.fermi_level_eV`, so the
calibration can be checked or reused later.

## kz -> momentum

**Functions → Convert to momentum (kz, k∥)...** on a kz map's viewer converts it to `(k_z, k_par, E)`
in Å⁻¹, with k_z along the first axis.

The conversion is a coordinate change:

```
k_par = sqrt(A E_kin) sin(alpha)
k_z   = sqrt(A m* (E_kin + V0) - k_par^2)      A = 0.262466 Å⁻² eV⁻¹
```

It is done by **inverting** the map rather than projecting it: for each
target `(k_z, k_par)` there is exactly one photon energy and one emission
angle that produced it, so the conversion is a resampling of the original
regular grid — one `map_coordinates` call per energy slice. Exact for every
slice, and about 120× faster than interpolating a scattered cloud (2.8 s for
a real 41×664×602 cube). Points the measurement never reached come back as
NaN rather than as the nearest sample smeared outwards.

**Flatten the Fermi surface first.** The conversion reads the energy axis
literally, so an edge that bends across the analyser angle becomes a bend in
k_z that looks like dispersion. Open the slit cut, use Fermi-surface
correction, then convert. The dialog measures the bend and says so; it is a
reminder, not a refusal.

**Choosing V₀.** The inner potential is not measured by the experiment, so
the window is built around choosing it:

- the preview reconverts one energy slice as you drag, in milliseconds;
- the zone boundaries are drawn over it from the space group, lattice
  constants, the surface normal and the in-plane direction the slit ran
  along: the real cut of the 3-D zone by the (k_z, k_par) plane, for any
  surface, not only (001). The converted map's **Brillouin zone...** entry
  draws the same cut on its contour;
- **Scan V0** converts one slice at each of a range of inner potentials,
  measures the k_z period each gives, and reports where it equals the
  lattice's. It also reports how well that pins V₀ — for a scan covering a
  couple of zones, only to a few eV, which is the measurement's limit rather
  than the method's;
- **Measure a k_z period...** lets you click two points that are the same
  feature one zone apart. Their separation is matched against every
  reciprocal lattice vector within 15% and the matching planes are listed —
  which says how the crystal cleaved. Centring is handled: in a body-centred
  lattice the period along [001] is 4π/a, not 2π/a, and the answer says so.
