#!/usr/bin/env julia
# Generates a synthetic Gaussian seamount bathymetry for CI testing.
# Output: bathymetry/balanus-GMRT-bathymetry-preprocessed.nc

using NCDatasets

# Parameters — must match seamount.jl CLI defaults: --H=100, --FWHM=500
const H      = 100.0    # peak height (m)
const FWHM   = 500.0    # target FWHM (m)
const σ      = FWHM / (2 * sqrt(2 * log(2)))   # Gaussian σ derived from FWHM

# Grid: fine enough for accurate discrete FWHM measurement (rtol < 1e-3)
const dx     = 15.0     # grid spacing (m)
const extent = 1500.0   # half-extent of domain (m)

x = collect(-extent:dx:extent)
y = collect(-extent:dx:extent)

# Gaussian seamount centred at origin: elevation[i,j] = f(x[i], y[j])
elevation = [H * exp(-(xi^2 + yj^2) / (2σ^2)) for xi in x, yj in y]

# Measure FWHM with the same formula used by measure_FWHM() in utils.jl
Δx = diff(x)[1]
Δy = diff(y)[1]
area_at_HM    = sum(elevation .> H/2) * Δx * Δy
measured_FWHM = 2 * sqrt(area_at_HM / π)

@assert isapprox(measured_FWHM, FWHM, rtol=0.01) "FWHM mismatch: got $measured_FWHM, expected $FWHM"
@info "Synthetic Gaussian seamount ready" H measured_FWHM

outpath = joinpath(@__DIR__, "..", "bathymetry", "balanus-GMRT-bathymetry-preprocessed.nc")
NCDatasets.NCDataset(outpath, "c") do ds
    NCDatasets.defDim(ds, "x", length(x))
    NCDatasets.defDim(ds, "y", length(y))

    xv = NCDatasets.defVar(ds, "x", Float64, ("x",))
    yv = NCDatasets.defVar(ds, "y", Float64, ("y",))
    ev = NCDatasets.defVar(ds, "periodic_elevation", Float64, ("x", "y"))

    xv[:] = x
    yv[:] = y
    ev[:, :] = elevation

    # Attributes read by seamount.jl
    ds.attrib["FWHM"] = measured_FWHM
    ds.attrib["H"]    = H
end

@info "Saved synthetic bathymetry to $outpath"
