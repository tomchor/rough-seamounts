@info "Starting to plot video..."

#+++ Setup Makie backend based on environment
is_headless() = !haskey(ENV, "DISPLAY") || isempty(ENV["DISPLAY"])
has_opengl() = try !isempty(read(`glxinfo -B`, String)) catch; false end
if is_headless() || !has_opengl()
    @info "Either headless environment, or environment without openGL detected. Loading CairoMakie"
    using CairoMakie
    get!(ENV, "GKSwstype", "nul")
    Mk = CairoMakie
else
    @info "Interactive environment. Loading GLMakie"
    using GLMakie
    Mk = GLMakie
end
#---

#+++ Load required packages
using Printf: @sprintf
using Oceananigans
using Oceananigans.Units
using Oceananigans: prettytime
import NCDatasets
#---

#+++ Read datasets
# Get main dataset paths
simname_fallback = "balanus_Ro_b0.05_Fr_b0.3_L0.8_dz4_T_adv_spinup12"
fpath_xyii = (@isdefined simulation) ? simulation.output_writers[:nc_xyii].filepath : "data/xyii.$simname_fallback.nc"
fpath_xizi = (@isdefined simulation) ? simulation.output_writers[:nc_xizi].filepath : "data/xizi.$simname_fallback.nc"
#---

#+++ Get parameters
if !((@isdefined params) && (@isdefined simulation))
    # Read metadata from NetCDF file
    NCDatasets.NCDataset(fpath_xyii) do ds
        global params = (; (Symbol(k) => ds.attrib[k] for k in keys(ds.attrib))...)
    end
end
#---

#+++ Setup animation parameters
# Get times from the first FieldTimeSeries
times = FieldTimeSeries(fpath_xyii, "u", architecture=CPU()).times
n_times = length(times)
max_frames = 200
frame_step = max(1, floor(Int, n_times / max_frames))
frames = 1:frame_step:n_times

@info "Animation setup: $n_times time steps → $(length(frames)) frames (step = $frame_step)"
#---

#+++ Define plotting parameters
# Color ranges for each variable
color_ranges = Dict(
    :u  => (range=(-params.U∞, +params.U∞) .* 1.2, colormap=:balance),
    :v  => (range=(-params.U∞, +params.U∞) .* 1.2, colormap=:balance),
    :PV => (range=params.N²∞ * abs(params.f₀) * [-5, +5], colormap=:seismic),
    :εₖ => (range=(1e-10, 1e-7), colormap=:inferno, colorscale=log10),
    :Ro => (range=(-4, +4), colormap=:balance)
)

# Layout parameters
layout_params = (
    title_height = 8,
    panel_width = 300,
    cbar_height = 8,
    column_gap = 20,
    row_gap = -50,
    title_row_gap = -50
)
#---

#+++ Create figure and setup
fig = Figure(figure_padding = (10, 30, 10, 10))
n = Observable(1)

# Create title in two lines within one row
title = @lift "α = $(@sprintf "%.2g" params.α), Frₕ = $(@sprintf "%.2g" params.Fr_b), Roₕ = $(@sprintf "%.2g" params.Ro_b), Sᴮᵘ = $(@sprintf "%.2g" params.Slope_Bu), Δz = $(@sprintf "%.2g" params.Δz_min) m;   Time = $(@sprintf "%s" prettytime(times[$n])) = $(@sprintf "%.3g" times[$n]/params.T_adv) adv periods = $(@sprintf "%.3g" times[$n]/params.T_inertial) Inertial periods"

# Create single title row with two lines
fig[1, 1:3] = Label(fig, title, fontsize=18, tellwidth=false, height=layout_params.title_height)

colgap!(fig.layout, layout_params.column_gap)
rowgap!(fig.layout, layout_params.row_gap)

# Configure column widths: xyii plot column wide, xizi plot column wide, colorbar column narrow
colsize!(fig.layout, 1, Auto(1.0))  # xyii plot column (main width)
colsize!(fig.layout, 2, Auto(1.0))  # xizi plot column (main width)
colsize!(fig.layout, 3, Auto(0.6))  # Colorbar column (wider for better spacing)

panel_width = layout_params.panel_width
panel_height = 150  # Fixed height for all panels
#---

#+++ Create axes and plots explicitly
#+++ Row 1
@info "Creating panel: u"
ax_u_xyii = Axis(fig[2, 1]; ylabel="y (m)", width=panel_width, height=panel_height)
u_xyii = FieldTimeSeries(fpath_xyii, "u", architecture=CPU())
u_xyiiₙ = @lift u_xyii[$n]
hidexdecorations!(ax_u_xyii)
hm_u_xyii = heatmap!(ax_u_xyii, u_xyiiₙ; colorrange=color_ranges[:u].range, colormap=color_ranges[:u].colormap)

ax_u_xizi = Axis(fig[2, 2]; ylabel="z (m)", width=panel_width, height=panel_height)
u_xizi = FieldTimeSeries(fpath_xizi, "u", architecture=CPU())
u_xiziₙ = @lift u_xizi[$n]
hidexdecorations!(ax_u_xizi)
hm_u_xizi = heatmap!(ax_u_xizi, u_xiziₙ; colorrange=color_ranges[:u].range, colormap=color_ranges[:u].colormap, interpolate=false)

Colorbar(fig[2, 3], hm_u_xyii; label="u", vertical=true, width=layout_params.cbar_height, height=panel_height, ticklabelsize=12)
#---

#+++ Row 2
@info "Creating panel: PV"
ax_PV_xyii = Axis(fig[3, 1]; ylabel="y (m)", width=panel_width, height=panel_height)
PV_xyii = FieldTimeSeries(fpath_xyii, "PV", architecture=CPU())
PV_xyiiₙ = @lift PV_xyii[$n]
hidexdecorations!(ax_PV_xyii)
hm_PV_xyii = heatmap!(ax_PV_xyii, PV_xyiiₙ; colorrange=color_ranges[:PV].range, colormap=color_ranges[:PV].colormap, interpolate=false)

ax_PV_xizi = Axis(fig[3, 2]; ylabel="z (m)", width=panel_width, height=panel_height)
PV_xizi = FieldTimeSeries(fpath_xizi, "PV", architecture=CPU())
PV_xiziₙ = @lift PV_xizi[$n]
hidexdecorations!(ax_PV_xizi)
hm_PV_xizi = heatmap!(ax_PV_xizi, PV_xiziₙ; colorrange=color_ranges[:PV].range, colormap=color_ranges[:PV].colormap, interpolate=false)

Colorbar(fig[3, 3], hm_PV_xyii; label="PV", vertical=true, width=layout_params.cbar_height, height=panel_height, ticklabelsize=12)
#---

#+++ Row 3
@info "Creating panel: εₖ"
ax_εₖ_xyii = Axis(fig[4, 1]; ylabel="y (m)", width=panel_width, height=panel_height)
εₖ_xyii = FieldTimeSeries(fpath_xyii, "εₖ", architecture=CPU())
εₖ_xyiiₙ = @lift εₖ_xyii[$n]
hidexdecorations!(ax_εₖ_xyii)
hm_εₖ_xyii = heatmap!(ax_εₖ_xyii, εₖ_xyiiₙ; colorrange=color_ranges[:εₖ].range, colormap=color_ranges[:εₖ].colormap, interpolate=false)

ax_εₖ_xizi = Axis(fig[4, 2]; ylabel="z (m)", width=panel_width, height=panel_height)
εₖ_xizi = FieldTimeSeries(fpath_xizi, "εₖ", architecture=CPU())
εₖ_xiziₙ = @lift εₖ_xizi[$n]
hidexdecorations!(ax_εₖ_xizi)
hm_εₖ_xizi = heatmap!(ax_εₖ_xizi, εₖ_xiziₙ; colorrange=color_ranges[:εₖ].range, colormap=color_ranges[:εₖ].colormap, interpolate=false)

Colorbar(fig[4, 3], hm_εₖ_xyii; label="εₖ", vertical=true, width=layout_params.cbar_height, height=panel_height, ticklabelsize=12)
#---

#+++ Row 4
@info "Creating panel: Ro"
ax_Ro_xyii = Axis(fig[5, 1]; xlabel="x (m)", ylabel="y (m)", width=panel_width, height=panel_height)
Ro_xyii = FieldTimeSeries(fpath_xyii, "Ro", architecture=CPU())
Ro_xyiiₙ = @lift Ro_xyii[$n]
hm_Ro_xyii = heatmap!(ax_Ro_xyii, Ro_xyiiₙ; colorrange=color_ranges[:Ro].range, colormap=color_ranges[:Ro].colormap, interpolate=false)

ax_Ro_xizi = Axis(fig[5, 2]; xlabel="x (m)", ylabel="z (m)", width=panel_width, height=panel_height)
Ro_xizi = FieldTimeSeries(fpath_xizi, "Ro", architecture=CPU())
Ro_xiziₙ = @lift Ro_xizi[$n]
hm_Ro_xizi = heatmap!(ax_Ro_xizi, Ro_xiziₙ; colorrange=color_ranges[:Ro].range, colormap=color_ranges[:Ro].colormap, interpolate=false)

Colorbar(fig[5, 3], hm_Ro_xyii; label="Ro", vertical=true, width=layout_params.cbar_height, height=panel_height, ticklabelsize=12)
#---
#---

#+++ Adjust figure and record animation
@info "Recording animation with $(length(frames)) frames"
resize_to_layout!(fig)

Mk.record(fig, "$(@__DIR__)/../anims/$(params.simname).mp4", frames,
          framerate=14, compression=30, px_per_unit=1) do frame
    @info "Frame $frame / $(frames[end])"
    n[] = frame
end

@info "Animation saved successfully!"
#---