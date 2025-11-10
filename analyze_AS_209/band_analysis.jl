using Interpolations
using Statistics
using Trapz

function load_magnitude_data(filename)
    wl = Float64[]
    flux_density = Float64[]
    open(filename, "r") do file
        for line in eachline(file)
            if startswith(line, "#") || isempty(strip(line))
                continue
            end
            parts = split(strip(line))
            push!(wl, parse(Float64, parts[1]))
            push!(flux_density, parse(Float64, parts[2]))
        end
    end
    return wl, flux_density
end


"""
    l2_on_uniform_grid(x1,y1,x2,y2; N=2000)

Approximate ∫(f-g)^2 dx by resampling both curves to a uniform grid
over their overlap using linear interpolation.
"""
function l2_on_uniform_grid(x1,y1,x2,y2; N=2000)
    xlo = max(x1[1], x2[1]); xhi = min(x1[end], x2[end])
    xlo < xhi || return 0.0
    xu = range(xlo, xhi, length=N)
    itp1 = LinearInterpolation(x1, y1, extrapolation_bc=Flat())
    itp2 = LinearInterpolation(x2, y2, extrapolation_bc=Flat())
    d = itp1.(xu) .- itp2.(xu)
    # trapezoidal rule
    dx = step(xu)
    return dx * (sum(d.^2) - 0.5*(first(d)^2 + last(d)^2))
end

"""
    area_nonoverlap(x1, y1, x2, y2; N=4000, domain=:overlap, zero_outside=false)

Return ∫ |f(x) - g(x)| dx using linear interpolation and a uniform grid.

Arguments:
- x1,y1 : first curve (x strictly increasing)
- x2,y2 : second curve (x strictly increasing)

Keywords:
- N            : number of grid points
- domain       : :overlap (default) or :union
- zero_outside : if true and domain=:union, treat curve values outside their
                 native x-range as 0. If false, they are held flat at edges.

Also returns a named tuple with `total`, `area1_only`, `area2_only`.
"""
function area_nonoverlap(x1, y1, x2, y2; N=4000, domain=:overlap, zero_outside=false)
    @assert length(x1) == length(y1) && length(x2) == length(y2) "x/y length mismatch"
    @assert issorted(x1) && issorted(x2) "x arrays must be strictly increasing"

    xlo = domain === :overlap ? max(first(x1), first(x2)) : min(first(x1), first(x2))
    xhi = domain === :overlap ? min(last(x1),  last(x2))  : max(last(x1),  last(x2))
    xlo < xhi || return (total=0.0, area1_only=0.0, area2_only=0.0)

    # Linear gridded interpolants + flat extrapolation
    itp1 = extrapolate(interpolate((x1,), y1, Gridded(Linear())), Flat())
    itp2 = extrapolate(interpolate((x2,), y2, Gridded(Linear())), Flat())

    # Optional zeroing outside native domain (only meaningful over :union)
    f1 = zero_outside && domain === :union ?
         (x -> (x < x1[1] || x > x1[end]) ? 0.0 : itp1(x)) :
         (x -> itp1(x))
    f2 = zero_outside && domain === :union ?
         (x -> (x < x2[1] || x > x2[end]) ? 0.0 : itp2(x)) :
         (x -> itp2(x))

    # Uniform grid and trapezoid rule
    xu = range(xlo, xhi; length=N)
    f  = f1.(xu)
    g  = f2.(xu)
    df = f .- g
    d  = abs.(df)

    dx   = (xhi - xlo) / (N - 1)
    trap(z) = dx * (sum(z) - 0.5*(first(z) + last(z)))

    total      = trap(d)
    area1_only = trap(max.(df, 0))   # where f > g
    area2_only = trap(max.(-df, 0))  # where g > f
    return (total=total, area1_only=area1_only, area2_only=area2_only)
end

@info "Integrated Bandpasses" integrated_2MASS=integrated_2MASS integrated_NIRC2=integrated_NIRC2

res = l2_on_uniform_grid(wl_2MASS_Ks, flux_density_2MASS_Ks, wl_NIRC2_Kp, flux_density_NIRC2_Kp)

@info "L2 Difference between 2MASS Ks and NIRC2 Kp bandpasses" L2_difference=res

area = area_nonoverlap(wl_2MASS_Ks, flux_density_2MASS_Ks, wl_NIRC2_Kp, flux_density_NIRC2_Kp; domain=:union, zero_outside=true)

@info "Area Non-overlap between 2MASS Ks and NIRC2 Kp bandpasses" area_nonoverlap=area