### A Pluto.jl notebook ###
# v0.20.13

using Markdown
using InteractiveUtils

# ╔═╡ 591b7532-d974-42b3-bda6-69c2e1b042ae
begin
	using Pkg
	using Plots
	using AstroImages
end

# ╔═╡ ebbdc3b6-a3ae-49b4-a398-a5547686c67d
Pkg.develop(path="/Users/jsn/landing/projects/AIR.jl")

# ╔═╡ d8e16029-6472-446d-bbc0-b622ea5c4653
using AIR

# ╔═╡ 964d7aa9-ac1d-4d30-8e0c-61a4faf24027
image_path = "/Users/jsn/landing/projects/AIR.jl/data/2002-06-16/sequences/as209_median.fits"

# ╔═╡ 8ae98f0c-14c6-4970-bac4-5f57ef8fe27b
begin
	frame = load(image_path)
	outer_circle_mask = make_circle_mask(size(frame), 230)
	inner_circle_mask = make_circle_mask(size(frame), 60)
	frame[.!outer_circle_mask] .= NaN
	frame[inner_circle_mask] .= NaN
end

# ╔═╡ 81f8369e-6ea6-4c60-8b50-897dd0cc9ea9
begin

	function circle(x, y, r=1; n=30)
    θ = 0:2π/n:2π
    Plots.Shape(r*sin.(θ) .+ x, r*cos.(θ) .+ y)
	end
	
	heatmap(range(-2.5, 2.5, length=500),
			range(-2.5, 2.5, length=500),
			frame.data,
			aspect_ratio=:equal,
			xlims=(-2.5,2.5),
			ylims=(-2.5,2.5),
			title="AS 209",
			xlabel="[arcsec]",
			ylabel="[arcsec]",
			cmap=:matter,
			colorbar=true,
			clim=(-30, 30),
		    colorbar_title="[e-/s]",
		   	legend=false,
		   	size=(400,400))
	
	plot!(circle(-0.88, -1.21, 0.25), seriestype=[:shape], c=:white, fillalpha=0, width=4, linestyle=:dash, linecolor=:white)

end

# ╔═╡ Cell order:
# ╠═591b7532-d974-42b3-bda6-69c2e1b042ae
# ╠═ebbdc3b6-a3ae-49b4-a398-a5547686c67d
# ╠═d8e16029-6472-446d-bbc0-b622ea5c4653
# ╠═964d7aa9-ac1d-4d30-8e0c-61a4faf24027
# ╠═8ae98f0c-14c6-4970-bac4-5f57ef8fe27b
# ╠═81f8369e-6ea6-4c60-8b50-897dd0cc9ea9
