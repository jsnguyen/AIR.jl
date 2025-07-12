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

# ╔═╡ 109b615c-1151-4a83-9ed4-162006796dc8
dates = ["2002-06-16", "2002-08-02", "2002-08-21","2005-07-27"]

# ╔═╡ 170c8ca7-a1bf-48e0-9e19-ea6b53d25bd2
function circle(x, y, r=1; n=30)
	θ = 0:2π/n:2π
	Plots.Shape(r*sin.(θ) .+ x, r*cos.(θ) .+ y)
end

# ╔═╡ 81f8369e-6ea6-4c60-8b50-897dd0cc9ea9
function make_plot(img,date)
	height, width = size(img)
	println(height, " ", width, " ", size(img))
	scale = 100
	height /= scale
	width /= scale
	px, py = width/2, height/2
	xlims = (-2,2)
	ylims = (-2,2)
	heatmap(range(-px, px, length=size(img,2)),
			range(-py, py, length=size(img,1)),
			img,
			aspect_ratio=:equal,
			xlims=xlims,
			ylims=ylims,
			title="$(date) AS 209",
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

# ╔═╡ 96fb5436-e812-4f5d-856c-a9c2ce4ba1e8
plots = []

# ╔═╡ 8ae98f0c-14c6-4970-bac4-5f57ef8fe27b
for date in dates
	image_path = "/Users/jsn/landing/projects/AIR.jl/data/$(date)/sequences/as209_median.fits"
	frame = load(image_path)
	outer_circle_mask = make_circle_mask(size(frame), 200 )
	inner_circle_mask = make_circle_mask(size(frame), 60)
	frame[.!outer_circle_mask] .= NaN
	frame[inner_circle_mask] .= NaN
	p = make_plot(frame.data, date)
	push!(plots,p)
	savefig(p,"$(date)_plot.png")
end

# ╔═╡ 3f81bab7-5f80-42b2-9e4b-24583559f5f6
plots

# ╔═╡ 6096a7e2-e509-4e6a-8b02-2f5ec43e6f1f
for date in dates
	image_path = "/Users/jsn/landing/projects/AIR.jl/data/$(date)/sequences/as209_cropped.fits"
	frame = load(image_path)
	height, width = size(frame)
	println(height, " ", width, " ", size(frame))
	scale = 100
	height /= scale
	width /= scale
	px, py = width/2, height/2
	p = heatmap(
			frame,
			aspect_ratio=:equal,
			title="$(date) AS 209",
			xlabel="[arcsec]",
			ylabel="[arcsec]",
			cmap=:matter,
			colorbar=true,
			clim=(-30, 30),
		    colorbar_title="[e-/s]",
		   	legend=false,
		   	size=(400,400))
	savefig(p,"$(date)_gauss.png")
end

# ╔═╡ Cell order:
# ╠═591b7532-d974-42b3-bda6-69c2e1b042ae
# ╠═ebbdc3b6-a3ae-49b4-a398-a5547686c67d
# ╠═d8e16029-6472-446d-bbc0-b622ea5c4653
# ╠═109b615c-1151-4a83-9ed4-162006796dc8
# ╠═170c8ca7-a1bf-48e0-9e19-ea6b53d25bd2
# ╠═81f8369e-6ea6-4c60-8b50-897dd0cc9ea9
# ╠═96fb5436-e812-4f5d-856c-a9c2ce4ba1e8
# ╠═8ae98f0c-14c6-4970-bac4-5f57ef8fe27b
# ╠═3f81bab7-5f80-42b2-9e4b-24583559f5f6
# ╠═6096a7e2-e509-4e6a-8b02-2f5ec43e6f1f
