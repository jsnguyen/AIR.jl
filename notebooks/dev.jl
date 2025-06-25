### A Pluto.jl notebook ###
# v0.20.9

using Markdown
using InteractiveUtils

# ╔═╡ 65f23f24-da36-438f-84bb-197070abcd84
using Pkg

# ╔═╡ a22e5e38-93a6-446d-9c78-6e0f493643dc
using Statistics

# ╔═╡ b8d05614-a814-499a-bd7e-3fba5c02907b
using AstroImages

# ╔═╡ d8e16029-6472-446d-bbc0-b622ea5c4653
using AIR

# ╔═╡ ebbdc3b6-a3ae-49b4-a398-a5547686c67d
# ╠═╡ disabled = true
#=╠═╡
Pkg.develop(path="/Users/jsn/landing/projects/AIR.jl")
  ╠═╡ =#

# ╔═╡ ac922521-c501-4eec-a8dd-4f6c9215ac06
obslog_filename = "/Users/jsn/landing/projects/AIR.jl/reductions/obslogs/2005-07-27_AS_209.toml"

# ╔═╡ cdab589c-0c40-4ba3-b697-bd35c9ff2820
obslog = load_obslog(obslog_filename)

# ╔═╡ 0cfdf547-93da-47e5-bf01-3e4a94c00575
darks_frames = load_frames(obslog, "darks")

# ╔═╡ cee7c328-6b61-41cb-840e-37b17659eb22
darks_dict = match_keys(darks_frames)

# ╔═╡ 4479ecb7-94b0-45cc-8a4a-95715eb35e67
flats_frames = load_frames(obslog, "flats_lampon")

# ╔═╡ f00bff5b-dc8c-4bef-a637-b59b0b165ce3
flats_frames[68]["FILENAME"]

# ╔═╡ 596a2357-3308-4095-a2d2-e0ca79ad0a48
flats_dict = match_keys(flats_frames)

# ╔═╡ ce498cde-40c8-454a-80b4-92bc081c7e89
for key in keys(darks_dict)
    println("Dark frame key: $(key), count: $(length(darks_dict[key]))")
end

# ╔═╡ 673a8e66-1627-43bd-8b14-f77c9ab1005d
for key in keys(flats_dict)
    println("Flat frame key: $(key), count: $(length(flats_dict[key]))")
end

# ╔═╡ 8743466b-eda6-4db3-9464-8111990c5940
master_darks = Dict{Any, AstroImage}()

# ╔═╡ 36116556-8ef7-4d7a-a823-6889c839407c
for key in keys(darks_dict)
    master_darks[key] = mean(darks_dict[key])
end

# ╔═╡ c3e7718e-2a34-4e4e-bdff-cfa7b0aea25f
master_flats = Dict{Any, AstroImage}()

# ╔═╡ e4f01cfa-f4a0-44f1-afaa-a8cf607e464b
for key in keys(flats_dict)
    master_flats[key] = mean(flats_dict[key])
end

# ╔═╡ 5e99257f-32b0-4161-96b6-0c5df8b00d2f
iv = imview(master_flats[collect(keys(master_flats))[1]]; clims=Percent(99), cmap=:magma, stretch=identity, contrast=1.0, bias=0.5)

# ╔═╡ Cell order:
# ╠═a22e5e38-93a6-446d-9c78-6e0f493643dc
# ╠═b8d05614-a814-499a-bd7e-3fba5c02907b
# ╠═65f23f24-da36-438f-84bb-197070abcd84
# ╠═ebbdc3b6-a3ae-49b4-a398-a5547686c67d
# ╠═d8e16029-6472-446d-bbc0-b622ea5c4653
# ╠═ac922521-c501-4eec-a8dd-4f6c9215ac06
# ╠═cdab589c-0c40-4ba3-b697-bd35c9ff2820
# ╠═0cfdf547-93da-47e5-bf01-3e4a94c00575
# ╠═cee7c328-6b61-41cb-840e-37b17659eb22
# ╠═4479ecb7-94b0-45cc-8a4a-95715eb35e67
# ╠═f00bff5b-dc8c-4bef-a637-b59b0b165ce3
# ╠═596a2357-3308-4095-a2d2-e0ca79ad0a48
# ╠═ce498cde-40c8-454a-80b4-92bc081c7e89
# ╠═673a8e66-1627-43bd-8b14-f77c9ab1005d
# ╠═8743466b-eda6-4db3-9464-8111990c5940
# ╠═36116556-8ef7-4d7a-a823-6889c839407c
# ╠═c3e7718e-2a34-4e4e-bdff-cfa7b0aea25f
# ╠═e4f01cfa-f4a0-44f1-afaa-a8cf607e464b
# ╠═5e99257f-32b0-4161-96b6-0c5df8b00d2f
