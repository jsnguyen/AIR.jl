using Printf
using Statistics
using AstroImages
using Photometry
using Plots

using AIR
import AIR.crop 

#psf = load("AS_209_data/2002-06-16/sequences/AS_209_1_template_psf.fits")
#psf = load("AS_209_data/2002-08-21/sequences/AS_209_2_template_psf_big.fits")
psf = load("AS_209_data/2002-08-21/sequences/AS_209_2_template_psf.fits")

aperture = CircularAperture(size(psf, 2)/2 + 0.5, size(psf, 1)/2 + 0.5, size(psf, 1)/2)
circle_photometry = photometry(aperture, psf)
total = circle_photometry.aperture_sum

sums = []
range = 20:2:size(psf, 1)/2
for i in range
    aperture = CircularAperture(size(psf, 2)/2 + 0.5, size(psf, 1)/2 + 0.5, i)
    circle_photometry = photometry(aperture, psf)
    push!(sums, circle_photometry.aperture_sum)
end

println(size(psf))

for (r,s) in zip(range, sums)
    @printf("%3d px: %8.3e (%.4f)\n", r, s, s/total)
end
plot(range, sums; title="$(size(psf,1)/2)", ylims=(0, 1.2e6))