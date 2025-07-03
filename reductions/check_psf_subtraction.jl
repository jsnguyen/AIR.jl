using AIR
using AstroImages
using ImageTransformations
using CoordinateTransformations
using Plots
using Statistics

# Load the results from PSF subtraction
sequences_folder = "data/2002-06-16/sequences"

# Load the original and shifted PSFs and residual
as209_median = load(joinpath(sequences_folder, "as209_median.fits"))
hbc650_median = load(joinpath(sequences_folder, "hbc650_median.fits"))
hbc650_shifted = load(joinpath(sequences_folder, "hbc650_median_shifted.fits"))
hbc650_scaled_shifted = load(joinpath(sequences_folder, "hbc650_median_scaled_shifted.fits"))
residual = load(joinpath(sequences_folder, "psf_subtraction_residual.fits"))

# Clean any NaN values for statistics and plotting
hbc650_shifted_clean = replace(hbc650_shifted, NaN => 0.0, Inf => 0.0, -Inf => 0.0)
hbc650_scaled_shifted_clean = replace(hbc650_scaled_shifted, NaN => 0.0, Inf => 0.0, -Inf => 0.0)

println("Data loaded successfully:")
println("AS209 median size: $(size(as209_median))")
println("HBC650 median size: $(size(hbc650_median))")
println("HBC650 shifted size: $(size(hbc650_shifted))")
println("HBC650 scaled+shifted size: $(size(hbc650_scaled_shifted))")
println("Residual size: $(size(residual))")

# Print statistics
println("\nStatistics:")
println("AS209 median: min=$(round(minimum(as209_median), digits=3)), max=$(round(maximum(as209_median), digits=3)), std=$(round(std(as209_median), digits=3))")
println("HBC650 median: min=$(round(minimum(hbc650_median), digits=3)), max=$(round(maximum(hbc650_median), digits=3)), std=$(round(std(hbc650_median), digits=3))")
println("HBC650 shifted: min=$(round(minimum(hbc650_shifted_clean), digits=3)), max=$(round(maximum(hbc650_shifted_clean), digits=3)), std=$(round(std(hbc650_shifted_clean), digits=3))")
println("HBC650 scaled+shifted: min=$(round(minimum(hbc650_scaled_shifted_clean), digits=3)), max=$(round(maximum(hbc650_scaled_shifted_clean), digits=3)), std=$(round(std(hbc650_scaled_shifted_clean), digits=3))")
println("Residual: min=$(round(minimum(residual), digits=3)), max=$(round(maximum(residual), digits=3)), std=$(round(std(residual), digits=3))")

# Create plots showing the PSF subtraction
center = size(as209_median) .÷ 2
crop_size = 100
y_range = (center[1]-crop_size÷2):(center[1]+crop_size÷2)
x_range = (center[2]-crop_size÷2):(center[2]+crop_size÷2)

# Plot the original PSFs and results (now with 6 panels)
p1 = heatmap(as209_median[y_range, x_range], title="AS209 (Target)", aspect_ratio=:equal, c=:hot)
p2 = heatmap(hbc650_median[y_range, x_range], title="HBC650 (Reference)", aspect_ratio=:equal, c=:hot)  
p3 = heatmap(hbc650_shifted_clean[y_range, x_range], title="HBC650 Shifted", aspect_ratio=:equal, c=:hot)
p4 = heatmap(hbc650_scaled_shifted_clean[y_range, x_range], title="HBC650 Scaled+Shifted", aspect_ratio=:equal, c=:hot)
p5 = heatmap(residual[y_range, x_range], title="Residual", aspect_ratio=:equal, c=:RdBu)

plot_combined = plot(p1, p2, p3, p4, p5, layout=(2,3), size=(1200,800))
savefig(plot_combined, "data/2002-06-16/plots/psf_subtraction_comparison.png")

println("\nPSF subtraction comparison plot saved to: data/2002-06-16/plots/psf_subtraction_comparison.png")

# Create a cross-section plot through the center
center_row = center[1]
cross_section_as209 = as209_median[center_row, x_range]
cross_section_hbc650 = hbc650_median[center_row, x_range] 
cross_section_shifted = hbc650_shifted_clean[center_row, x_range]
cross_section_scaled_shifted = hbc650_scaled_shifted_clean[center_row, x_range]
cross_section_residual = residual[center_row, x_range]

x_coords = x_range .- center[2]  # Center coordinates around 0

p_cross = plot(x_coords, cross_section_as209, label="AS209", linewidth=2)
plot!(x_coords, cross_section_hbc650, label="HBC650 Original", linewidth=2)
plot!(x_coords, cross_section_shifted, label="HBC650 Shifted", linewidth=2)
plot!(x_coords, cross_section_scaled_shifted, label="HBC650 Scaled+Shifted", linewidth=2)
plot!(x_coords, cross_section_residual, label="Residual", linewidth=2)
xlabel!("Pixel offset from center")
ylabel!("Intensity")
title!("Cross-section through PSF centers")

savefig(p_cross, "data/2002-06-16/plots/psf_cross_section.png")
println("Cross-section plot saved to: data/2002-06-16/plots/psf_cross_section.png")
