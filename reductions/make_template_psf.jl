using OrderedCollections
using Glob
using AstroImages
using ImageTransformations
using ImageFiltering
using Statistics
using ProgressMeter
using CoordinateTransformations
using LinearAlgebra
using Optim
using ForwardDiff

using AIR

"""
    fit_gaussian_center_lstsq(data; sigma=5.0)

Fit a 2D Gaussian using least squares to find the PSF center with sub-pixel precision.
Evaluates the Gaussian continuously on the pixel grid.
"""
function fit_gaussian_center_lstsq(data; sigma=5.0)
    rows, cols = size(data)
    
    # Initial guess: brightest pixel (but we'll optimize from here)
    _, max_idx = findmax(data)
    initial_cy, initial_cx = Float64.(Tuple(max_idx))
    
    # Set up the optimization problem
    function objective(center_params)
        cy, cx = center_params[1], center_params[2]
        
        # Create design matrix for current center position
        n_points = rows * cols
        A = zeros(n_points, 2)  # [gaussian_terms, background_terms]
        b = zeros(n_points)     # pixel values
        
        idx = 1
        for i in 1:rows, j in 1:cols
            # Gaussian term: exp(-0.5 * ((i-cy)²+(j-cx)²) / σ²)
            r_squared = (i - cy)^2 + (j - cx)^2
            gaussian_term = exp(-0.5 * r_squared / sigma^2)
            
            A[idx, 1] = gaussian_term  # amplitude coefficient
            A[idx, 2] = 1.0           # background coefficient
            b[idx] = data[i, j]       # pixel value
            idx += 1
        end
        
        # Solve least squares for amplitude and background
        try
            params = A \ b
            amplitude, background = params[1], params[2]
            
            # Calculate and return residual (what we want to minimize)
            predicted = A * params
            residual = sum((b - predicted).^2)
            
            # Add penalty if amplitude is negative
            if amplitude <= 0
                residual += 1e10
            end
            
            return residual
        catch
            # Return large value if least squares fails
            return 1e10
        end
    end
    
    # Optimize to find best center
    try
        result = optimize(objective, 
                         [initial_cy, initial_cx],
                         NelderMead(),
                         Optim.Options(iterations = 100))
        
        optimized_center = Optim.minimizer(result)
        return optimized_center[1], optimized_center[2]  # cy, cx
    catch
        # Fallback to initial guess if optimization fails
        @warn "Gaussian optimization failed, using brightest pixel"
        return initial_cy, initial_cx
    end
end

autolog("$(@__FILE__).log") do

    sequence_obslog_folder = "reductions/obslogs"
    sequence_obslog_path = joinpath(sequence_obslog_folder, "2002-06-16_sequences.toml")

    @info "Loading sequence_obslog from" sequence_obslog_path
    sequence_obslog = load_obslog(sequence_obslog_path)
    sequences = load_sequences(sequence_obslog)
    sequences_folder = joinpath(sequence_obslog["data_folder"], "sequences")

    unsaturated_sequences = ["as209_1", "hbc650_1", "hbc630_1"]

    coarse_size = 320
    fine_size = 256

    for key in keys(sequences)
        if !(key in unsaturated_sequences)
            @info "Skipping saturated sequence $key"
            continue
        end

        cropped_frames = AstroImage[]
        for frame in sequences[key]
            _, coarse_center = findmax(frame.data)
            cropped = crop(frame, (coarse_size, coarse_size), center=Tuple(coarse_center))
            push!(cropped_frames, cropped)
        end

        # Center each frame individually using Gaussian fitting
        centered_frames = AstroImage[]
        @info "Centering frames using Gaussian fitting..."
        
        for (i, frame) in enumerate(cropped_frames)
            # Fit Gaussian to find PSF center
            cy_fit, cx_fit = fit_gaussian_center_lstsq(frame.data, sigma=5.0)
            
            # Calculate shift to center the PSF
            rows, cols = size(frame.data)
            target_cy, target_cx = rows/2, cols/2
            shift = (cy_fit-target_cy, cx_fit-target_cx)
            
            # Apply centering shift
            centered_data = warp(frame.data,
                                Translation(shift...),
                                axes(frame.data),
                                fill=0.0)
            
            # Crop to final size
            centered_cropped = crop(AstroImage(centered_data, frame.header), 
                                    (fine_size, fine_size))
            push!(centered_frames, centered_cropped)
            
            @info "Frame $i" gaussian_center=(cy_fit, cx_fit) shift=shift
                
        end
        
        @info "Successfully centered $(length(centered_frames)) frames"
        
        # Create template PSF from centered frames
        if length(centered_frames) > 0
            template_stack = framelist_to_cube(centered_frames)
            template_psf = median(template_stack, dims=3) |> x->dropdims(x, dims=3)
            
            # Save results
            save(joinpath(sequences_folder, "$(key)_cropped_sequence.fits"), framelist_to_cube(cropped_frames))
            save(joinpath(sequences_folder, "$(key)_centered_sequence.fits"), framelist_to_cube(centered_frames))
            save(joinpath(sequences_folder, "$(key)_template_psf.fits"), template_psf)
        end

    end

end