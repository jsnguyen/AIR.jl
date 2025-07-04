using OrderedCollections
using Glob
using AstroImages
using ImageTransformations
using ImageFiltering
using Statistics
using ProgressMeter
using CoordinateTransformations
using LinearAlgebra


using AIR



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
            cy_fit, cx_fit, _, _, _ = fit_gaussian_center_lstsq(frame.data, sigma=5.0)
            
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
            circle_mask = make_circle_mask(size(template_psf), 10)
            template_psf_cored = copy(template_psf)  # Create a copy to modify
            template_psf_cored[circle_mask] .= 0.0

            
            # Save results
            save(joinpath(sequences_folder, "$(key)_cropped_sequence.fits"), framelist_to_cube(cropped_frames))
            save(joinpath(sequences_folder, "$(key)_centered_sequence.fits"), framelist_to_cube(centered_frames))
            save(joinpath(sequences_folder, "$(key)_template_psf.fits"), template_psf)
            save(joinpath(sequences_folder, "$(key)_template_psf_cored.fits"), template_psf_cored)
        end

    end

end