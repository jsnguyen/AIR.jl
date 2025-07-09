function make_masters(frames, keylist; n_sigma::Float64=6.0, median_size::Int64=5, method::Function=median)

    frame_dict = match_keys(frames, keylist)

    for key in keys(frame_dict)
        @printf "key -> %s, count -> %d\n " key length(frame_dict[key])
    end

    master_frames = Dict{Any, AstroImage}()
    master_frames_masks = Dict{Any, BitMatrix}()
    for key in keys(frame_dict)

        # applies method (mean, median) to frame list
        stack = framelist_to_cube(frame_dict[key])
        method_stack = method(stack, dims=3) |> x -> dropdims(x, dims=3)
        mf = AstroImage(method_stack)

        # take each frame, make a sigma clip mask, and then combine them or-wise
        masks = [make_sigma_clip_mask(frame.data, n_sigma) for frame in frame_dict[key]]
        sigma_clip_mask = reduce(.|, masks)

        # crop the original bad pixel mask to the size of the median frame
        bad_pixel_mask = copy(NIRC2_bad_pixel_mask)
        if size(bad_pixel_mask) != size(mf)
            bad_pixel_mask, _, _ = crop(NIRC2_bad_pixel_mask, size(mf))
        end

        # combine masks
        mask = bad_pixel_mask .| sigma_clip_mask

        # make the median frame to do pixel replacement
        # not super efficient, but it works
        median_mf = mapwindow(median, mf.data, (median_size, median_size))

        # finally, assign the median values to the masked pixels
        mf.data[mask] .= median_mf[mask]

        # repopulate the header with the key values so we can find the keys later
        for (i,k) in enumerate(keylist)
            mf[k] = key[i]
        end

        mf["NFRAMES"] = length(frame_dict[key])
        mf["NSIGMASK"] = n_sigma
        mf["MEDSIZE"] =  median_size
        mf["NPIXMASK"] = sum(mask)
        mf["MAMEDIAN"] = median(mf.data[.!mask])
        mf["MAMEAN"] = mean(mf.data[.!mask])
        mf["MASTD"] = std(mf.data[.!mask])

        master_frames[key] = mf
        master_frames_masks[key] = mask
    end

    return master_frames, master_frames_masks

end

function find_matching_master(frame, masters, keylist)

    matches = [all_header_keywords_match(frame, m, keylist) for m in masters]

    if any(matches)
        ind = findfirst(matches)
        matched = masters[ind]
    else
        return nothing
    end

    return matched

end



"""
Guarantees finding a flat frame that matches at least the FILTER. Frame should be cropped as well.
"""
function find_closest_flat(frame, master_flats, flats_keylist=["FILTER"])
    matched_flat = find_matching_master(frame, master_flats, flats_keylist)

    if matched_flat !== nothing
        if (size(matched_flat) != size(frame))
            if image_is_larger(matched_flat, frame)
                cropped_flat,_ ,_ = crop(matched_flat.data, size(frame))
                matched_flat = AstroImage(cropped_flat, matched_flat.header)
            else
                @warn "Flat frame is smaller than the target frame, not cropping: $(frame["FILENAME"])"
                return nothing
            end
        end
    else
        @warn "No matching flat found for $(frame["FILENAME"])"
    end

    return matched_flat

end

"""
Guarantees finding a dark frame that matches at least the ITIME. Frame should be cropped as well.

# Arguments
- `frame`: The frame for which we want to find a matching dark.
- `master_darks`: A list of master dark frames.
- `ranked_darks_keylist`: A list of keylists in order of preference for matching dark frames. Each keylist is a list of header keywords that should match. At minimum we want to match the ITIME.
"""
function find_closest_dark(frame, master_darks, ranked_darks_keylist = [["NAXIS1", "NAXIS2", "ITIME", "COADDS"], ["NAXIS1", "NAXIS2", "ITIME"], ["ITIME"]])

    matched_dark = nothing
    for (i,keylist) in enumerate(ranked_darks_keylist)
        matched_dark = find_matching_master(frame, master_darks, keylist)
        if matched_dark !== nothing

            # for the 2nd and 3rd case, rescale by coadds
            if (i==2) || (i==3)
                @warn "Rescaling dark frame by COADDS $(matched_dark["COADDS"]) -> $(frame["COADDS"])"
                matched_dark = matched_dark ./ matched_dark["COADDS"] .* frame["COADDS"]
            end

            if i == 3
                if  (size(matched_dark) != size(frame))
                    if image_is_larger(matched_dark, frame)
                        cropped_dark, _, _ = crop(matched_dark.data, size(frame))
                        matched_dark= AstroImage(cropped_dark, matched_dark.header)
                    else
                        @warn "Dark frame is smaller than the target frame, not cropping: $(frame["FILENAME"])"
                        return nothing
                    end
                end
            end

            break
        end
    end

    if matched_dark === nothing
        @warn "No matching dark found"
    end

    return matched_dark

end