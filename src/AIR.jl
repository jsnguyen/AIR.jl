module AIR

using Logging
using LoggingExtras

using Glob
using TOML
using OrderedCollections
using Printf
using Dates

using AstroImages
using SkyCoords
using AstroAngles
using Plots
using Statistics
using ImageFiltering
using LsqFit
using CoordinateTransformations
import CoordinateTransformations: recenter
using Rotations
using AstroLib
using Optim
using ForwardDiff
using LACosmic
using ImageTransformations

include("angles.jl")
include("constants.jl")
include("fitting.jl")
include("io.jl")
include("logging.jl")
include("NIRC2.jl")
include("reduction.jl")
include("utils.jl")

# angles.jl
export calculate_north_angle, par_angle, rotate_image_center

# constants.jl
export observatory_lat, observatory_lon

# fitting.jl
export optimal_subtract_target, fit_generic_kernel, fit_2d_gaussian, gaussian_2d, gaussian_2d_rotated, fit_and_crop, cross_correlation_center, cross_correlate_align

# io.jl
export ObslogPaths, Obslog, load_obslog, load_frames, load_master, load_masks, load_sequences, load_rejects, write_toml

# logging.jl
export autolog, @logname, @autolog

# NIRC2.jl
export _mask_dir, NIRC2_plate_scale, NIRC2_bad_pixel_mask, get_NIRC2_gain, get_NIRC2_readnoise

# reduction.jl
export make_masters, find_matching_master, find_closest_flat, find_closest_dark

# utils.jl
export small_angle_distance, deg2arcsec, arcsec2deg, pixel_center_coordinates, cartesian_coordinates, remove_nan!, all_header_keywords_match, pretty_print_toml, framelist_to_cube, match_keys, make_and_clear, make_circle_mask, make_annulus_mask, make_sigma_clip_mask, crop, subpixel_crop, measure_background, make_frametable, argquantile

end
