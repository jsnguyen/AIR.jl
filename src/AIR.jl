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

export small_angle_distance, pretty_print_toml, deg2arcsec, arcsec2deg, all_header_keywords_match, load_obslog, load_frames, framelist_to_cube, match_keys, crop, subpixel_crop, make_sigma_clip_mask, make_masters, autolog, @logname, @autolog, NIRC2_bad_pixel_mask, find_matching_master, find_closest_flat, find_closest_dark, get_NIRC2_gain, make_and_clear, make_circle_mask, calculate_north_angle, par_angle, load_sequences, fit_gaussian_center_variable_sigma, NIRC2_plate_scale, fit_generic_kernel, fit_2d_gaussian, gaussian_2d, fit_and_crop, pixel_center_coordinates, cartesian_coordinates, cross_correlate_align, measure_background, rotate_image_center, remove_nan!, optimal_subtract_target, load_master, load_masks
end
