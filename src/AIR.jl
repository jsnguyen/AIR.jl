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
using AstroLib
using Optim
using ForwardDiff
using LACosmic
using ImageTransformations

include("angles.jl")
include("constants.jl")
include("fitting.jl")
include("images.jl")
include("io.jl")
include("logging.jl")
include("NIRC2.jl")
include("reduction.jl")
include("utils.jl")

export small_angle_distance, pretty_print_toml, deg2arcsec, arcsec2deg, all_header_keywords_match, load_obslog, load_frames, framelist_to_cube, match_keys, crop, make_sigma_clip_mask, make_masters, autolog, NIRC2_bad_pixel_mask, find_matching_master, find_closest_flat, find_closest_dark, get_NIRC2_gain, make_and_clear, gaussian2d_fit, gaussian2d_fixedwidth_fit, make_circle_mask, calculate_north_angle, par_angle, load_sequences, fit_gaussian_center_lstsq, fit_gaussian_center_variable_sigma, subtract_psf_with_shift, fit_gaussian_center_variable_sigma, NIRC2_plate_scale



end
