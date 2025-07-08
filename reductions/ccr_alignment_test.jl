using AstroImages
using ImageTransformations
using ImageFiltering
using Statistics
using ProgressMeter
using CoordinateTransformations
using Plots
using Optim
using ForwardDiff

using AIR

# make template

template_len = 10
imsize = (template_len, template_len)
xs = (1:imsize[2])'
ys = 1:imsize[1]

center_x = imsize[2]/ 2
center_y = imsize[1]/ 2
center_x = pixel_center_coordinates(center_x)
center_y = pixel_center_coordinates(center_y)

template = gaussian_2d.(xs, ys, 1.0, center_x, center_y, 2, 2, 0)

# make image offset by some small amount

randx = 2*rand() - 1
randy = 2*rand() - 1
@info "Random injection" randx=randx randy=randy

img = 30
imsize = (img, img)
xs = (1:imsize[2])'
ys = 1:imsize[1]

center_x = imsize[2]/ 2 + randx
center_y = imsize[1]/ 2 + randy
center_x = pixel_center_coordinates(center_x)
center_y = pixel_center_coordinates(center_y)

@info "True center" center_x=center_x center_y=center_y
image = gaussian_2d.(xs, ys, 1.0, center_x, center_y, 2, 2, 0)

aligned = cross_correlate_align(image, template)
#@info "Final center" final_cx=final_cx final_cy=final_cy

p = heatmap(aligned; aspect_ratio=:equal, c=:hot, size=(600,600))

cy, cx = size(aligned)./2
cy = pixel_center_coordinates(cy)
cx = pixel_center_coordinates(cx)

hline!(p, [cy], color=:red, width=4, linestyle=:dash)
vline!(p, [cx], color=:red, width=4, linestyle=:dash)

p