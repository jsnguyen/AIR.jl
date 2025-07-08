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


for i in 30:40
    imsize = (i,i)
    xs = (1:imsize[2])'
    ys = 1:imsize[1]

    for j in 11:25

        @info "Testing crop size $(j) with image size $(i)"

        randx = 2*rand() - 1
        randy = 2*rand() - 1

        center_x = imsize[2]/ 2 + randx
        center_y = imsize[1]/ 2 + randy

        center_x = pixel_center_coordinates(center_x)
        center_y = pixel_center_coordinates(center_y)

        @info "True center" center_x=center_x center_y=center_y
        a = gaussian_2d.(xs, ys, 1.0, center_x, center_y, 2, 2, 0)

        crop_size = (j,j)
        initial_guess = [1.0, center_x, center_y, 2.0, 2.0, 0.0]
        cropped_a, final_cx, final_cy, oy, ox = fit_and_crop(a, crop_size, initial_guess)

        p1 = heatmap(cropped_a; aspect_ratio=:equal, c=:hot, size=(600,600))
        center_x = pixel_center_coordinates(size(cropped_a,2)/2)
        center_y = pixel_center_coordinates(size(cropped_a,1)/2)
        hline!(p1, [center_y], color=:red, width=4, linestyle=:dash)
        vline!(p1, [center_x], color=:red, width=4, linestyle=:dash)

        p2 = heatmap(a; aspect_ratio=:equal, c=:hot, size=(600,600))
        hline!(p2, [final_cy], color=:blue, width=4, linestyle=:dash)
        vline!(p2, [final_cx], color=:blue, width=4, linestyle=:dash)

        @info "Cropped center" final_cx=final_cx-ox final_cy=final_cy-oy
        @info "Final center" final_cx=final_cx final_cy=final_cy

        plot(p1,p2, layout=(1,2), size=(1200,600))
        savefig("reductions/crop_test/crop_test_$(i)_$(j).png")

    end
end