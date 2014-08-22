#=useful source: https://github.com/JuliaLang/julia/blob/master/test/perf/kernel/raytracer.jl=#
using ImmutableArrays # useful library for fast matrix manipulations

type Ray # the ray type
    orig::Vector3 # has an origin, or location
    dir::Vector3 # and a direction, which has to be unitized
end

const step = 1e-6 # step size

function initiate(o::Array{Float64},d::Array{Float64})
    # just a utility function to declare a ray
    r = Ray(Vector3(o),unit(Vector3(d)))
end

function advance!(r::Ray)
    # simple but useful
    r.orig += step*r.dir
end

n = int(1e5) # let's go for n steps
r1 = initiate([0.,0.,0.],[0.,1.,1.]) # start a ray at 0,0,0 and give it some direction
for i = 1:n # trace!
    advance!(r1)
end
