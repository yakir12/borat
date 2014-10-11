module Shapes # like your Geometry module

export Sphere, Lens, Retina

using Vector3D, RefractiveIndexGradients

abstract Sphere

type Lens <: Sphere # the lens type
    c::Vec # center
    r::Float64 # radius
    rig::RadialRIG
    ε::Float64 # the step
end
Lens(a::Vector{Float64},b::Float64,c::Function,d::Float64) = Lens(Vec(a),b,RadialRIG(c,Vec(a)),d) # this is a constructor function

type Retina <: Sphere # the retina type
    c::Vec # retina center
    r::Float64 # radius
end
Retina(a::Vector{Float64},b::Float64) = Retina(Vec(a),b)

end
