module Shapes # like your Geometry module

export Sphere, Lens, Retina

using Vector3D

abstract Sphere

type Lens <: Sphere # the lens type
    c::Vec # center
    r::Float64 # radius
    ri::Function # the Refractive Index Gradient (RIG) that describes the distribution of refractive indices inside the lens. We restrict ourselves to spherically symmetrical distributions: so at a set distance form the lens center, the refractive index is equal in all directions (kind of like an onion?).
    ε::Float64 # the step
end
Lens(a::Vector{Float64},b::Float64,c::Function,d::Float64) = Lens(Vec(a),b,c,d) # this is a constructor function

type Retina <: Sphere # the retina type
    c::Vec # retina center
    r::Float64 # radius
end
Retina(a::Vector{Float64},b::Float64) = Retina(Vec(a),b)

end
