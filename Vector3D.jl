module Vector3D # this is simpler and faster than the ImmutableArrays

export Vec, +, *, -, /, dot, norm, unitize, myconvert

immutable Vec # a three dimensional vector
    x::Float64
    y::Float64
    z::Float64
end
Vec(xyz::Vector{Float64}) = Vec(xyz...)
# various functions for manipulatin the vector type:
+(a::Vec, b::Vec) = Vec(a.x+b.x, a.y+b.y, a.z+b.z)
*(a::Float64, b::Vec) = Vec(a*b.x, a*b.y, a*b.z)
-(a::Vec, b::Vec) = Vec(a.x-b.x, a.y-b.y, a.z-b.z)
-(a::Vec) = *(-1.,a)
/(a::Vec, b::Float64) = Vec(a.x/b, a.y/b, a.z/b)
# like in linear algebra:
dot(a::Vec, b::Vec) = (a.x*b.x + a.y*b.y + a.z*b.z)
norm(a::Vec) = sqrt(dot(a, a))
unitize(a::Vec) = (a/norm(a))
myconvert(a::Vec) = [a.x,a.y,a.z] # just a helper function to convert this to a matrix
end
