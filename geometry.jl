module Geometry

using ImmutableArrays, LightRay

export Sphere, Lens, Retina
export intersect!, advance!

################### Sphere #########################
# this is going to be used by both the lens and the retina, since they are both sphere shaped
type Sphere 
    org::Vector3
    rad::Real
    # Easy to avoid, but why not?
    Sphere(org,rad) = rad < 0. ? error("negative sphere radius") : new(org,rad)
end

# functions useful in ray tracing
# Detect the intersection of the ray with the geometry
# This function finds the intersection point of the light as it starts from its 
# initial point. It's a basic vector-sphere intersection thing. It works. 
# I copied it from the ray_sphere function in 
# https://github.com/JuliaLang/julia/blob/master/test/perf/kernel/raytracer.jl
function intersect!(r::Ray, s::Sphere)
    v = s.org - r.pos
    b = dot(v, r.dir)
    disc = b*b - dot(v, v) + s.rad*s.rad
    if disc >= 0.
        d = sqrt(disc)
        t2 = b + d
        if t2 >= 0.
            t1 = b - d
            r.pos += (t1 > 0. ? t1 : t2)*r.dir
        end
    end
end

########################## Lens #######################
# A simple spherical Lens type, should be more generic in the future
type Lens 
    s::Sphere
    
    # the Refractive Index Gradient (RIG) that describes the distribution of refractive 
    # indices inside the lens. We restrict ourselves to spherically symmetrical 
    # distributions: so at a set distance form the lens center, the refractive index is 
    # equal in all directions (kind of like an onion?).
    rig::Function
end

# this is a constructor function, I learned about it from the Julia manual. 
# Useful for constructing instances of Lens
Lens(o,r,f) = Lens(Sphere(Vector3(o),r),f) 

# redefine this for the Lens
intersect!(r::Ray, l::Lens) = intersect!(r, l.s)

# This is the main function that advances the ray along its path inside the lens. 
# Everything happens here, and in fact too much does. It would be better if we'd brake it 
# apart into logical segments. The algorithm is taken directly from page 74 in the 
# S1-rt.pdf file I got from http://fileadmin.cs.lth.se/cs/Education/EDAN30/lectures/S1-rt.pdf 
# in that website your graphics guy gave you 
# (i.e. http://cs.lth.se/english/course/edan30-photorealistic-computer-graphics/assignments/). 
function advance!(r::Ray,l::Lens, step::Real) 
    # first I find the refractive index at the current position
    n1 = l.rig(norm(r.pos - l.s.org))

    # then I advance the ray one step forward with its original direction
    r.pos += step*r.dir

    # useful
    po = r.pos - l.s.org

    # now I'm in the new position where the ray will refract. So there's a new refractive index here:
    n2 = l.rig(norm(po))

    # here I find out if the ray is on it's way in or way out. 
    # If its distance to the lens center is diminishing then the normal to the refractive 
    # surface is going outwards. But if the ray is leaving the lens the normal to 
    # the surface is pointing inwards.
    costheta = dot(r.dir,po)
    
    # is it going out? If so make out equal to +1, otherwise make it -1
    out = -sign(costheta)
    
    # this is the normal to the refractive surface 
    N = out*unit(po)
    
    # refraction, i.e. changing the direction of the ray, occurs with this new function
    refract!(r,N,n1/n2)
end



##################### Retina #########################
# the retina type, simple for now
type Retina 
    s::Sphere
end

Retina(o,r) = Retina(Sphere(Vector3(o),r))

# redefine this for the retina
intersect!(r::Ray, ret::Retina) = intersect!(r, ret.s)

end
