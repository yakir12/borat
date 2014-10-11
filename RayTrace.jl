module RayTrace # main model that calls on all the other modules

using Shapes, Vector3D, RefractiveIndexGradients

export trace! # nice, just one funnction to export

type Ray # the ray type
    pos::Vec # has a position, a location
    dir::Vec # and a direction, which has to be unitized, norm(dir) = 1
end

type Bull # this type is really just a container for the following fields. It can be thought of as a ray in the gradient lens. 
    ray::Ray # the ray, note its origin is now the center of the lens
    N::Vec # the normal to the refractive surface
    dis::Float64 # the distance to the lens center
    n1::Float64 # the refractive index at the previous position
end

# functions:
# initiates the Bull type.
function Bull(ray::Ray,lens::Lens) 
    ray.pos -= lens.c # set ray's origin on the lens' center
    return Bull(ray,Vec(0.,0.,0.),0.,0.)
end

# initiates the array of Rays. It starts them off from directly above the lens at distance L. Their directions are uniformly and randomly distributed across the lens.
function letBlight(lens::Lens,L::Float64)
    if isinf(L) # if the light source is infinitely far away
        # just ditribute the rays uniformly on the disk that the sphere makes:
        sqrtr = sqrt(rand()) # this is taken from http://mathworld.wolfram.com/DiskPointPicking.html
        θ = rand()*2π
        a = Vec([sqrtr*cos(θ),sqrtr*sin(θ),2lens.r]) # initiate the Ray with those starting points at distance 2*rad (could be any distance larger than rad really
        pos = lens.c + a # adjust in case the lens center is not at 0,0,0
        dir = Vec([0.,0.,-1.]) # just point them directly down
    else
        cosα = cos(asin(lens.r/L)) # theta is the angle between the optical axis and the ray that touches the surface of the sphere

        # z and α are randomly sampled from those ranges to generate the uniformly spaced points on the sphere
        z = rand()*(cosα - 1.) - cosα
        θ = rand()*2π

        # the x, y, and z that lie on the sphere's surface
        x = sqrt(1. - z*z)*cos(θ)
        y = sqrt(1. - z*z)*sin(θ)

        # this is the starting position of all the rays -- point source
        a = Vec([0.,0.,L]) 
        pos = lens.c + a
        dir = unitize(Vec([x,y,z]))
    end
    return Ray(pos,dir) # the ray type
end

# This function finds the intersection point of the light as it starts from its initial point. It's a basic vector-sphere intersection thing. It works. I copied it from the ray_sphere function in https://github.com/JuliaLang/julia/blob/master/test/perf/kernel/raytracer.jl
function intersect!(ray::Ray,s::Sphere)
    v = s.c - ray.pos
    b = dot(v, ray.dir)
    disc = b*b - dot(v, v) + s.r*s.r
    if disc >= 0.
        d = sqrt(disc)
        t2 = b + d
        if t2 >= 0.
            t1 = b - d
            ray.pos += (t1 > 0. ? t1 : t2)*ray.dir
        end
    end
end

# the refraction/reflection function, ray is ray, N is the normal to the refractive durface, n is the ratio between the refractive index before and afer the surface.
function refract!(ray::Ray,N::Vec,n::Float64)
    # this is the 'r' in page 74 in that pdf I just mentioned
    a = -dot(ray.dir,N)
    
    # this is the 'c' in page 74
    c = 1. - n*n*(1. - a*a)
    
    # and finally the 't' in page 74
    # added reflection as well... See page 54 (v is -ray.dir)
    t = c>0 ? n*ray.dir + (n*a - sqrt(c))*N : 2* dot(N,-ray.dir) *N + ray.dir

    # unitize the direction vector for the ray
    ray.dir = unitize(t)
end

# a step function. moves the ray forward one ε. Notice that the ray is now packed inside Bull, and is normalized to the lens center
function step!(bull::Bull,lens::Lens)
    bull.n1 = getRefractiveIndex(lens.rig, bull.ray.pos) # first I find the refractive index at the current position

    bull.ray.pos += lens.ε*bull.ray.dir # then I move the ray one ε forward with its original direction

    bull.dis = norm(bull.ray.pos) # get the distance to the lens center from the new location
    bull.N = bull.ray.pos/bull.dis # find the new location's normal 
end

# bend the ray accordingly to the refraction/reflection. 
function bend!(bull::Bull,lens::Lens) 
    n2 = getRefractiveIndex(lens.rig, bull.ray.pos) # new refractive index here

    # here I find out if the ray is on it's way in or way out. If its distance to the lens center is diminishing then the normal to the refractive surface is going outwards. But if the ray is leaving the lens the normal to the surface is pointing inwards.
    costheta = dot(bull.ray.dir,bull.ray.pos)
    
    out = -sign(costheta) # is it going out? If so make out equal to +1, otherwise make it -1
    
    bull.N = out*bull.ray.pos/bull.dis # this is the normal to the refractive surface 
    
    refract!(bull.ray,bull.N,bull.n1/n2) # refraction

end

# main tracing function. p is a matrix where all the ray's locations are saved. lens is lens, retina is retina, L is the distance to the light source, n_medium you know, maxiterATION, diagnose is a Dict that contains some diagnostics.
function trace!(p::Array{Float64,2},lens::Lens,retina::Retina,L::Float64,n_medium::Float64,maxiter::Int,diagnose::Dict)
    ray = letBlight(lens,L) # initiate all the rays

    p[:,1] = myconvert(ray.pos) # populate with the start position

    intersect!(ray,lens) # find the intersection point of the ray with the lens and update ray's position

    p[:,2] = myconvert(ray.pos) # add said intersection to the position matrix

	# made this clearer...
	lensRI = getRefractiveIndex(lens.rig, ray.pos)

    refract!(ray,unitize(ray.pos - lens.c),n_medium/lensRI) # refract at the surface of the lens due to the medium's refractive index (that might be different from that of the lens' periphery)
    bull = Bull(ray,lens) # make a bull instance
    for i = 3:maxiter-1 # go for maximum maxiterATIONS
        step!(bull,lens) # advance one ε
        if bull.dis > lens.r # is the ray exiting the lens? If so, intersect it with the lens (which is therefore closer)
            ray = bull.ray # extract the ray from bull
            ray.pos -= lens.ε*ray.dir-lens.c # correct for the mistake of exiting the lens AND the lens center
            intersect!(ray,lens) # intersect with the lens
            p[:,i] = myconvert(ray.pos) # add position to the position matrix
            diagnose[:Δ] = sqrt(sum((diagnose[:goal]-vec(p[:,i])).^2)) # find distance to the goal

            break # and break the loop, cause we're done with the lens
        end # if the ray didn't exit the lens: 
        p[:,i] = myconvert(bull.ray.pos + lens.c) # save the new position in the positions matrix
        bend!(bull,lens) # and refract 
    end
    
    # again refractive index computation
    lensRI = getRefractiveIndex(lens.rig, ray.pos)
    
    refract!(ray,-unitize(ray.pos - lens.c),lensRI/n_medium) # refract when exiting the lens and entering the surrounding medium
    intersect!(ray,retina) # intersection point with the retina
    i += 1 # just for the position matrix
    p[:,i] = myconvert(ray.pos) # add said position to the position matrix
    diagnose[:len] = i # and set the length of the ray refraction points in the positions matrix
end

end
