#=module Borat=#

#=export Vec, +, -, *, /, dot, norm, unitize, Ray, Sphere, Lens, Retina, initiate, intersect!, refract!,advance!, poke=#
#=first we define a few useful types. Since Julia is a highly (but not strictly) typed language, this is an awesome feature that speeds up the code and makes it clean and easy to follow=#

using PyPlot

immutable Vec{Float64}
    x::Float64
    y::Float64
    z::Float64
end
Vec(xyz::Vector{Float64}) = Vec(xyz...)

+(a::Vec, b::Vec) = Vec(a.x+b.x, a.y+b.y, a.z+b.z)
*(a::Float64, b::Vec) = Vec(a*b.x, a*b.y, a*b.z)
-(a::Vec, b::Vec) = Vec(a.x-b.x, a.y-b.y, a.z-b.z)
-(a::Vec) = *(-1.,a)
/(a::Vec, b::Float64) = Vec(a.x/b, a.y/b, a.z/b)
dot(a::Vec, b::Vec) = (a.x*b.x + a.y*b.y + a.z*b.z)
norm(a::Vec) = sqrt(dot(a, a))
unitize(a::Vec) = (a/norm(a))


type Ray # the ray type
    pos::Vec # has a position, a location
    dir::Vec # and a direction, which has to be unitized, norm(dir) = 1
end
#=Ray(p::Vector{Float64},d::Vector{Float64}) = Ray(Vec(p),Vec(d))=#
#=Ray(p::Vec,d::Vector{Float64}) = Ray(p,Vec(d))=#
#=Ray(p::Vector{Float64},d::Vec) = Ray(Vec(p),d)=#

type Bull
    ray::Ray
    N::Vec
    dis::Float64
    n1::Float64
end

abstract Sphere

type Lens <: Sphere# the lens type
    c::Vec
    r::Float64
    ri::Function # the Refractive Index Gradient (RIG) that describes the distribution of refractive indices inside the lens. We restrict ourselves to spherically symmetrical distributions: so at a set distance form the lens center, the refractive index is equal in all directions (kind of like an onion?).
    ε::Float64 # the step
end
Lens(a::Vector{Float64},b::Float64,c::Function,d::Float64) = Lens(Vec(a),b,c,d) # this is a constructor function, I learned about it from the Julia manual. Useful for constructing instances of Lens

type Retina <: Sphere# the retina type, simple for now
    c::Vec
    r::Float64
end
Retina(a::Vector{Float64},b::Float64) = Retina(Vec(a),b)
#=We should probably define a few more helpful types in the near future...=#

#=Now come the functions=#

#=this function initiates the array of Rays. It starts them off from directly above the lens at distance L. Their directions are uniformly and randomly distributed across the lens.=#

function Bull(ray::Ray,lens::Lens)

    ray.pos -= lens.c
    return Bull(ray,Vec(0.,0.,0.),0.,0.)
#=    n1 = lens.ri(norm(ray.pos))
    ray.pos += lens.ε*ray.dir
    N = unitize(ray.pos)
    dis = norm(ray.pos)
    return Bull(ray,N,dis,n1)=#
end


function letBlight(lens::Lens,L::Float64)
    if isinf(L) # if the light source is infinitely far away
        #=just ditribute the rays uniformly on the disk that the sphere makes=#
        sqrtr = sqrt(rand()) 
        θ = rand()*2π
        a = Vec([sqrtr*cos(θ),sqrtr*sin(θ),2lens.r])
        pos = lens.c + a
        dir = Vec([0.,0.,-1])
            #=initiate the Ray with those starting points at distance 2*rad (could be any distance larger than rad really, and point them directly down at the lens=#
    else
        #=theta is the angle between the optical axis and the ray that touches the surface of the sphere=#
        cosα = cos(asin(lens.r/L))

        #=z and α are randomly sampled from those ranges to generate the uniformly spaced points on the sphere=#
        z = rand()*(cosα - 1.) - cosα
        θ = rand()*2π

        #=now come the x, y, and z that lie on the sphere's surface=#
        x = sqrt(1. - z*z)*cos(θ)
        y = sqrt(1. - z*z)*sin(θ)

        #=this is the starting position of all the rays -- point source=#
        a = [0.,0.,L] 
        pos = lens.c + a
        dir = unitize(Vec([x,y,z]))
    #=initiate the Ray with its initial position, s, and it's unitized direction that will end up on the sphere=#
    end
    return Ray(pos,dir)
end

#=This function finds the intersection point of the light as it starts from its initial point. It's a basic vector-sphere intersection thing. It works. I copied it from the ray_sphere function in https://github.com/JuliaLang/julia/blob/master/test/perf/kernel/raytracer.jl...=#
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

#=I found it useful to separate that advance function. The main reason is that refraction occurs at the surface between the medium and the lens, not only in the lens. So this function can be used in both scenarios. ray is the ray, N is the normal to the refractive surface, and n is the ratio between the refractive index before the surface and the refractive index after the surface. =#
function refract!(ray::Ray,N::Vec,n::Float64)
    # this is the 'r' in page 74 in that pdf I just mentioned
    a = -dot(ray.dir,N)
    
    # this is the 'c' in page 74
    c = 1. - n*n*(1. - a*a)
    
    # and finally the 't' in page 74
    # added reflection as well... See page 54 (v is -ray.dir)
    t = c>0 ? n*ray.dir + (n*a - sqrt(c))*N : 2* dot(N,-ray.dir) *N + ray.dir

    # I unitize the direction vector for the ray, updating the ray's direction
    ray.dir = unitize(t)
end

function poke!(bull::Bull,lens::Lens)
 # first I find the refractive index at the current position
    bull.n1 = lens.ri(norm(bull.ray.pos))

    # then I advance the ray one ε forward with its original direction
    bull.ray.pos += lens.ε*bull.ray.dir

    # useful
    bull.dis = norm(bull.ray.pos)
    bull.N = bull.ray.pos/bull.dis
end

# This is the main function that advances the ray along its path inside the lens. 
# Everything happens here, and in fact too much does. It would be better if we'd brake it 
# apart into logical segments. The algorithm is taken directly from page 74 in the 
# S1-rt.pdf file I got from http://fileadmin.cs.lth.se/cs/Education/EDAN30/lectures/S1-rt.pdf 
# in that website your graphics guy gave you 
# (i.e. http://cs.lth.se/english/course/edan30-photorealistic-computer-graphics/assignments/). 
function advance!(bull::Bull,lens::Lens) 
    # now I'm in the new position where the ray will refract. So there's a new refractive index here:
    n2 = lens.ri(bull.dis)

    # here I find out if the ray is on it's way in or way out. 
    # If its distance to the lens center is diminishing then the normal to the refractive 
    # surface is going outwards. But if the ray is leaving the lens the normal to 
    # the surface is pointing inwards.
    costheta = dot(bull.ray.dir,bull.ray.pos)
    
    # is it going out? If so make out equal to +1, otherwise make it -1
    out = -sign(costheta)
    
    # this is the normal to the refractive surface 
    bull.N = out*bull.ray.pos/bull.dis
    
    # refraction, i.e. changing the direction of the ray, occurs with this new function
    refract!(bull.ray,bull.N,bull.n1/n2)

end


const ε = 1e-4 # the ε size with which the ray advances every iteration
const lens_r = 1. # the lens radius, this is only useful cause there are a couple of calculations that depend on this variable
const c = [0.,0.,0.] # this is the location of the center of the lens
const retina_r = 1.2*lens_r # this is the retina's radius, good practice to make it a function of the lens radius 
const L = Inf#1e3*lens_r # this is the distance between the center of the lens and the source light
const nrays = 100 # the number of discrete rays I'll be tracing
const n_medium = 1. # this is the refractive index of the medium surrounding the lens
#=this is the RIG function. Notice that its argument is now Real -- this is in the interest of generality... I set this RIG to be equal to the Luneburg lens, so that we can check the results real quick. see: http://en.wikipedia.org/wiki/Luneburg_lens=#
lens_r2 = lens_r*lens_r
ri(r::Float64) = sqrt(2. - r*r/lens_r2)

#=OK, now comes the actual calculations=#
lens = Lens(c,lens_r,ri,ε) # initiate the lens, this is how we initiate an instance of a type, it's like a function!
retina = Retina(c,retina_r) # initiate the retina

p = [Array(Vec,1)::Array{Vec,1} for i in 1:nrays] # a cell array of arrays of positions. Defined as an array who's elements are Vec
ray = Ray(Vec(zeros(3)),Vec(zeros(3)))
nsteps = 0
#=Now we iterate through all the rays, tracing them=#
for i = 1:nrays
    
    ray = letBlight(lens,L) # initiate all the rays
    p[i][1] = ray.pos # populate with the start position

    #=find the intersection point of the ray with the lens and update ray's position=#
    intersect!(ray,lens)
    push!(p[i],ray.pos) # add said position to the position arrays
    refract!(ray,unitize(ray.pos - lens.c),n_medium/lens.ri(lens.r)) # refract at the surface of the lens due to the medium's refractive index (that might be different from that of the lens' periphery)
    bull = Bull(ray,lens)
    nsteps = 0
    #=advance!(ray,lens) # advance once, just to penetrate the lens and be able to start that while loop below=#
    #=push!(p[i],ray.pos) # add said position to the position arrays=#
    while nsteps < 1e6 #length(p[i]) < 1e5 # until you looped more than # times (totally arbitrary, and useful only in cases where the rays are trapped in the lens)


        poke!(bull,lens) # advance one ε
        if bull.dis > lens.r # is the ray GOING TO exit the lens in the next iteration? If so, intersect it with the lens (which is therefore closer)
            ray = bull.ray
            ray.pos -= lens.ε*ray.dir-lens.c
            intersect!(ray,lens) # find intersection point with the lens, update Ray ray
            push!(p[i],ray.pos) # add said position to the position arrays
            break # and break the loop
        end
        advance!(bull,lens)
        push!(p[i],bull.ray.pos + lens.c) # save the new position in the positions array
        nsteps += 1
    end
    refract!(ray,-unitize(ray.pos - lens.c),lens.ri(lens.r)/n_medium) # refract when exiting the lens and entering the surrounding medium
    intersect!(ray,retina) # find intersection point with the retina, update Ray ray
    push!(p[i],ray.pos) # add said position to the position arrays
end

#=plot all the ray's trajectories<]=#
for i = 1:nrays
  x = map(t -> t.x,p[i])
  y = map(t -> t.y,p[i])
  z = map(t -> t.z,p[i])

  shift!(x) # Now I'm removing the point source from the plotting so that it won't look silly when we have L > say 3 (like in the case of L = 1e30 or something)
  shift!(y)
  shift!(z)

  plot3D(x,y,z,color="red") # plot the ray trajectory in red
end

#>now comes the shitty sphere plot. I think I improved it a bit :)<] =#
n = 100 # number of points to plot=#
u = linspace(0., π, n) # the two angles that define the location of the points on the sphere=#
v = linspace(0., 2π, n)' # notice that I transpose this one with '=#
u .+= 0*v # the transposition is "carried" and expanded through the .+ operator to make u a n by n matrix=#
v .+= 0*u # same here=#
x = lens.c.x + lens.r*cos(u).*sin(v) # then standard vector math=#
y = lens.c.y + lens.r*sin(u).*sin(v)
z = lens.c.z + lens.r*cos(v)
plot_surface(x, y, z, color="blue",alpha=.5,linewidth=0,rstride=4,cstride=4) # and finally plotting the sphere with some transparency so that we'll be able to see the ray=#

#>This is for plotting a part of the retina, like a hemisphere, but not hemi... This code is exactly the same as for the initiation of the rays, when they are gonna hit the lens<]=#
θ = .5 # this is the angle of the sphere-retina that will be plotted=#

#>z and h are uniformly sampled from those ranges to generate the uniformly spaced points on the sphere<]=#
z = linspace(cos(θ),1.,n)
h = linspace(0.,2π,n)'

#>now come the x, y, and z that lie on the retina's surface<]=#
x = retina.c.x + retina.r*sqrt(1. - z.^2).*cos(h)
y = retina.c.y + retina.r*sqrt(1. - z.^2).*sin(h)
z = retina.c.z - retina.r*z .+ 0.*h # I flip the z so it'll point downwards=#

plot_surface(x, y, z, color="green",alpha=.25,linewidth=0,rstride=4,cstride=4) # and finally plotting the sphere with some transparency so that we'll be able to see the ray=#

#>Thought I'd assess the accuracy of the focal point. The points in p one before last are all at the lens surface, and because the rig is a Luneburg lens, then they should all be lens radius below the lens origin. Feel free to improve on the following:<]=#
map(x -> [x[end-1]],p)
a = [x[end-1] for x in p]
a = [[x.x,x.y,x.z] for x in a]
a = cat(2,a...)
goal = [lens.c.x,lens.c.y,lens.c.z] .- [0.,0.,lens.r] # this is where they should all be at=#
σ = sum((a .- repmat(goal,1,nrays)).^2)/nrays

#=end=#
