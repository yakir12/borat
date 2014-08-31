using ImmutableArrays, PyPlot # useful library for fast matrix manipulations, plotting in 3D

#=first we define a few useful types. Since Julia is a highly (but not strictly) typed language, this is an awesome feature that speeds up the code and makes it clean and easy to follow=#
type Ray # the ray type
    pos::Vector3 # has a position, a location
    dir::Vector3 # and a direction, which has to be unitized, norm(dir) = 1

    #=see 1.13.2 Inner Constructor Methods. I'm tempted to just do this:=#
    #=Ray(pos,dir) = new(pos,unit(dir))=#
    #=but they say it's bad form???=#
    Ray(pos,dir) = norm(dir) != 1. ? error("directional vector not unitized") : new(pos,dir) 
end

type Sphere # this is going to be used by both the lens and the retina, since they are both sphere shaped
    org::Vector3
    rad::Number
    #=Easy to avoid, but why not?=#
    Sphere(org,rad) = rad < 0. ? error("negative sphere radius") : new(org,rad)
end

type Lens # the lens type
    s::Sphere
    rig::Function # the Refractive Index Gradient (RIG) that describes the distribution of refractive indices inside the lens. We restrict ourselves to spherically symmetrical distributions: so at a set distance form the lens center, the refractive index is equal in all directions (kind of like an onion?).
end
Lens(o,r,f) = Lens(Sphere(Vector3(o),r),f) # this is a constructor function, I learned about it from the Julia manual. Useful for constructing instances of Lens

type Retina # the retina type, simple for now
    s::Sphere
end
Retina(o,r) = Retina(Sphere(Vector3(o),r))
#=We should probably define a few more helpful types in the near future...=#

#=Now come the constants, there is no real reason to define the following variables as constants, other than that it makes sense that these don't change in value throughout the program. =#
const step = 1e-4 # the step size with which the ray advances every iteration
const lens_r = 1. # the lens radius, this is only useful cause there are a couple of calculations that depend on this variable
const L = Inf#1e3*lens_r # this is the distance between the center of the lens and the source light
const nrays = 100 # the number of discrete rays I'll be tracing
const c = Vector3([0.,0.,0.]) # this is the location of the center of the lens
const n_medium = 1. # this is the refractive index of the medium surrounding the lens
const retina_r = 1.2*lens_r # this is the retina's radius, good practice to make it a function of the lens radius 

#=Now come the functions=#
#=this is the RIG function. Notice that its argument is now Number -- this is in the interest of generality... I set this RIG to be equal to the Luneburg lens, so that we can check the results real quick. see: http://en.wikipedia.org/wiki/Luneburg_lens=#
rig(r::Number) = sqrt(2. - (r/lens_r)^2)

#=this function initiates the array of Rays. It starts them off from directly above the lens at distance L. Their directions are uniformly and randomly distributed across the lens.=#
function initiate(l::Lens,L::Number,n::Int)
    if isinf(L) # if the light source is infinitely far away
        #=just ditribute the rays uniformly on the disk that the sphere makes=#
        sqrtr = sqrt(rand(n)) 
        θ = rand(n)*2π
        x = l.s.org[1] + sqrtr .* cos(θ)
        y = l.s.org[2] + sqrtr .* sin(θ)
        #=an empty container for all the Rays=#
        r = Array(Ray,n)
        for i = 1:n
            #=initiate the Ray with those starting points at distance 2*rad (could be any distance larger than rad really, and point them directly down at the lens=#
            r[i] = Ray(Vector3([x[i],y[i],l.s.org[3] + 2l.s.rad]),Vector3([0.,0.,-1.]))
        end
    else
        #=theta is the angle between the optical axis and the ray that touches the surface of the sphere=#
        θ = asin(l.s.rad/L)

        #=epsilon is used only for the range in the rand function below, you'll see=#
        ε = eps(typeof(θ))

        #=z and h are randomly sampled from those ranges to generate the uniformly spaced points on the sphere=#
        z = rand(cos(θ):ε:1.,n)
        h = rand(0.:ε:2π,n)

        #=z = linspace(cos(θ),1.,n)=#
        #=h = linspace(0.,2π,n)=#

        #=now come the x, y, and z that lie on the sphere's surface=#
        x = sqrt(1. - z.^2).*cos(h)
        y = sqrt(1. - z.^2).*sin(h)
        z *= -1. # I flip the z so it'll point downwards

        #=this is the starting position of all the rays -- point source=#
        s = Vector3([0.,0.,L] .+ l.s.org)

        #=an empty container for all the Rays=#
        r = Array(Ray,n)
        for i = 1:n
            #=initiate the Ray with its initial position, s, and it's unitized direction that will end up on the sphere=#
            r[i] = Ray(s,unit(Vector3([x[i],y[i],z[i]])))
        end
    end
    return r
end

#=This function finds the intersection point of the light as it starts from its initial point. It's a basic vector-sphere intersection thing. It works. I copied it from the ray_sphere function in https://github.com/JuliaLang/julia/blob/master/test/perf/kernel/raytracer.jl...=#
function intersect!(r::Ray,s::Sphere)
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

#=I found it useful to separate that advance function. The main reason is that refraction occurs at the surface between the medium and the lens, not only in the lens. So this function can be used in both scenarios. r is the ray, N is the normal to the refractive surface, and n is the ratio between the refractive index before the surface and the refractive index after the surface. =#
function refract!(r::Ray,N::Vector3,n::Number)
    # this is the 'r' in page 74 in that pdf I just mentioned
    a = -dot(r.dir,N)
    
    # this is the 'c' in page 74
    c = 1. - n^2*(1. - a^2)
    
    # and finally the 't' in page 74
    # added reflection as well... See page 54 (v is -r.dir)
    t = c>0 ? n*r.dir + (n*a - sqrt(c))*N : 2* dot(N,-r.dir) *N + r.dir

    # I unitize the direction vector for the ray, updating the ray's direction
    r.dir = unit(t)
end


# This is the main function that advances the ray along its path inside the lens. 
# Everything happens here, and in fact too much does. It would be better if we'd brake it 
# apart into logical segments. The algorithm is taken directly from page 74 in the 
# S1-rt.pdf file I got from http://fileadmin.cs.lth.se/cs/Education/EDAN30/lectures/S1-rt.pdf 
# in that website your graphics guy gave you 
# (i.e. http://cs.lth.se/english/course/edan30-photorealistic-computer-graphics/assignments/). 
function advance!(r::Ray,l::Lens) 
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

#=OK, now comes the actual calculations=#
l = Lens(c,lens_r,rig) # initiate the lens, this is how we initiate an instance of a type, it's like a function!
retina = Retina(c,retina_r) # initiate the retina
r0 = initiate(l,L,nrays) # initiate all the rays

p = [Array(Vector3,1)::Array{Vector3,1} for i in 1:nrays] # a cell array of arrays of positions. Defined as an array who's elements are Vector3

#=Now we iterate through all the rays, tracing them=#
for i = 1:nrays
    
    #=I first copy the i^th ray's initial position to r, for convenience=#
    r = deepcopy(r0[i])
    p[i][1] = r.pos # populate with the start position

    #=find the intersection point of the ray with the lens and update r's position=#
    intersect!(r,l.s)
    push!(p[i],r.pos) # add said position to the position arrays
    refract!(r,unit(r.pos - l.s.org),n_medium/l.rig(l.s.rad)) # refract at the surface of the lens due to the medium's refractive index (that might be different from that of the lens' periphery)
    advance!(r,l) # advance once, just to penetrate the lens and be able to start that while loop below
    push!(p[i],r.pos) # add said position to the position arrays
    while length(p[i]) < 1e5 # until you looped more than # times (totally arbitrary, and useful only in cases where the rays are trapped in the lens)
        advance!(r,l) # advance one step
        push!(p[i],r.pos) # save the new position in the positions array
        if norm(r.pos + step*r.dir-l.s.org) > l.s.rad # is the ray GOING TO exit the lens in the next iteration? If so, intersect it with the lens (which is therefore closer)
            intersect!(r,l.s) # find intersection point with the lens, update Ray r
            push!(p[i],r.pos) # add said position to the position arrays
            break # and break the loop
        end
    end
    refract!(r,-unit(r.pos - l.s.org),l.rig(l.s.rad)/n_medium) # refract when exiting the lens and entering the surrounding medium
    intersect!(r,retina.s) # find intersection point with the retina, update Ray r
    push!(p[i],r.pos) # add said position to the position arrays
end

#=plot all the ray's trajectories=#
for i = 1:nrays
    x = map(t -> t[1],p[i])
    y = map(t -> t[2],p[i])
    z = map(t -> t[3],p[i])
  
    shift!(x) # Now I'm removing the point source from the plotting so that it won't look silly when we have L > say 3 (like in the case of L = 1e30 or something)
    shift!(y)
    shift!(z)

    plot3D(x,y,z,color="red") # plot the ray trajectory in red
end

#=now comes the shitty sphere plot. I think I improved it a bit :)=# 
n = 100 # number of points to plot
u = linspace(0., π, n) # the two angles that define the location of the points on the sphere
v = linspace(0., 2π, n)' # notice that I transpose this one with '
u .+= 0*v # the transposition is "carried" and expanded through the .+ operator to make u a n by n matrix
v .+= 0*u # same here
x = l.s.org.e1 + l.s.rad*cos(u).*sin(v) # then standard vector math
y = l.s.org.e2 + l.s.rad*sin(u).*sin(v)
z = l.s.org.e3 + l.s.rad*cos(v)
plot_surface(x, y, z, color="blue",alpha=.5,linewidth=0,rstride=4,cstride=4) # and finally plotting the sphere with some transparency so that we'll be able to see the ray

#=This is for plotting a part of the retina, like a hemisphere, but not hemi... This code is exactly the same as for the initiation of the rays, when they are gonna hit the lens=#
θ = .5 # this is the angle of the sphere-retina that will be plotted

#=epsilon is used only for the range in the rand function below, you'll see=#
ε = eps(typeof(θ))

#=z and h are uniformly sampled from those ranges to generate the uniformly spaced points on the sphere=#
z = linspace(cos(θ),1.,n)
h = linspace(0.,2π,n)'

#=now come the x, y, and z that lie on the retina's surface=#
x = retina.s.org.e1 + retina.s.rad*sqrt(1. - z.^2).*cos(h)
y = retina.s.org.e1 + retina.s.rad*sqrt(1. - z.^2).*sin(h)
z = retina.s.org.e1 - retina.s.rad*z .+ 0.*h # I flip the z so it'll point downwards

plot_surface(x, y, z, color="green",alpha=.25,linewidth=0,rstride=4,cstride=4) # and finally plotting the sphere with some transparency so that we'll be able to see the ray

#=Thought I'd assess the accuracy of the focal point. The points in p one before last are all at the lens surface, and because the rig is a Luneburg lens, then they should all be lens radius below the lens origin. Feel free to improve on the following:=#
map(x -> [x[end-1]],p)
a = [[x[end-1]] for x in p]
a = cat(2,a...)
goal = [l.s.org] .- [0.,0.,l.s.rad] # this is where they should all be at
σ = sum((a .- repmat(goal,1,nrays)).^2)/nrays
