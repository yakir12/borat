#= This is what I need to do in the Julia runtime environment to get things running:
cd("Documents/Work/RayTracing/borat")
include("main.jl")
=#

require("lightray.jl")
require("geometry.jl")

# useful library for fast matrix manipulations, plotting in 3D
@everywhere using ImmutableArrays, PyPlot 

# our own stuff
@everywhere using Geometry, LightRay

@everywhere begin
#=Now come the constants, there is no real reason to define the following variables as constants, other than that it makes sense that these don't change in value throughout the program. =#
const step = 1e-4 # the step size with which the ray advances every iteration
const lens_r = 1. # the lens radius, this is only useful cause there are a couple of calculations that depend on this variable
const L = Inf#1e3*lens_r # this is the distance between the center of the lens and the source light
const nrays = 100 # the number of discrete rays I'll be tracing
const c = Vector3([0.,0.,0.]) # this is the location of the center of the lens
const c_light = c + Vector3([0.,0.,L]) # light source location
const n_medium = 1. # this is the refractive index of the medium surrounding the lens
const retina_r = 1.2*lens_r # this is the retina's radius, good practice to make it a function of the lens radius 

const procsToUse = 8
# Now come the functions

# This is the RIG (refractive index gradient) function. 
# Notice that its argument is now Real -- this is in the interest of generality... 
# I set this RIG to be equal to the Luneburg lens, so that we can check the results 
# real quick. see: http://en.wikipedia.org/wiki/Luneburg_lens
rig(r::Real) = sqrt(2. - (r/lens_r)^2)

# The actual calculations 
# initiate the lens, this is how we initiate an instance of a type, it's like a function!
const l = Lens(c,lens_r,rig)

# initiate the retina
const retina = Retina(c,retina_r) 

end

# helper, plots one ray
function plotOneRay(oneray::Array{Vector3{Float64},1})
    x = map(t -> t[1],oneray)
    y = map(t -> t[2],oneray)
    z = map(t -> t[3],oneray)
  
    shift!(x) # Now I'm removing the point source from the plotting so that it won't look silly when we have L > say 3 (like in the case of L = 1e30 or something)
    shift!(y)
    shift!(z)

    plot3D(x,y,z,color="red") # plot the ray trajectory in red
end

# define the function for tracing one ray
@everywhere function traceOneRay(l::Lens,L::Real) #:Array{Vector3, 1}
    # initiate this ray
    r = initiate(c_light, L, lens_r)
    
    # initialise position array
    oldpi = Array(Vector3{Float64},1)
	oldpi[1] = r.pos
	
    # intersect with lens
    intersect!(r,l)
    push!(oldpi, r.pos)
    # refract at the surface of the lens due to the medium's refractive index 
    # (that might be different from that of the lens' periphery)
    refract!(r,unit(r.pos - l.s.org),n_medium/l.rig(l.s.rad))
    
    # advance once, just to penetrate the lens and be able to start that while loop below     
    advance!(r,l,step)
	push!(oldpi, r.pos)

	# until you looped more than # times 
	# (totally arbitrary, and useful only in cases where the rays are trapped in the lens)
    while length(oldpi) < 1e5 
    	# advance one step and save
        advance!(r,l,step) 
        push!(oldpi, r.pos)
        
        # is the ray GOING TO exit the lens in the next iteration? 
        # If so, intersect it with the lens (which is therefore closer)
        if norm(r.pos + step*r.dir-l.s.org) > l.s.rad
            # find intersection point with the lens and save position
            intersect!(r,l)
            push!(oldpi, r.pos)
            
            # the ray got out, stop tracing it
            break 
        end
    end
    
    # refract when exiting the lens and entering the surrounding medium
    refract!(r,-unit(r.pos - l.s.org),l.rig(l.s.rad)/n_medium)
    
    # find intersection point with the retina, and save pos    
    intersect!(r,retina) 
    push!(oldpi, r.pos)

	oldpi
end


# set this up ONCE for a parallel environment
np = nprocs()
if np >= procsToUse
	println("Already running on ",np," processors.")
else
	deltaproc = procsToUse-np
	println("Starting ",deltaproc," more processors.")
	addprocs(deltaproc)
	println("Now running on ",nprocs()," processors.")
end

# testing pmap
rayindices = [i::Int for i in 1:nrays]

# helper fun
@everywhere function tOR(j::Int)
#	id = myid()
#	println(id)
	traceOneRay(l,L)
end

# start timer
tic()

# run ray tracing in parallel
pm = pmap(tOR, rayindices)

# print time from tic
toc()

# do the plot, pmap?
map(plotOneRay,pm)

#=
#@sync @parallel for i = 1:nrays
for i = 1:nrays
	ar = tOR(i)
	plotOneRay(ar)
end
=#

# now comes the shitty sphere plot. I think I improved it a bit :)=# 
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

p = pm
#=Thought I'd assess the accuracy of the focal point. The points in p one before last are all at the lens surface, and because the rig is a Luneburg lens, then they should all be lens radius below the lens origin. Feel free to improve on the following:=#
map(x -> [x[end-1]],p)
a = [[x[end-1]] for x in p]
a = cat(2,a...)
goal = [l.s.org] .- [0.,0.,l.s.rad] # this is where they should all be at
σ = sum((a .- repmat(goal,1,nrays)).^2)/nrays
