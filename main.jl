using PyPlot # would be nice to not need to load this on every proccessor
@everywhere using Shapes, RayTrace # PyPlot is slow and in the ultimate version will be removed
@everywhere begin
# constants:
const ε = 1e-4 # the ε size with which the ray advances every iteration
const lens_r = 1. # the lens radius
const maxiter = int(2π*lens_r/ε) # maximum iteration allowed. This value is based on the arc-length of an imaginary ray that travels along the surface of the lens, times two! So should be safe.
const c = [100.,-20.,330.] # location of the center of the lens. It has no real use here, but good to have it.
const retina_r = 1.2*lens_r # the retina's radius, good practice to make it a function of the lens radius 
const L = Inf # the distance between the center of the lens and the source light
const nrays = 1000 # the number of discrete rays
const n_medium = 1. # the refractive index of the medium surrounding the lens
const lens_r2 = lens_r*lens_r # just to make the refractive index gradient function faster
const goal = c .- [0.,0.,lens_r] # this is where all rays that start at L=Inf should end up at.
ri(r::Float64) = sqrt(2. - r*r/lens_r2) # this is the refractive index gradient function. I set this RIG to be equal to the Luneburg lens, so that we can check the results real quick.
const lens = Lens(c,lens_r,ri,ε) # initiate the lens
const retina = Retina(c,retina_r) # initiate the retina
#const procsToUse = 2

# actual calculations

end

@everywhere function oneIter(i)
# iterate through all the rays, ploting them
    p = zeros(3,maxiter) # ray locations
    diagnose = [:Δ => 0., :len => 0, :goal => goal] # diagnosys, including the distance to the goal, the number of refraction points, the goal
    trace!(p,lens,retina,L,n_medium,maxiter,diagnose) # main function that ray-traces the system
    return (p,diagnose)
    # plotting:
    #ind = round(linspace(2,diagnose[:len],10)) # let's plot just 10 points inside the lens, no need to go crazy
    #plot3D(vec(p[1,ind]),vec(p[2,ind]),vec(p[3,ind]),color="red") # plotting it in red...
end

x = pmap(oneIter,1:nrays)


Δ = zeros(nrays) # the matrix of euclidean distance between where the ray exsits the lens and the goal
for i = 1:nrays
    p = x[i][1]
    diagnose = x[i][2]
    Δ[i] = diagnose[:Δ] # put that distance in the matrix
    ind = round(linspace(2,diagnose[:len],10)) # let's plot just 10 points inside the lens, no need to go crazy
    plot3D(vec(p[1,ind]),vec(p[2,ind]),vec(p[3,ind]),color="red") # plotting it in red...
end

# set this up ONCE for a parallel environment
#=np = nprocs()
if np >= procsToUse
	println("Already running on ",np," processors.")
else
	deltaproc = procsToUse-np
	println("Starting ",deltaproc," more processors.")
	addprocs(deltaproc)
	println("Now running on ",nprocs()," processors.")
end=#

# sphere plot
n = 100 # number of points to plot=#
u = linspace(0., π, n) # the two angles that define the location of the points on the sphere
v = linspace(0., 2π, n)' # notice that I transpose this one with '
u .+= 0*v # the transposition is "carried" and expanded through the .+ operator to make u a n by n matrix
v .+= 0*u # same here
x = lens.c.x + lens.r*cos(u).*sin(v) # then standard vector math
y = lens.c.y + lens.r*sin(u).*sin(v)
z = lens.c.z + lens.r*cos(v)
plot_surface(x, y, z, color="blue",alpha=.5,linewidth=0,rstride=4,cstride=4) # and finally plotting the sphere with some transparency so that we'll be able to see the ray

# plotting a part of the retina, like a hemisphere, but not hemi... This code is exactly the same as for the initiation of the rays, when they are gonna hit the lens
θ = .5 # this is the (arbitrary) angle of the sphere-retina that will be plotted

# z and h are uniformly sampled from these ranges to generate the uniformly spaced points on the sphere
z = linspace(cos(θ),1.,n)
h = linspace(0.,2π,n)'

# now come the x, y, and z that lie on the retina's surface
x = retina.c.x + retina.r*sqrt(1. - z.^2).*cos(h)
y = retina.c.y + retina.r*sqrt(1. - z.^2).*sin(h)
z = retina.c.z - retina.r*z .+ 0.*h # I flip the z so it'll point downwards

plot_surface(x, y, z, color="green",alpha=.25,linewidth=0,rstride=4,cstride=4) # and finally plotting the sphere with some transparency

σ = mean(Δ) # mean error
