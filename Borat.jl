using PyPlot, Shapes, RayTrace

const maxiter = int(1e5)
const ε = 1e-4 # the ε size with which the ray bend every iteration
const lens_r = 1. # the lens radius, this is only useful cause there are a couple of calculations that depend on this variable
const c = [0.,0.,0.] # this is the location of the center of the lens
const retina_r = 1.2*lens_r # this is the retina's radius, good practice to make it a function of the lens radius 
const L = Inf#1e3*lens_r # this is the distance between the center of the lens and the source light
const nrays = 100 # the number of discrete rays I'll be tracing
const n_medium = 1. # this is the refractive index of the medium surrounding the lens
const lens_r2 = lens_r*lens_r
const goal = [0.,0.,-lens_r]
ri(r::Float64) = sqrt(2. - r*r/lens_r2) # this is the RIG function. Notice that its argument is now Real -- this is in the interest of generality... I set this RIG to be equal to the Luneburg lens, so that we can check the results real quick. see: http://en.wikipedia.org/wiki/Luneburg_lens

#=OK, now comes the actual calculations=#
const lens = Lens(c,lens_r,ri,ε) # initiate the lens, this is how we initiate an instance of a type, it's like a function!
const retina = Retina(c,retina_r) # initiate the retina

P = zeros(nrays)
nsteps = 0
#=Now we iterate through all the rays, tracing them=#
for i = 1:nrays
    p = zeros(3,maxiter)
    Δ = [0.]
    iter = [0]
    trace!(p,lens,retina,L,n_medium,maxiter,goal,Δ,iter)
    P[i] = Δ[1]
    plot3D(vec(p[1,1:iter[1]]),vec(p[2,1:iter[1]]),vec(p[3,1:iter[1]]),color="red")
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
σ = mean(P)
