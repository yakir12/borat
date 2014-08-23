#=useful source: https://github.com/JuliaLang/julia/blob/master/test/perf/kernel/raytracer.jl=#
using ImmutableArrays, PyPlot # useful library for fast matrix manipulations, ploting in 3D

type Ray # the ray type
    pos::Vector3 # has an pos, or location
    dir::Vector3 # and a direction, which has to be unitized
end

type Lens
    orig::Vector3
    rad::Float64
    rig::Function
end

const step = 1e-5 # step size
const rad = 1.
const n_periphery = 1.34
const n_center = 1.55

rig(r::Float64) = (n_center-n_periphery)*(rad-r).*(rad+r)/rad^2+n_periphery

function initiate(o::Array{Float64},d::Array{Float64})
    # just a utility function to declare a ray
    r = Ray(Vector3(o),unit(Vector3(d)))
end

#=function refract!(r::Ray,l::Lens)=#
    
 
function advance!(r::Ray,l::Lens)
    N = unit(r.pos)
    a = dot(-r.dir,N)
    n1 = l.rig(norm(r.pos))
    n2 = l.rig(norm(r.pos + step*r.dir))
    n = n1/n2
    c = 1. - n^2(1.-a^2)
    if c < 0
        c = abs(c)
    end

    t = n*r.dir + (n*a - sqrt(c))*N
    r.dir = unit(t)
    r.pos += step*r.dir
end

r1 = initiate([.5,.5,.5],[0.,0.,-1.]) # start a ray at 0,0,0 and give it some direction
p1 = {r1.pos}
l = Lens(Vector3([0.,0.,0.]),rad,rig)
advance!(r1,l)
push!(p1,r1.pos)
while (norm(r1.pos) < l.rad) & (length(p1) < 1e6)
    advance!(r1,l)
    push!(p1,r1.pos)
end

plot3D(map(x->x.e1,p1),map(x->x.e2,p1),map(x->x.e3,p1),color="red")

n = 100
u = linspace(0., pi, n)
v = linspace(0., 2pi, n)'
u .+= 0*v
v .+= 0*u

x = l.orig.e1 + l.rad*cos(u).*sin(v)
y = l.orig.e2 + l.rad*sin(u).*sin(v)
z = l.orig.e3 + l.rad*cos(v)
plot_surface(x, y, z, color="blue",alpha=.5,linewidth=0)
