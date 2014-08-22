#=useful source: https://github.com/JuliaLang/julia/blob/master/test/perf/kernel/raytracer.jl=#
using ImmutableArrays, PyPlot # useful library for fast matrix manipulations, ploting in 3D

type Ray # the ray type
    pos::Vector3 # has an pos, or location
    dir::Vector3 # and a direction, which has to be unitized
end

type Lens
    orig::Vector3
    rad::Float64
end

const step = 1e-3 # step size

function initiate(o::Array{Float64},d::Array{Float64})
    # just a utility function to declare a ray
    r = Ray(Vector3(o),unit(Vector3(d)))
end

function advance!(r::Ray)
    # simple but useful
    r.pos += step*r.dir
end

r1 = initiate([0.,0.,1.],[0.,0.,-1.]) # start a ray at 0,0,0 and give it some direction
p1 = {r1.pos}
l = Lens(Vector3([0.,0.,0.]),1.)
advance!(r1)
push!(p1,r1.pos)
while norm(r1.pos) < l.rad
    advance!(r1)
    push!(p1,r1.pos)
end

plot3D(map(x->x.e1,p1),map(x->x.e2,p1),map(x->x.e3,p1))

n = 100
u = linspace(0., pi, n)
v = linspace(0., 2pi, n)'
u .+= 0*v
v .+= 0*u

x = l.orig.e1 + l.rad*cos(u).*sin(v)
y = l.orig.e2 + l.rad*sin(u).*sin(v)
z = l.orig.e3 + l.rad*cos(v)
plot_surface(x, y, z, color="green",alpha=.5,linewidth=0)
