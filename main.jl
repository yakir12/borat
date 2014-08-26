using ImmutableArrays, PyPlot # useful library for fast matrix manipulations, ploting in 3D

#=first we define a few useful types. Since Julia is a highly (but not stricktly) typed language, this is an awesome feature that speeds up the code and makes it clean and easy to follow=#
type Ray # the ray type
    pos::Vector3 # has a position, a location
    dir::Vector3 # and a direction, which has to be unitized, norm(dir) = 1
end

type Lens # the lens type, we restrict ourselves to sphere shaped lenses
    orig::Vector3 # the location of the sphere center
    rad::Float64 # the radius of the sphere
    rig::Function # the Refractive Index Gradient (RIG) that describes the distribution of refractive indices inside the lens. We restrict ourdelves to spherically symmetrical distributions: so at a set distance form the lens center, the refractive index is equal in all directions (kind of like an onion?).
end
#=We should probably define a few more helpful types in the near future...=#

#=Now come the constants, there is no real reason to define the following variables as constants, other than that it makes sense that these don't change in value throughout the program. =#
const step = 1e-2 # the step size withwhich the ray advances every iterration
const rad = 1. # the lens radius, this is only usfule cause there are a couple of calculations that depend on htis variable
const n_periphery = 1.34 # the refractive index at the periphery of the lens, so at the very surface of the lens, when the distance from the lens center is equal to rad
const n_center = 1.35 # the refractive index at the center of the lens, when the distance from the lens center is equal to zero

#=Now come the functions=#
rig(r::Float64) = (n_center-n_periphery)*(rad-r)*(rad+r)/rad^2+n_periphery # this is the RIG function. Notice that its argument is of type Float64. It's just a parabula. Could be anything else really. It looks the way it does to make sure that it's a concave parabula with a n_center at r=0 and n_periphery at r=rad. this can be optimized a bit by pre calculating parts of that equation

function initiate(o::Array{Float64},d::Array{Float64}) # just a utility function to declare a ray. This is useful cause we can declare rays using a 3 element arrays (which isn't YET a Vector3D type), plus it takes care of the unitization of the direction vector
    r = Ray(Vector3(o),unit(Vector3(d)))
end

function advance!(r::Ray,l::Lens) # this is the main function that advances the ray along its path inside the lens. Everything happens here, and in fact too much does. It would be better if we'd brake it apart into logical segments. The algorithm is taken directly from page 74 in the S1-rt.pdf file I got from http://fileadmin.cs.lth.se/cs/Education/EDAN30/lectures/S1-rt.pdf in that website your graphics guy gave you (i.e. http://cs.lth.se/english/course/edan30-photorealistic-computer-graphics/assignments/). I'll refer to that page during the comments below. Notice that I didn't shifted the system to the lens center as I should. Since it's at 0,0 it doens't matter, but this should be included.
    costheta = dot(r.dir,r.pos) # here I find out if the ray is entering or exiting the lens (is it on it's way in or way out). If its distance to the lens center is diminishing then the normal to the refractive surface is going outwards. but if the ray is leaving the lens the normal to the surface is pointing inwards.
    out = -sign(costheta) # is it going out? if so make out equal to +1, otherwise make it -1
    N = out*unit(r.pos) # this is the normal to the refractive surface
    a = dot(-r.dir,N) # this is the 'r' in page 74 in that pdf I just mentioned
    n1 = l.rig(norm(r.pos + out*step/2*unit(r.pos))) # this is the refractive index before the current position of the ray
    n2 = l.rig(norm(r.pos - out*step/2*unit(r.pos))) # and this is the refractive index after. this is problematic...
    n = n1/n2 # this is the ration between the two and the 'η' in page 74
    c = 1. - n^2*(1. - a^2) # this is the 'c' in page 74
    t = n*r.dir + (n*a - sqrt(c))*N # and finally the 't' in page 74
    r.dir = unit(t) # I unitize the direction vector for the ray, updatiung the ray's direction
    r.pos += step*r.dir # and step the ray forward along said direction, updatiung the ray's position
end

#=OK, now comes the actual calculations=#
x = y = .9*rad/sqrt(2) # I define the x and y for the location of the ray, arbitrary choosing them to be equal to each other (in the near future we'll define a start point for the light to shoot at the lens, and find where in the lens surface the light intersect with the lens)
z = sqrt(rad^2-x^2-y^2) # z is calculated so that the ray will be on the lens surface
r = initiate([x,y,z],[0.,0.,-1.]) # start a ray at x,y,z and give it a downward direction
p1 = Array(Vector3,1) # an array of positions. Defined as an array who's elements are Vector3
p1[1] = r.pos # populate with the start position
l = Lens(Vector3([0.,0.,0.]),rad,rig) # initiate the lens, this is how we initiate an instance of a type, it's like a function!
advance!(r,l) # advance once, just to penetrate the lens and be able to start that while loop below
push!(p1,r.pos) # add said position to the position arrays
while (norm(r.pos) < l.rad) & (length(p1) < 2e2) # loop until the distance of the ray form the lens center is larger or equal to rad OR until you looped more than 200 times (totally arbitrary)
    advance!(r,l) # advance one step
    push!(p1,r.pos) # save the new position in the positions array
end

#=next we extract all the x,y,z values from the locations, maybe there's a more Julainian way to do this.=#
xyz = cat(2,p1...) # here I use the ... operator to extract all the vectors from the locations and concatinate them along the second dimension
x = vec(xyz[1,:]) # I convert each to a vector (in Julia, a vector is a stricktly one dimentional vector, which is different than a matirx with dimensions [1,2] or [2,1])
y = vec(xyz[2,:])
z = vec(xyz[3,:])
plot3D(x,y,z,color="red") # plot the ray trajectory in red

#=now comes the shitty sphere plot. surely there is a better way?=#
n = 100 # number of points to plot
u = linspace(0., π, n) # the two angles that define the location of the points on the sphere
v = linspace(0., 2π, n)' # notice that I transpose this one with '
u .+= 0*v # the transposition is "carried" and expanded through the .+ operator to make u a n by n matrix
v .+= 0*u # same here
x = l.orig.e1 + l.rad*cos(u).*sin(v) # then standard vector math
y = l.orig.e2 + l.rad*sin(u).*sin(v)
z = l.orig.e3 + l.rad*cos(v)
plot_surface(x, y, z, color="blue",alpha=.5,linewidth=0) # and finally ploting the sphere with some transparency so that we'll be able to see the ray

