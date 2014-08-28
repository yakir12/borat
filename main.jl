using ImmutableArrays, PyPlot # useful library for fast matrix manipulations, ploting in 3D

#=first we define a few useful types. Since Julia is a highly (but not stricktly) typed language, this is an awesome feature that speeds up the code and makes it clean and easy to follow=#
type Ray # the ray type
    pos::Vector3 # has a position, a location
    dir::Vector3 # and a direction, which has to be unitized, norm(dir) = 1
end

type Lens # the lens type, we restrict ourselves to sphere shaped lenses
    org::Vector3 # the location of the sphere center
    rad::Number # the radius of the sphere
    rig::Function # the Refractive Index Gradient (RIG) that describes the distribution of refractive indices inside the lens. We restrict ourdelves to spherically symmetrical distributions: so at a set distance form the lens center, the refractive index is equal in all directions (kind of like an onion?).
end
#=We should probably define a few more helpful types in the near future...=#

#=Now come the constants, there is no real reason to define the following variables as constants, other than that it makes sense that these don't change in value throughout the program. =#
const step = 1e-2 # the step size withwhich the ray advances every iterration
const rad = 1. # the lens radius, this is only usfule cause there are a couple of calculations that depend on htis variable
const n_periphery = 1.34 # the refractive index at the periphery of the lens, so at the very surface of the lens, when the distance from the lens center is equal to rad
const n_center = 1.45 # the refractive index at the center of the lens, when the distance from the lens center is equal to zero
const L = 3.*rad # this is the distance between the center of the lens and the source light
const nrays = 10 # the number of discrete rays I'll be tracing
const c = Vector3([0.,0.,0.]) # this is the location of the center of the lens

#=Now come the functions=#
#=this is the RIG function. Notice that its argument is now Number -- this is in the interest of generality... I set this RIG to be equal to the Luneburg lens, so that we can check the results real quick. see: http://en.wikipedia.org/wiki/Luneburg_lens=#
rig(r::Number) = sqrt(2. - (r/rad)^2)

#=this function initiates the array of Rays. It starts them off from directly above the lens at distance L. Thier directions are uniformly and randomly distributed across the lens.=#
function initiate(l::Lens,L::Number,n::Int)
    #=theta is the angle between the optical axis and the ray that touches the surface of the sphere=#
    θ = asin(l.rad/L)

    #=epsilon is used only for the range in the rand function below, you'll see=#
    ε = eps(typeof(θ))

    #=z and h are randomly sampled from those ranges to generate the uniformaly spaced points on the sphere=#
    z = rand(cos(θ):ε:1.,n)
    h = rand(0.:ε:2π,n)

    #=now come the x, y, and z that lie on the sphere's surface=#
    x = sqrt(1. - z.^2).*cos(h)
    y = sqrt(1. - z.^2).*sin(h)
    z *= -1. # I flip the z so it'll point downwards

    #=this is the starting position of all the rays -- point source=#
    s = Vector3([0.,0.,L] .+ l.org)

    #=an empty container for all the Rays=#
    r = Array(Ray,n)
    for i = 1:n

        #=initiate the Ray with its initial position, s, and it's unitized direction that will end up on the sphere=#
        r[i] = Ray(s,unit(Vector3([x[i],y[i],z[i]])))
    end
    return r
end

#=This function finds the intersection point of the light as it starts from its enitial point. It's a basic vector-sphere intersection thing. It works. Sorry I didn't comment it :P =#
function hitlens!(r::Ray,l::Lens)
    a = r.pos - l.org
    b = dot(r.dir,a)^2 - norm(a)^2 + l.rad
    c = -dot(r.dir,a) - sqrt(b)
    r.pos += c*r.dir
end

# This is the main function that advances the ray along its path inside the lens. 
# Everything happens here, and in fact too much does. It would be better if we'd brake it 
# apart into logical segments. The algorithm is taken directly from page 74 in the 
# S1-rt.pdf file I got from http://fileadmin.cs.lth.se/cs/Education/EDAN30/lectures/S1-rt.pdf 
# in that website your graphics guy gave you 
# (i.e. http://cs.lth.se/english/course/edan30-photorealistic-computer-graphics/assignments/). 
# I'll refer to that page during the comments below. Notice that I didn't shifted 
# the system to the lens center as I should. Since it's at 0,0 it doens't matter, 
# but this should be included.

function advance!(r::Ray,l::Lens) 
    # here I find out if the ray is on it's way in or way out. 
    # If its distance to the lens center is diminishing then the normal to the refractive 
    # surface is going outwards. but if the ray is leaving the lens the normal to 
    # the surface is pointing inwards.
    costheta = dot(r.dir,r.pos)
    
    # is it going out? if so make out equal to +1, otherwise make it -1
    out = -sign(costheta)
    
    # this is the normal to the refractive surface 
    N = out*unit(r.pos)
    
    # this is the 'r' in page 74 in that pdf I just mentioned
    a = -dot(r.dir,N)
  
    # save the old dir
    olddir = r.dir
    
    # this is the refractive index before the current position of the ray
    # n1 = l.rig(norm(r.pos + out*(step/2)*unit(r.pos)))
    # why not just at pos?
    n1 = l.rig(norm(r.pos))
    # and this is the refractive index after. this is problematic...
    # n2 = l.rig(norm(r.pos - out*(step/2)*unit(r.pos)))
    # why not just at the new position, after the step?
    n2 = l.rig(norm(r.pos + step*olddir))
    
    # this is the ratio between the two and the 'η' in page 74
    n = n1/n2
    
    # this is the 'c' in page 74
    c = 1. - n^2*(1. - a^2)
    
    # and finally the 't' in page 74
    # added reflection as well... see page 54 (v is -r.dir)
    t = c>0 ? n*r.dir + (n*a - sqrt(c))*N : 2* dot(N,-r.dir) *N + r.dir
        
    # I unitize the direction vector for the ray, updatiung the ray's direction
    r.dir = unit(t)
    
    # and step the ray forward along said direction, updatiung the ray's position
    # r.pos += (step/2)*olddir + (step/2)*r.dir
    r.pos += step*olddir
end

#=OK, now comes the actual calculations=#
l = Lens(c,rad,rig) # initiate the lens, this is how we initiate an instance of a type, it's like a function!
r0 = initiate(l,L,nrays) # initiate all the rays

p = {Array(Vector3,1) for i in 1:nrays} # a cell array of arrays of positions. Defined as an array who's elements are Vector3

#=Now we itirate through all the rays, tracing them=#
for i = 1:nrays
    
    #=I first copy the i^th ray's initial position to r, for conviniance=#
    r = deepcopy(r0[i])
    p[i][1] = r.pos # populate with the start position

    #=find the intersection point and update r's position=#
    hitlens!(r,l)
    push!(p[i],r.pos) # add said position to the position arrays
    advance!(r,l) # advance once, just to penetrate the lens and be able to start that while loop below
    push!(p[i],r.pos) # add said position to the position arrays
    while (norm(r.pos) < l.rad) & (length(p[i]) < 4e3) # loop until the distance of the ray form the lens center is larger or equal to rad OR until you looped more than 200 times (totally arbitrary)
    #    try
            advance!(r,l) # advance one step
    #    catch
    #    	println("got some domain error... we have only this many points:",length(p))
    #    	break
    #    end
        push!(p[i],r.pos) # save the new position in the positions array
    end
end

#=plot all the ray's trajectories=#
for i = 1:nrays
      xyz = cat(2,p[i]...) # here I use the ... operator to extract all the vectors from the locations and concatinate them along the second dimension
    x = vec(xyz[1,:]) # I convert each to a vector (in Julia, a vector is a stricktly one dimentional vector, which is different than a matirx with dimensions [1,2] or [2,1])
    y = vec(xyz[2,:])
    z = vec(xyz[3,:])

    plot3D(x,y,z,color="red") # plot the ray trajectory in red
end

#=now comes the shitty sphere plot. I think I improved it a bit :)=# 
n = 100 # number of points to plot
u = linspace(0., π, n) # the two angles that define the location of the points on the sphere
v = linspace(0., 2π, n)' # notice that I transpose this one with '
u .+= 0*v # the transposition is "carried" and expanded through the .+ operator to make u a n by n matrix
v .+= 0*u # same here
x = l.org.e1 + l.rad*cos(u).*sin(v) # then standard vector math
y = l.org.e2 + l.rad*sin(u).*sin(v)
z = l.org.e3 + l.rad*cos(v)
plot_surface(x, y, z, color="blue",alpha=.5,linewidth=0,rstride=4,cstride=4) # and finally ploting the sphere with some transparency so that we'll be able to see the ray

