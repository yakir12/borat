module LightRay

using ImmutableArrays

export Ray
export initiate, refract!

# the light ray type
type Ray 
	# has a (current) position, a location
    pos::Vector3 
    
    # and a direction, which has to be unitized, norm(dir) = 1
    dir::Vector3 

    # see 1.13.2 Inner Constructor Methods. I'm tempted to just do this:
    # Ray(pos,dir) = new(pos,unit(dir))
    # but they say it's bad form???
    Ray(pos,dir) = norm(dir) != 1. ? error("directional vector not unitized") : new(pos,dir) 
end

# some methods associated with rays

# This function initiates one ray. Each starts from a given position rOrg,  
# and is randomly directed through a circle of radius rad. The circle is assumed to have
# the center exactly under the ray origin, at distance L on the vertical.
# The distribution is uniform.
# NOTE: this is still very much focused on spherical lenses - could not change this
#       completely, but it must be done, to make it more general
function initiate(rOrg::Vector3, L::Real, rad::Real): Ray
    if isinf(L) # if the light source is infinitely far away
        #=just ditribute the rays uniformly on the disk that the sphere makes=#
        sqrtr = sqrt(rand()) 
        θ = rand()*2π
        x = rOrg[1] + sqrtr * cos(θ)
        y = rOrg[2] + sqrtr * sin(θ)
        # the ray actually points vertically down, 
        # starting at the computed point on the disc
        Ray(Vector3([x,y,2rad]),Vector3([0.,0.,-1.]))
    else
        # theta is the angle between the optical axis and 
        # the ray that touches the surface of the sphere
        θ = asin(rad/L)

        # epsilon is used only for the range in the rand function below, you'll see
        ε = eps(typeof(θ))

        # z and h are randomly sampled from those ranges to generate the uniformly spaced 
        # points on the sphere
        z = rand(cos(θ):ε:1.)
        h = rand(0.:ε:2π)

        # now come the x, y, and z that lie on the sphere's surface
        x = sqrt(1. - z.^2).*cos(h)
        y = sqrt(1. - z.^2).*sin(h)
        
        # I flip the z so it'll point downwards
        z *= -1. 
		
		# now build the ray that starts at origin and is directed towards
		# the point computed previously
        Ray(rOrg,unit(Vector3([x,y,z])))
    end
 
end



# Refraction occurs at the surface between the medium and the lens, not only in the lens. 
# So this function can be used in both scenarios. 
# r is the ray, N is the normal to the refractive surface, and n is the ratio between 
# the refractive index before the surface and the refractive index after the surface.
function refract!(r::Ray,N::Vector3,n::Real)
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

end
