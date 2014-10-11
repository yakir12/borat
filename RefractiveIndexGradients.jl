module RefractiveIndexGradients

using Vector3D

export RefractiveIndexGradient, RadialRIG, getRefractiveIndex

abstract RefractiveIndexGradient

type RadialRIG <: RefractiveIndexGradient
	# the Refractive Index Gradient (RIG) that describes the distribution of refractive 
	# indices inside the lens. We restrict ourselves to spherically 
	# symmetrical distributions: so at a set distance form the lens center, 
	# the refractive index is equal in all directions (kind of like an onion?).
	
	# example: ri(r::Float64) = sqrt(2. - r*r/lens_r2)
	# the position needs to be adjusted to the center of symmetry
	rig_rel::Function
	
	# RIG center of symmetry (if the same as lens center, then the lens is symmetrical)
    rig_c::Vec
end

function getRefractiveIndex(rr::RadialRIG, p::Vec)
 	rr.rig_rel(norm(p-rr.rig_c))
end

end