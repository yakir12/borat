module Tetrahedrons

using ImmutableArrays

export Tetrahedron
export isInside

# helper, computes the determinant of 4 vector matrix
function det4Vm(v1::Vector3,v2::Vector3,v3::Vector3,v4::Vector3)
	m = [[v1 v2 v3 v4]' [1.,1.,1.,1.]]
	det(m)
end
 

immutable Tetrahedron
	# the 4 points define the tetrahedron
	corners::Array{Vector3{Float64},1}
	
	# determinant 0
	d0::Float64
	
	# create determinant d0
	function Tetrahedron(va::Array{Vector3,1}) 
		if length(va) != 4
			error("Four points needed to describe a tetrahedron.")
		end
		d0 = det4Vm(va[1],va[2],va[3],va[4])
		
		# this comparison might need some error margin
		if d0 == 0.
			error("Degenerated tetrahedron (all corners coplanar).") 
		end
		new(va,d0)
	end
	
end

# a constructor
Tetrahedron(v1::Vector3,v2::Vector3,v3::Vector3,v4::Vector3) = Tetrahedron(Vector3[v1,v2,v3,v4])
# t = Tetrahedron(Vector3([0.,0.,0.]),Vector3([1.,0.,0.]),Vector3([0.,1.,0.]),Vector3([0.,0.,1.]))

# test for a point inside a tetrahedron
function isInside(p::Vector3, t::Tetrahedron):Bool

	sign(t.d0) == sign(det4Vm(p,t.corners[2],t.corners[3],t.corners[4])) && 	
	sign(t.d0) == sign(det4Vm(t.corners[1],p,t.corners[3],t.corners[4])) &&
	sign(t.d0) == sign(det4Vm(t.corners[1],t.corners[2],p,t.corners[4])) &&
	sign(t.d0) == sign(det4Vm(t.corners[1],t.corners[2],t.corners[3],p))	
end

end
