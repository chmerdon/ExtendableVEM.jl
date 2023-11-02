"""
````
function solve_poisson_serendipity(
	xgrid::ExtendableGrid;
	rhs::Function,
	eu::Function,
	g::Function,
	stab_coeff = 1
	kwargs...)
````

2nd order S-VEM solver for the Poisson problem with right-hand side rhs, exact solution eu
and full Dirichlet boundary data g, which can be functions of the form
(x) -> some value depending on x.

The main assembly loop is inspired by

L. Beirao da Veiga, et al. "Serendipity Nodal VEM spaces." Computers & Fluids (2016):2-12. https://doi.org/10.1016/j.compfluid.2016.02.015

"""
function solve_poisson_serendipity(xgrid::ExtendableGrid{Tc, Ti}; penalty = 1e30, rhs = x -> 2, eu = x -> 1, eu_dx = x -> 1, g = x -> 0, stab_coeff = 1) where {Tc, Ti}
	ndofs_nodes = num_nodes(xgrid)
	ndofs_edges = size(xgrid[FaceNodes],2)
	ncells = num_cells(xgrid)
	ndofs = ndofs_nodes + ndofs_edges
	@info "Assembling S2-VEM (ndofs = $ndofs)..."
	coordinates::Matrix{Tc} = xgrid[Coordinates]
	cellnodes::VariableTargetAdjacency{Ti} = xgrid[CellNodes]
	celldiameters::Vector{Tc} = xgrid[CellDiameters]
	cellcenters::Matrix{Tc} = xgrid[CellCenters]
	celledges = xgrid[CellFaces]
	#edges_normals = xgrid[FaceNormals]
	npolys = 3    # k*(k+1)/2
	npolys_s = 6  # (k+1)*(k+2)/2
	A = ExtendableSparseMatrix{Tc, Int64}(ndofs, ndofs)
	b = zeros(ndofs) # right-hand side vector
	x = zeros(ndofs) # solution vector

	# max_vemdofs = 2*max_num_targets_per_source(cellnodes)
	# D_max = zeros(max_vemdofs, npolys_s)
	# B_max = zeros(npolys_s, 2*max_vemdofs)
	# P_max = zeros(npolys_s, max_vemdofs)
	# S_max = zeros(max_vemdofs, max_vemdofs)
	# Aloc_max = zeros(max_vemdofs, max_vemdofs)
	# ele_node_max = zeros(Int32, max_vemdofs)

	normal1 = zeros(2) # normal vert
	normal2 = zeros(2) # normal edge

	# # exact function
	# eu1 = zeros(ndofs_nodes)
	# for i in 1:ndofs_nodes
	# 	eu1[i] = eu(coordinates[:,i])
	# end
	# eu2 = zeros(ndofs_edges)

	# # storage record matrix
	ele_nodes = Array{Vector,1}([])
	Ps = Array{Matrix,1}([])

	H = zeros(npolys_s, npolys_s)
	G_s = zeros(npolys_s, npolys_s)
	G = zeros(npolys_s, npolys_s)
	H_s = zeros(npolys_s, 1)
	fb = zeros(npolys_s, 1)


	for cell ∈ 1:ncells
		## collect information for polygon
		nvertices = num_targets(cellnodes, cell)
		vertices = view(cellnodes, :, cell)
		nvemdofs = 2*nvertices ## for S-VEM of degree 2
		center = view(cellcenters, :, cell)
		diameter = celldiameters[cell]
		edges = view(celledges, :, cell)
		coordinate = view(coordinates, :, vertices)
		#edge_normals = view(edges_normals, :, edges)
		
		## assembly auxiliary matrices for projector
		## of VEM functions v_j to polynomial basis p_m
		## where m is the monomial index
		D = zeros(nvemdofs, npolys_s)
		B_s = zeros(npolys_s, nvemdofs)
		B = zeros(npolys_s, nvemdofs)
		C = zeros(npolys, 2*nvemdofs)
		P = zeros(npolys_s, nvemdofs)
		P_s = zeros(npolys_s, nvemdofs)
		S = zeros(nvemdofs, nvemdofs)
		Aloc = zeros(nvemdofs, nvemdofs)
		ele_node = zeros(Int32, nvemdofs)
		
		## H[α,β] = (p_α, p_β)
		for i in 1:npolys_s
			for j in 1:npolys_s
				H[i,j] = Gauss_int_2D_tri((x,y)->poly_base(x,y,i,0,center,diameter)*poly_base(x,y,j,0,center,diameter),9,coordinate,center)
			end
		end

		for j ∈ 1:nvertices
			vert = view(coordinates, :, vertices[j])
			prev = view(coordinates, :, vertices[mod1(j - 1, nvertices)])
			next = view(coordinates, :, vertices[mod1(j + 1, nvertices)])
			verc = (vert + next)./2

			#eu2[edges[j]] = eu(verc)

			## calculate sum of normal vectors of two adjacent edges of vert
			## calculate normal vectors of edge
			normal1[1] = next[2] - prev[2]
			normal1[2] = prev[1] - next[1]
			normal2[1] = next[2] - vert[2]
			normal2[2] = vert[1] - next[1]

			for m in 1:npolys_s
				## D[j,m] = j-th dof evaluated for monomial p_m 
				##		  = p_m(ver)
				## ver include vert and verc
				D[j, m] = poly_base(vert[1],vert[2],m,0,center,diameter)
				D[j + nvertices, m] = poly_base(verc[1],verc[2],m,0,center,diameter)
			end

			for m in 1:npolys
				## C[m,[j j + nvemdofs]] = (p_m, u_h)_Γ
				##		                 = p_m ⋅ u_h ⋅ n
				##                       = [p_m ⋅ ver ⋅ n_1 p_m ⋅ ver ⋅ n_2]
				## ver include vert:1/6 and verc:2/3 using Simpson formula
				C[m, [j j + nvemdofs]] = 1/6 * poly_base(vert[1],vert[2],m,0,center,diameter) .* normal1'  #(edge_normals[:,j]' + edge_normals[:,mod1(j - 1, nvertices)]')
				C[m, [j + nvertices j + nvertices + nvemdofs]] = 2/3 * poly_base(verc[1],verc[2],m,0,center,diameter) .* normal2'  #edge_normals[:,j]'
			end
		end
		B .= D'    # B = D'
		## local projector matrix Π_k^S : S-VEM -> P_k
		mul!(G, B, D)
		P .= G \ B
		## local stabilisation matrix
		S .= (I - D * P)
		
		## calculate the (Π_k^S u_h, ∇ ⋅ p_m)
		## p_m ∈ P_1, ∇ ⋅ p_m ∈ P_0
		## (Π_k^S u_h, ∇ ⋅ p_m) = Π_k^S (p_α, p_1)
		for i in 1:npolys_s
			H_s[i] = Gauss_int_2D_tri((x,y)->1/diameter*poly_base(x,y,i,0,center,diameter),9,coordinate,center)
		end
		CC = H_s' * P_s

		## calculate the B[m,j] = - (Π_k^S u_j, ∇ ⋅ p_m) + (p_m, u_j)_Γ
		C[2,1:nvemdofs] -= CC'
		C[3,nvemdofs+1:end] -= CC'

		B_s[1:3,:] .= C[1:3,1:nvemdofs]
		B_s[4:6,:] .= C[1:3,nvemdofs+1:end]

		## calculate the G[i,j] = (p_i, p_j), p_i ∈ (P_k-1)^2	
		G_s[1:npolys,1:npolys] .= H[1:3,1:3]
		G_s[npolys+1:end,npolys+1:end] .= H[1:3,1:3]

		## local projector matrix Π_k-1^0 : ∇u -> P_k-1
		P_s .= G_s \ B_s 

		## local stiffness matrix for S-VEM
		Aloc .= P_s' * G_s * P_s .+ stab_coeff * S' * S

		## record the element dofs number
		ele_node .= [vertices; ndofs_nodes .+ edges] 

		## storage the element dofs number and element projector
		push!(Ps, P)
		push!(ele_nodes, ele_node)

		## add local matrix to global matrix
		for j ∈ 1:nvemdofs
			for k ∈ 1:nvemdofs
				A[ele_node[j], ele_node[k]] += Aloc[j, k]
			end
		end

		## assemble right-hand side
		## (f, v) = (f, Π_k^S v) = Π_k^S (f, p_m)
		for i in 1:npolys_s
			fb[i] = Gauss_int_2D_tri((x,y)->rhs([x y])*poly_base(x,y,i,0,center,diameter),9,coordinate,center)
		end
		view(b, ele_node) .+= P' * fb
	end

	## boundary data by penalization
	boundary_nodes = unique(view(xgrid[BFaceNodes], :))
	bfacenodes = xgrid[BFaceNodes]
	facenodes = xgrid[FaceNodes]
	boundary_edges = findm(bfacenodes,facenodes)
	for node in boundary_nodes
		x[node] = g(view(coordinates, :, node))
		A[node, node] = penalty
		b[node] = x[node] * penalty
		# other method for boundary condition
		# A[node,:] .= 0
		# A[node, node] = 1 
		# b[node] = x[node]
	end
	for boundary_edge in boundary_edges
		x[ndofs_nodes + boundary_edge] = g((view(coordinates, :, facenodes[1,boundary_edge]) .+ view(coordinates, :, facenodes[2,boundary_edge]))./2)
		A[ndofs_nodes + boundary_edge, ndofs_nodes + boundary_edge] = penalty
		b[ndofs_nodes + boundary_edge] = x[ndofs_nodes + boundary_edge] * penalty
		# other method for boundary condition
		# A[ndofs_nodes + boundary_edge,:] .= 0
		# A[ndofs_nodes + boundary_edge, ndofs_nodes + boundary_edge] = 1 
		# b[ndofs_nodes + boundary_edge] = x[ndofs_nodes + boundary_edge]
	end
	flush!(A)

	## solve linear system
	@info "Solving..."
	x = A \ b
	residual = norm(A*x - b)
	@info "...finished with residual = $residual"

	## exact solution 
	# EU =[eu1;eu2]

	## error 
	# energy norm
	#er_E = sqrt.(abs.((x-EU)'*A*(x-EU)))/abs.(x'*A*x);

	er_L2 = 0
	er_H1 = 0
	for cell in 1:ncells
		u = x[ele_nodes[cell]] # element dofs
		P = Ps[cell]           # element projector
		a = P * u              # element numerical solution coefficient
		center = view(cellcenters, :, cell)
		diameter = celldiameters[cell]
		vertices = view(cellnodes, :, cell)
		coordinate = view(coordinates, :, vertices)
		er_L2 += Gauss_int_2D_tri((x,y)->er_f([x y],eu,a,0,center,diameter,2),9,coordinate,center)
		er_H1 += Gauss_int_2D_tri((x,y)->er_f([x y],eu_dx,a,1,center,diameter,2),9,coordinate,center) + Gauss_int_2D_tri((x,y)->er_f([x y],eu_dy,a,2,center,diameter,2),9,coordinate,center)
	end
	er_L2 = sqrt(er_L2)
	er_H1 = sqrt(er_H1)
    @info "...finished with er_L2 = $er_L2"
	@info "...finished with er_H1 = $er_H1"

	return x
end