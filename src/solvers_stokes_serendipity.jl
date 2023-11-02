"""
````
function solve_stokes_serendipity(
	xgrid::ExtendableGrid;
	rhs::Function,
	eu::Function,
	g::Function,
	stab_coeff = 1
	kwargs...)
````

2nd order S-VEM solver for the Stokes problem with right-hand side rhs, exact solution eu
and full Dirichlet boundary data g, which can be functions of the form
(x) -> some value depending on x.

The main assembly loop is inspired by

L. Beirao da Veiga, et al. "Serendipity Nodal VEM spaces." Computers & Fluids (2016):2-12. https://doi.org/10.1016/j.compfluid.2016.02.015
P. Wriggers, M. L. D. Bellis and B. Hudobivnik . "A Taylor-Hood type virtual element formulations for large incompressible strains." Computer Methods in Applied Mechanics and Engineering 385.4(2021):114021. https://doi.org/10.1016/j.cma.2021.114021

"""
function solve_stokes_serendipity(xgrid::ExtendableGrid{Tc, Ti}; penalty = 1e30, rhs = x -> 2, eu = x -> 1, eu_dx = x -> 1, eu_dy = x -> 1, g = x -> 0, stab_coeff = 1) where {Tc, Ti}
	ndofs_nodes = num_nodes(xgrid)
	ndofs_edges = size(xgrid[FaceNodes],2)
	ncells = num_cells(xgrid)
	ndofs_u = 2 * (ndofs_nodes + ndofs_edges)
	ndofs_u_x = ndofs_nodes + ndofs_edges
	ndofs_p = ndofs_nodes
	ndofs = ndofs_u + ndofs_p
	@info "Assembling S2-S1-VEM (ndofs = $ndofs)..."
	coordinates::Matrix{Tc} = xgrid[Coordinates]
	cellnodes::VariableTargetAdjacency{Ti} = xgrid[CellNodes]
	celldiameters::Vector{Tc} = xgrid[CellDiameters]
	cellcenters::Matrix{Tc} = xgrid[CellCenters]
	celledges = xgrid[CellFaces]
	npolys = 3      # k*(k+1)/2
	npolys_s = 6	# (k+1)*(k+2)/2
	A = ExtendableSparseMatrix{Tc, Int64}(ndofs + 1, ndofs + 1) # add one dof to meet ∫pdx = 0 
	b = zeros(ndofs + 1) # right-hand side vector
	xx = zeros(ndofs + 1) # solution vector
	x = zeros(ndofs) # solution vector

	# max_vemdofs = 5*max_num_targets_per_source(cellnodes) + 1
	# D_max = zeros(max_vemdofs, npolys_s)
	# B_max = zeros(npolys_s, max_vemdofs)
	# P_max = zeros(2*npolys_s, max_vemdofs)
	# S_max = zeros(max_vemdofs, max_vemdofs)
	# Aloc_max = zeros(max_vemdofs, max_vemdofs)

	# exact function
	# eu1 = zeros(ndofs_nodes)
	# for i in 1:ndofs_nodes
	# 	eu1[i] = eu(coordinates[:,i],1)
	# end
	# eu2 = zeros(ndofs_edges)
	# eu3 = zeros(ndofs_nodes)
	# for i in 1:ndofs_nodes
	# 	eu1[i] = eu(coordinates[:,i],2)
	# end
	# eu4 = zeros(ndofs_edges)
	# eu5 = zeros(ndofs_nodes)
	# for i in 1:ndofs_nodes
	# 	eu1[i] = eu(coordinates[:,i],3)
	# end

	normal1 = zeros(2) # normal vert
	normal2 = zeros(2) # normal edge

	# # storage record matrix
	ele_nodes = Array{Vector,1}([])
	Pus = Array{Matrix,1}([])
	Pps = Array{Matrix,1}([])

	H = zeros(npolys_s, npolys_s)
	G = zeros(npolys_s, npolys_s)
	G_p = zeros(npolys, npolys)
	G_s = zeros(npolys_s, npolys_s)
	H_s = zeros(npolys_s, 1)
	fb = zeros(2*npolys_s, 1)

	for cell ∈ 1:ncells
		## collect information for polygon
		nvertices = num_targets(cellnodes, cell)
		vertices = view(cellnodes, :, cell)
		nvemdofs_u = 4*nvertices
		nvemdofs_u_x = 2*nvertices
		nvemdofs_p = nvertices
		nvemdofs = nvemdofs_u + nvemdofs_p ## for S2-S1-VEM
		center = view(cellcenters, :, cell)
		diameter = celldiameters[cell]
		edges = view(celledges, :, cell)
		coordinate = view(coordinates, :, vertices)

		## assembly auxiliary matrices for projector
		## of VEM functions v_j to polynomial basis p_m
		## where m is the monomial index
		D = zeros(nvemdofs_u_x, npolys_s)
		D_p = zeros(nvemdofs_p, npolys)
		B_p = zeros(npolys, nvemdofs_p)
		B_s = zeros(npolys_s, nvemdofs_u_x)
		B = zeros(npolys_s, nvemdofs_u_x)
		P = zeros(npolys_s, nvemdofs_u_x)
		P_p = zeros(npolys, nvemdofs_p)
		P_s = zeros(npolys_s, nvemdofs_u_x)
		P_u = zeros(2*npolys_s, nvemdofs_u)
		S = zeros(nvemdofs_u_x, nvemdofs_u_x)
		Aloc = zeros(nvemdofs_u, nvemdofs_u)
		Aloc_u = zeros(nvemdofs_u_x, nvemdofs_u_x)
		Bloc = zeros(nvemdofs_u, nvemdofs_p)
		C = zeros(npolys, nvemdofs_u)
		P_div = zeros(npolys, nvemdofs_u)
		PP = zeros(2*npolys_s, nvemdofs_u)
		ele_node = zeros(Int64,nvemdofs)
				
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

			# eu2[edges[j]] = eu(verc,1)
			# eu4[edges[j]] = eu(verc,2)

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
				C[m, [j j + nvemdofs_u_x]] = 1/6*poly_base(vert[1],vert[2],m,0,center,diameter) .* normal1'
				C[m,[j + nvertices j + nvertices + nvemdofs_u_x]] = 2/3*poly_base(verc[1],verc[2],m,0,center,diameter) .* normal2'
			end
		end
		B .= D'  # B = D'
		D_p = D[1:nvertices, 1:npolys]  # D_p = D [1:nvertices, 1:npolys]
		B_p .= D_p' # B_p = D_p'

		## local projector matrix Π_k^S : S-VEM -> P_k
		mul!(G, B, D)
		mul!(G_p, B_p, D_p)
		P .= G \ B
		P_p .= G_p \ B_p

		## local stabilisation matrix
		S .= (I - D * P)

		## local stiffness matrix for VEM
		## calculate the (Π_k^S u_h, ∇ ⋅ p_m) and (Π_k^S u_h, ∇ p_m)
		## p_m ∈ P_1, ∇ ⋅ p_m ∈ P_0, ∇ p_m ∈ (P_0)^2
		## (Π_k^S u_h, ∇ ⋅ p_m) = Π_k^S (p_α, p_1)
		## (Π_k^S u_h, ∇ p_m) = [Π_k^S (p_α, p_1) Π_k^S (p_α, p_1)]
		for i in 1:npolys_s
			H_s[i] = Gauss_int_2D_tri((x,y)->1/diameter*poly_base(x,y,i,0,center,diameter),9,coordinate,center)
		end
		CC = H_s' * P_s

		## calculate the B[m,j] = - (Π_k^S u_j, ∇ ⋅ p_m) + (p_m, u_j)_Γ
		C[2, 1:nvemdofs_u_x] .-= CC'
		C[3, nvemdofs_u_x + 1:end] .-= CC'

		B_s[1:npolys,:] .= C[1:3,1:nvemdofs_u_x]
		B_s[npolys + 1:npolys_s,:] .= C[1:3,nvemdofs_u_x+1:end]

		## local projector matrix Π_k-1^div : ∇⋅u -> P_k-1
		H_p = H[1:npolys, 1:npolys]

		P_div .= H_p \ C
		
		## local stiffness matrix (divu,p) for S-VEM 
		Bloc .= P_div' * H_p * P_p 

		## calculate the G[i,j] = (p_i, p_j), p_i ∈ (P_k-1)^2
		G_s[1:npolys,1:npolys] .= H_p
		G_s[npolys+1:end,npolys+1:end] .= H_p

		## local projector matrix Π_k-1^0 : ∇u -> P_k-1
		P_s .= G_s \ B_s

		## local stiffness matrix u for S-VEM 
		Aloc_u .= P_s' * G_s * P_s .+ stab_coeff * S' * S

		## local stiffness matrix [u_x u_y]
		Aloc[1:nvemdofs_u_x, 1:nvemdofs_u_x] .= Aloc_u
		Aloc[nvemdofs_u_x + 1:end, nvemdofs_u_x + 1:end] .= Aloc_u
		
		loc = [Aloc Bloc;Bloc' zeros(nvemdofs_p,nvemdofs_p)]

		## record the element dofs number
		ele_node .= [vertices; ndofs_nodes .+ edges; ndofs_u_x .+ vertices; ndofs_u_x + ndofs_nodes .+ edges; ndofs_u .+ vertices]

		## storage the element dofs number and element projector
		PP[1:npolys_s, 1:nvemdofs_u_x] .= P 
		PP[npolys_s + 1 : end, nvemdofs_u_x + 1 : end] .= P 
		push!(Pus, PP)
		push!(Pps, P_p)
		push!(ele_nodes, ele_node)

		## add local matrix to global matrix
		for j ∈ 1:nvemdofs
			for k ∈ 1:nvemdofs
				A[ele_node[j], ele_node[k]] += loc[j, k]
			end
		end

		# calculate ∫pdx = 0 
		L = H[1,1:npolys]

		L = P_p' * L

		## add L to global matrix
		for i ∈ 1:nvertices
			A[end, ndofs_u + vertices[i]] += L[i]
			A[ndofs_u + vertices[i], end] += L[i]
		end

		## record the [u_x, u_y] projector
		P_u[1:npolys_s, 1:nvemdofs_u_x] .= P
		P_u[1 + npolys_s:end, 1 + nvemdofs_u_x:end] .= P

		## assemble right-hand side
		## (f, v) = (f, Π_k^S v) = Π_k^S (f, p_m)
		for i in 1:npolys_s
			fb[i] = Gauss_int_2D_tri((x,y)->rhs([x y],1)*poly_base(x,y,i,0,center,diameter),9,coordinate,center)
			fb[i + npolys_s] = Gauss_int_2D_tri((x,y)->rhs([x y],2)*poly_base(x,y,i,0,center,diameter),9,coordinate,center)
		end
		view(b, ele_node[1:nvemdofs_u]) .+= P_u' * fb
	end

	## boundary data by penalization
	boundary_nodes = unique(view(xgrid[BFaceNodes], :))
	bfacenodes = xgrid[BFaceNodes]
	facenodes = xgrid[FaceNodes]
	boundary_edges = findm(bfacenodes,facenodes)
	for node in boundary_nodes
		x[node] = g(view(coordinates, :, node),1)
		A[node, node] = penalty
		b[node] = x[node] * penalty
		# other method for boundary condition
		# A[node,:] .= 0
		# A[node, node] = 1 #penalty
		# b[node] = x[node] #* penalty

		x[node + ndofs_u_x] = g(view(coordinates, :, node),2)
		A[node + ndofs_u_x, node + ndofs_u_x] = penalty
		b[node + ndofs_u_x] = x[node + ndofs_u_x] * penalty
		# other method for boundary condition
		# A[node + ndofs_u_x,:] .= 0
		# A[node + ndofs_u_x, node + ndofs_u_x] = 1 #penalty
		# b[node + ndofs_u_x] = x[node + ndofs_u_x]  #* penalty
	end
	for boundary_edge in boundary_edges
		x[ndofs_nodes + boundary_edge] = g((view(coordinates, :, facenodes[1,boundary_edge]) .+ view(coordinates, :, facenodes[2,boundary_edge]))./2,1)
		A[ndofs_nodes + boundary_edge, ndofs_nodes + boundary_edge] = penalty
		b[ndofs_nodes + boundary_edge] = x[ndofs_nodes + boundary_edge] * penalty
		# other method for boundary condition
		# A[ndofs_nodes + boundary_edge,:] .= 0
		# A[ndofs_nodes + boundary_edge, ndofs_nodes + boundary_edge] = 1 
		# b[ndofs_nodes + boundary_edge] = x[ndofs_nodes + boundary_edge]

		
		x[ndofs_u_x + ndofs_nodes + boundary_edge] = g((view(coordinates, :, facenodes[1,boundary_edge]) .+ view(coordinates, :, facenodes[2,boundary_edge]))./2,2)
		A[ndofs_u_x + ndofs_nodes + boundary_edge, ndofs_u_x + ndofs_nodes + boundary_edge] = penalty
		b[ndofs_u_x + ndofs_nodes + boundary_edge] = x[ndofs_u_x + ndofs_nodes + boundary_edge] * penalty
		# other method for boundary condition
		# A[ndofs_u_x + ndofs_nodes + boundary_edge,:] .= 0
		# A[ndofs_u_x + ndofs_nodes + boundary_edge, ndofs_u_x + ndofs_nodes + boundary_edge] = 1 
		# b[ndofs_u_x + ndofs_nodes + boundary_edge] = x[ndofs_u_x + ndofs_nodes + boundary_edge]
	end
	flush!(A)

	## solve linear system
	@info "Solving..."
	xx = A \ b
	residual = norm(A*xx - b)
	@info "...finished with residual = $residual"

	x = xx[1:end-1]

	## exact solution 
	# EU =[eu1;eu2;eu3;eu4;eu5]
	
	## error 
	# energy norm
	# er_E = sqrt.(abs.((x-EU)'*A[1:end-1,1:end-1]*(x-EU)))/abs.(x'*A[1:end-1,1:end-1]*x);
    # @info "...finished with er_E = $er_E"

	er_L2 = 0
	er_H1 = 0
	er_L2_p = 0
	for cell in 1:ncells
		nvertices = num_targets(cellnodes, cell)
		center = view(cellcenters, :, cell)
		diameter = celldiameters[cell]
		vertices = view(cellnodes, :, cell)
		coordinate = view(coordinates, :, vertices)
		nvemdofs_u = 4 * nvertices;
		u = x[ele_nodes[cell]]								# element dofs
		Pu = Pus[cell]										# element u projector
		Pp = Pps[cell]										# element p projector
		au = Pu * u[1:nvemdofs_u]                          # element u numerical solution coefficient
		ap = Pp * u[nvemdofs_u + 1:end]                    # element p numerical solution coefficient
		er_L2 += Gauss_int_2D_tri((x,y)->er_f_2([x y],eu,au,0,center,diameter,2),9,coordinate,center)
		er_H1 += Gauss_int_2D_tri((x,y)->er_f_2([x y],eu_dx,au,1,center,diameter,2),9,coordinate,center) + Gauss_int_2D_tri((x,y)->er_f_2([x y],eu_dy,au,2,center,diameter,2),9,coordinate,center)
		er_L2_p += Gauss_int_2D_tri((x,y)->er_fp([x y],eu,ap,0,center,diameter,2),9,coordinate,center)
	end
	er_L2 = sqrt(er_L2)
	er_H1 = sqrt(er_H1)
	er_L2_p = sqrt(er_L2_p)
    @info "...finished with er_L2 = $er_L2"
	@info "...finished with er_H1 = $er_H1"
	@info "...finished with er_L2_p = $er_L2_p"

	return x
end