"""
````
function solve_poisson(
	xgrid::ExtendableGrid;
	rhs::Function,
	g::Function,
	stab_coeff = 1
	kwargs...)
````

Lowest order VEM solver for the Poisson problem with right-hand side rhs
and full Dirichlet boundary data g, which can be functions of the form
(x) -> some value depending on x.

The main assembly loop is inspired by

Sutton, O.J. The virtual element method in 50 lines of MATLAB. Numer Algor 75, 1141–1159 (2017). https://doi.org/10.1007/s11075-016-0235-3

"""
function solve_poisson(xgrid::ExtendableGrid{Tc, Ti}; penalty = 1e30, rhs = x -> 1, g = x -> 0, stab_coeff = 1) where {Tc, Ti}
	ndofs = num_nodes(xgrid)
	ncells = num_cells(xgrid)
	coordinates::Matrix{Tc} = xgrid[Coordinates]
	cellnodes::VariableTargetAdjacency{Ti} = xgrid[CellNodes]
	celldiameters::Vector{Tc} = xgrid[CellDiameters]
	cellvolumes::Vector{Tc} = xgrid[CellVolumes]
	cellcenters::Matrix{Tc} = xgrid[CellCenters]
	monomials = [(0, 0), (1, 0), (0, 1)] # monomial multi-indices for VEM of degree 1
	npolys = length(monomials)
	A = ExtendableSparseMatrix{Tc, Int64}(ndofs, ndofs)
	b = zeros(ndofs) # right-hand side vector
	x = zeros(ndofs) # solution vector
	max_vemdofs = max_num_targets_per_source(cellnodes)
	D_max = zeros(max_vemdofs, npolys)
	B_max = zeros(npolys, max_vemdofs)
	P_max = zeros(npolys, max_vemdofs)
	S_max = zeros(max_vemdofs, max_vemdofs)
	Aloc_max = zeros(max_vemdofs, max_vemdofs)
	G = zeros(npolys, npolys)
	normal = zeros(2)
	for cell ∈ 1:ncells
		## collect information for polygon
		nvertices = num_targets(cellnodes, cell)
		vertices = view(cellnodes, :, cell)
		nvemdofs = nvertices ## for VEM of degree 1
		center = view(cellcenters, :, cell)

		## assembly auxiliary matrices for projector
		D = view(D_max, 1:nvemdofs, :)
		B = view(B_max, :, 1:nvemdofs)
		P = view(P_max, :, 1:nvemdofs)
		S = view(S_max, 1:nvemdofs, 1:nvemdofs)
		Aloc = view(Aloc_max, 1:nvemdofs, 1:nvemdofs)
		fill!(D, 0)
		fill!(B, 0)
		D[:, 1] .= 1
		B[1, :] .= 1 / nvemdofs
		for j ∈ 1:nvemdofs
			vert = view(coordinates, :, vertices[j])
			prev = view(coordinates, :, vertices[mod1(j - 1, nvemdofs)])
			nextv = view(coordinates, :, vertices[mod1(j + 1, nvemdofs)])
			normal[1] = nextv[2] - prev[2]
			normal[2] = prev[1] - nextv[1]
			for m in 2:npolys
				poly_degree = monomials[m]
				monomial_grad = poly_degree ./ celldiameters[cell]
				D[j, m] = dot(vert, monomial_grad) - dot(center, monomial_grad)
				B[m, j] = 0.5 * dot(monomial_grad, normal)
			end
		end
		## local projector matrix Π∇ : VEM -> monomials
		mul!(G, B, D)
		P .= G \ B
		## local stabilisation matrix
		S .= (I - D * P)
		## local stiffness matrix for VEM
		G[1, :] .= 0 # G is now local stiffness matrix for monomials
		Aloc .= P' * G * P .+ stab_coeff * S' * S

		## add local matrix to global matrix
		for j ∈ 1:nvertices
			for k ∈ 1:nvertices
				A[vertices[j], vertices[k]] += Aloc[j, k]
			end
		end

		## assemble right-hand side
		view(b, vertices) .+= rhs(center) * cellvolumes[cell] / nvemdofs
	end

	## boundary data by penalization
	boundary_nodes = unique(view(xgrid[BFaceNodes], :))
	for node in boundary_nodes
		x[node] = g(view(coordinates, :, node))
		A[node, node] = penalty
		b[node] = x[node] * penalty
	end
	flush!(A)

	## solve linear system
	@time x = A \ b

	return x
end
