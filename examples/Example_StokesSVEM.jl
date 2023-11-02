module Example_StokesSVEM

using ExtendableVEM
using ExtendableGrids
using GridVisualize

### exact solution eu, right-hand side rhs and full Dirichlet boundary data g
function eu(x,k) # k:1 u_x, k:2 u_y, k:3 p
	if k == 1
		eu = x[1].^2 .* x[2].^2 + exp(-x[2])
	elseif k == 2
		eu = -2/3 .* x[1] .* x[2].^3 .+ 2 .- pi * sin(π .* x[1])
	elseif k == 3
		eu = (2 .- pi .* sin(π .* x[1])) .* cos( 2 .* π .* x[2])
	end
	return eu
end

function eu_dx(x,k) # k:1 u_x, k:2 u_y, k:3 p
	if k == 1
		eu = 2 * x[1] .* x[2].^2
	elseif k == 2
		eu = -2/3 .* x[2].^3 .+ 2 .- pi^2 * cos(π .* x[1])
	elseif k == 3
		eu =  pi^2 .* cos(π .* x[1]) .* cos( 2 .* π .* x[2])
	end
	return eu
end

function eu_dy(x,k) # k:1 u_x, k:2 u_y, k:3 p
	if k == 1
		eu = 2 * x[1].^2 .* x[2] - exp(-x[2])
	elseif k == 2
		eu = -2 .* x[1] .* x[2].^2
	elseif k == 3
		eu = 2 .* π .* (2 .- pi .* sin(π .* x[1])) .* sin( 2 .* π .* x[2])
	end
	return eu
end

function rhs(x,k) # k:1 f_x, k:2 f_y
	if k == 1
		rhs  = -2 * x[1].^ 2 .- 2 * x[2] .^ 2 - exp( - x[2]) + pi^2 * cos( π * x[1]) .* cos( 2 .* π * x[2])
	elseif k == 2
		rhs  = 4 * x[1] .* x[2] - pi^3 .* sin(π * x[1]) + 2 * pi * (2 .- pi .* sin(π .* x[1])) .* sin( 2 * π * x[2])
	end
	return rhs
end

g(x,k) = eu(x,k)

function main(; nrefs = 3, deform = 0.2 * 2.0^(-nrefs), Plotter = nothing)
	## load grid with triangles and rectangles and frame them as polygons
	xgrid = uniform_refine(grid_unitsquare_mixedgeometries(), nrefs)
	xgrid[CellGeometries] = VectorOfConstants{ElementGeometries, Int}(Polygon2D, num_cells(xgrid))

	## deform coordinates a little to make it interesting
	nnodes = num_nodes(xgrid)
	coordinates = xgrid[Coordinates]
	for j ∈ 1:nnodes
		if any(abs.(coordinates[:, j]) .<= 0) || any(abs.(coordinates[:, j]) .>= 1)
		else
			coordinates[:, j] .+= deform * (rand() - 0.5)
		end
	end

	## calculate facenodes, cellcenters, celldiameters, cellvolumes, subtriangulation
	prepare_polygons!(xgrid)

	## solve VEM
	x = solve_stokes_serendipity(xgrid; rhs = rhs, eu = eu, eu_dx = eu_dx, eu_dy = eu_dy, g = g)

	# ## calculate average values at polygon centers for plot on subtriangulation
	# subtriangulation = xgrid[SubTriangulation]
	# ncells = num_cells(xgrid)
	# append!(x, zeros(Float64, ncells))
	# cellnodes = xgrid[CellNodes]
	# for cell ∈ 1:num_cells(xgrid)
	# 	x[nnodes+cell] = sum(view(x, view(cellnodes, :, cell))) / num_targets(cellnodes, cell)
	# end

	# ## plot polygons
	# f = plot_polygons(xgrid)
	# #matplotlib.use("TKAgg")

	# ## plot subtriangulations and solution
	# pl = GridVisualizer(; Plotter = Plotter, layout = (1, 2), clear = true, size = (800, 400))
	# gridplot!(pl[1, 1], subtriangulation, linewidth = 1)
	# scalarplot!(pl[1, 2], subtriangulation, x; Plotter = Plotter)
	# #reveal(pl)

	return x
end

end #module