module Example_PoissonSVEM

using ExtendableVEM
using ExtendableGrids
using GridVisualize

### right-hand side rhs, exact solution eu eu_dx eu_dy and full Dirichlet boundary data g
rhs(x) = pi^2 * x[2]^2 * sin(π * x[1]) - 2 * sin(π * x[1])
eu(x) = x[2]^2 * sin(π * x[1])
eu_dx(x) = x[2]^2 * pi * cos(π * x[1])
eu_dy(x) = 2 * x[2] * sin(π * x[1])
g(x) = eu(x)

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
	x = solve_poisson_serendipity(xgrid; rhs = rhs, eu = eu, eu_dx = eu_dx, g = g)

	x0 = solve_poisson_new(xgrid; rhs = rhs, eu = eu, eu_dx = eu_dx, g = g)

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

	return x, x0
end

end #module