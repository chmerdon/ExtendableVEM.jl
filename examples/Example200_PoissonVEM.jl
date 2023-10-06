module Example200_PoissonVEM

using ExtendableVEM
using ExtendableGrids
using GridVisualize

# forcing function
rhs(x) = sin(4 * π * x[1]) * sin(2 * π * x[2])
# Boundary condition
g(x) = 0

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
	x = solve_poisson(xgrid; rhs = rhs, g = g)

	## calculate average values at polygon centers for plot on subtriangulation
	subtriangulation = xgrid[SubTriangulation]
	ncells = num_cells(xgrid)
	append!(x, zeros(Float64, ncells))
	cellnodes = xgrid[CellNodes]
	for cell ∈ 1:num_cells(xgrid)
		x[nnodes+cell] = sum(view(x, view(cellnodes, :, cell))) / num_targets(cellnodes, cell)
	end

	## plot polygons
	f = plot_polygons(xgrid)

	## plot subtriangulations and solution
	pl = GridVisualizer(; Plotter = Plotter, layout = (1, 2), clear = true, size = (800, 400))
	gridplot!(pl[1, 1], subtriangulation, linewidth = 1)
	scalarplot!(pl[1, 2], subtriangulation, x; Plotter = Plotter)

	return f
end

end #module