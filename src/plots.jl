function plot_polygons(xgrid)
    f = Figure()
    cellnodes = xgrid[CellNodes]
    coordinates = xgrid[Coordinates]
	Axis(f[1, 1])
	for cell ∈ 1:num_cells(xgrid)
		poly!(Makie.Polygon([Point2f(view(coordinates, :, n)) for n in view(cellnodes, :, cell)]), strokewidth = 1, strokecolor = :black, color = :gray)
	end
    return f
end