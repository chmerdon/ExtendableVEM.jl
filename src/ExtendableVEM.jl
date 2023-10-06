module ExtendableVEM

using ExtendableGrids
export Polygon2D
using ExtendableSparse
using LinearAlgebra
using CairoMakie
using Makie.GeometryBasics
using ColorSchemes

## some of this should go to ExtendableGrids.jl at some point
include("gridstuff.jl")
export prepare_edges_polygons!, prepare_polygons!
export get_subtriangulation!
export grid_unitsquare
export CellDiameters, CellCenters, SubCells, SubTriangulation

include("solvers.jl")
export solve_poisson ## lowest order VEM following Sutton's 50 lines of Matlab

## very basic grid plotter for polygons
include("plots.jl")
export plot_polygons

end # module ExtendableVEM
