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

include("solvers_new.jl")
export solve_poisson_new ## lowest order VEM with error estimate

include("solvers_serendipity.jl")
export solve_poisson_serendipity, er_f ## 2nd order Serendipity VEM for Poisson problem

include("solvers_stokes_serendipity.jl")
export solve_stokes_serendipity, er_f_2, er_fp ## 2nd order Serendipity VEM (S2-S1) for Stokes problem

include("Gauss_int_2D.jl")
export Gauss_int_2D_tri ## Gauss integration for 2D element

include("poly_base.jl")
export poly_base ## basic 2D function for P_k in each element

include("findm.jl")
export findm ## find the bedge number

include("er.jl")
export er_f, er_f_2, er_fp ## find the bedge number

## very basic grid plotter for polygons
include("plots.jl")
export plot_polygons

end # module ExtendableVEM
