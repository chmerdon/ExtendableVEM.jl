# ExtendableVEM.jl
Lowest order virtual Element Method on polygons with ExtendableGrids and ExtendableSparse for the 2D Poisson problem.


The main assembly loop is inspired by:\
Sutton, O.J. The virtual element method in 50 lines of MATLAB. Numer Algor 75, 1141–1159 (2017). https://doi.org/10.1007/s11075-016-0235-3


There is an example that can be run in the julia REPL by:
```julia
julia> using PyPlot
julia> include("examples/Example200_PoissonVEM.jl")
julia> Example200_PoissonVEM.main(; Plotter = PyPlot)
```





