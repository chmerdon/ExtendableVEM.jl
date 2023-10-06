# ExtendableVEM.jl
Lowest order virtual Element Method on polygons for the 2D Poisson problem
based on the grid manager [ExtendableGrids](https://github.com/j-fu/ExtendableGrids.jl) and the sparse matrix manager [ExtendableSparse](https://github.com/j-fu/ExtendableSparse.jl).


The main assembly loop is inspired by:\
Sutton, O.J. The virtual element method in 50 lines of MATLAB. Numer Algor 75, 1141–1159 (2017). https://doi.org/10.1007/s11075-016-0235-3


There is an example that can be run in the julia REPL by:
```julia
julia> using PyPlot
julia> include("examples/Example200_PoissonVEM.jl")
julia> Example200_PoissonVEM.main(; Plotter = PyPlot)
```





