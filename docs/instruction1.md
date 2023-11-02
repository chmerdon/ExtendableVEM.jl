new: 
    src/
        er.jl          # calculate error function
        findm.jl       # find the boundary edge number 
        Gauss_int_2D.jl  # Gauss integration
        poly_base.jl     # base function
        solvers_new.jl      # 1-VEM for Poisson with H1 and L2 error     right
        solvers_seredipity.jl   # 2-SVEM for Poisson with H1 and L2 error    not all right but small
        solvers_stokes_seredipity.jl  # 2-SVEM for Stokes with H1 and L2 error    not right
    examples/
        Example_PoissonSVEM.jl   # show the Poisson problem by solvers_new.jl and solvers_seredipity.jl
        Example_StokesSVEM.jl    # show the Stokes problem by solvers_stokes_seredipity.jl
update:
    src/
        ExtendableVEM   # to include the new .jl 

addition:
    using the solvers_new.jl can get the right error order so the er.jl, Gauss_int_2D.jl, poly_base.jl functions should be right
    using the solvers_seredipity.jl can get the small error but not the right error order 
    using the solvers_stokes_seredipity.jl get the bad error 

    guess there are some little parts wrong in the code "solvers_seredipity.jl" may due to the misunderstanding of Julia programming
        have implemented the method in one element to check the projection matrix and stiffness matrix, they seemed all right.

    also tried re-deduce the construction of the projection matrix under the vector, and no problems were found.
