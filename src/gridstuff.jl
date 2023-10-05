
abstract type CellDiameters <: AbstractGridFloatArray1D end
abstract type CellCenters <: AbstractGridFloatArray2D end
abstract type SubCells <: AbstractGridAdjacency end
abstract type SubTriangulation <: AbstractGridComponent end

function prepare_subtriangulation!(xgrid::ExtendableGrid{Tc,Ti}) where {Tc,Ti}
    xCellNodes = xgrid[CellNodes]
    xCoordinates = xgrid[Coordinates]
    ncells = num_sources(xCellNodes)
    nnodes = num_nodes(xgrid)
    ntriangles = 0
    for cell = 1 : ncells
        ntriangles += num_targets(xCellNodes, cell)
    end  

    cellregions = xgrid[CellRegions]
    cellnodes_triangles = zeros(Ti, 3, ntriangles)
    cellregions_triangles = zeros(Ti, ntriangles)
    subcells = deepcopy(xCellNodes)
    ntriangles = 0
    for cell = 1 : ncells
        nvertices = num_targets(xCellNodes, cell)
        for j = 1 : nvertices
            ntriangles += 1
            cellnodes_triangles[1, ntriangles] = xCellNodes[j, cell]
            cellnodes_triangles[2, ntriangles] = xCellNodes[mod1(j + 1, nvertices), cell]
            cellnodes_triangles[3, ntriangles] = nnodes + cell
            subcells.colentries[subcells.colstart[cell]+j-1] = ntriangles
            cellregions_triangles[ntriangles] = cellregions[cell]
        end 
    end
    subtriangulation = ExtendableGrid{Tc,Ti}()
    subtriangulation[Coordinates] = [xCoordinates xgrid[CellCenters]]
    subtriangulation[CellNodes] = cellnodes_triangles
    subtriangulation[CellRegions] = cellregions_triangles
    subtriangulation[BFaceNodes] = xgrid[BFaceNodes]
    subtriangulation[BFaceRegions] = xgrid[BFaceRegions]
    subtriangulation[CellGeometries] = VectorOfConstants{ElementGeometries,Int}(Triangle2D, ntriangles)
    subtriangulation[CoordinateSystem] = Cartesian2D


    triangle_volumes = subtriangulation[CellVolumes]
    cellvolumes = zeros(Tc, ncells)
    ntriangles = 0
    for cell = 1 : ncells
        nvertices = num_targets(xCellNodes, cell)
        for k = 1 : nvertices
            ntriangles += 1
            cellvolumes[cell] += triangle_volumes[ntriangles]
        end
    end

    xgrid[SubCells] = subcells
    xgrid[CellVolumes] = cellvolumes
    xgrid[SubTriangulation] = subtriangulation
end

function prepare_polygons!(xgrid::ExtendableGrid{Tc,Ti}) where {Tc,Ti}

    ## first enumerate edges
	prepare_edges_polygons!(xgrid)

    xCellNodes = xgrid[CellNodes]
    xCoordinates = xgrid[Coordinates]
    ncells = num_sources(xCellNodes)
    dim = size(xCoordinates,1)

    xCellDiameters = zeros(Tc, ncells)
    xCellCenters = zeros(Tc, dim, ncells)

    max = 0.0
    length = 0.0
    for cell = 1 : ncells
        nvertices = num_targets(xCellNodes, cell)
        max = 0.0
        for j = 1 : nvertices
            node_j = xCellNodes[j, cell]
            for k = j+1 : nvertices
                length = sqrt((xCoordinates[1, node_j] - xCoordinates[1, xCellNodes[k, cell]]).^2 + (xCoordinates[2, node_j] - xCoordinates[2, xCellNodes[k, cell]]).^2)
                if length > max
                    max = length
                end
            end
            for d = 1 : dim
                xCellCenters[d,cell] += xCoordinates[d, node_j]
            end
        end
        xCellCenters[:, cell] ./= nvertices
        xCellDiameters[cell] = max
    end  
    xgrid[CellDiameters] = xCellDiameters
    xgrid[CellCenters] = xCellCenters

    prepare_subtriangulation!(xgrid)
end

function prepare_edges_polygons!(xgrid::ExtendableGrid{Tc,Ti}) where {Tc,Ti}

    xCellNodes = xgrid[CellNodes]
    ncells = num_sources(xCellNodes)
    nnodes = num_sources(xgrid[Coordinates])

    # transpose CellNodes to get NodeCells
    xNodeCells = atranspose(xCellNodes)

    # instantiate new empty adjacency fields
    xFaceCells = zeros(Ti,0) # cells are appended and at the end rewritten into 2,nfaces array
    xFaceNodes = zeros(Ti,0)
    xCellFaces::Union{VariableTargetAdjacency{Ti}, Matrix{Ti}} = VariableTargetAdjacency(Ti)
    xCellFaceSigns::Union{VariableTargetAdjacency{Ti}, Matrix{Ti}} = VariableTargetAdjacency(Ti)
    # pre-allocate xCellFaces, each polygon has nnodes many edges
    for cell = 1 : ncells
        nedges = num_targets(xCellNodes, cell)
        append!(xCellFaces,zeros(Ti,nedges))
        append!(xCellFaceSigns,zeros(Ti,nedges))
    end  

    # temporary variables
    face::Ti = 0
    cell::Ti = 0
    cell2::Ti = 0
    nneighbours::Ti = 0
    faces_per_cell::Ti = 0
    faces_per_cell2::Ti = 0
    current_item::Array{Ti,1} = zeros(Ti,2)
    current_item2::Array{Ti,1} = zeros(Ti,2)
    flag4item::Array{Bool,1} = zeros(Bool, nnodes)
    no_neighbours_found::Bool = true
    same_face::Bool = false
    
    # loop over cells
    for cell = 1 : ncells
        faces_per_cell = num_targets(xCellNodes, cell)

        # loop over cell faces
        for k = 1 : faces_per_cell

            # check if face is already known to cell
            if xCellFaces[k, cell] > 0
                continue;
            end    

            current_item[1] = xCellNodes[k, cell] 
            current_item[2] = k == faces_per_cell ? xCellNodes[1, cell] : xCellNodes[k+1, cell]
            
            # flag face nodes and commons4cells
            for j = 1:2
                flag4item[current_item[j]] = true
            end

            # loop over neighbours of first node
            nneighbours = num_targets(xNodeCells, current_item[1])
            no_neighbours_found = true
            for n = 1 : nneighbours
                cell2 = xNodeCells[n, current_item[1]]

                # skip if cell2 is the same as cell
                if (cell == cell2) 
                    continue; 
                end

                # find face enumeration rule
                faces_per_cell2 = num_targets(xCellNodes, cell2)

                # loop over faces face2 of adjacent cell2
                for f2 = 1 : faces_per_cell2

                    # check if face f2 is already known to cell2
                    if xCellFaces[f2,cell2] != 0
                        continue;
                    end    

                    # otherwise compare nodes of face and face2
                    current_item2[1] = xCellNodes[f2, cell2] 
                    current_item2[2] = f2 == faces_per_cell2 ? xCellNodes[1, cell2] : xCellNodes[f2+1, cell2]
                    same_face = (current_item[1] == current_item2[2]) && (current_item[2] == current_item2[1])
                    
                    # if all nodes are the same, register face
                    if (same_face)
                        no_neighbours_found = false
                        face += 1
                        push!(xFaceCells,cell)
                        push!(xFaceCells,cell2)
                        xCellFaces.colentries[xCellFaces.colstart[cell]+k-1] = face
                        xCellFaces.colentries[xCellFaces.colstart[cell2]+f2-1] = face
                        xCellFaceSigns.colentries[xCellFaceSigns.colstart[cell]+k-1] = 1
                        xCellFaceSigns.colentries[xCellFaceSigns.colstart[cell2]+f2-1] = -1
                        append!(xFaceNodes, current_item)
                        break;
                    end
                end
            end

            # if no common neighbour cell is found, register face (boundary faces)
            if no_neighbours_found == true
                face += 1
                push!(xFaceCells,cell)
                push!(xFaceCells,0)
                xCellFaces.colentries[xCellFaces.colstart[cell]+k-1] = face
                xCellFaceSigns.colentries[xCellFaceSigns.colstart[cell]+k-1] = 1
                append!(xFaceNodes, current_item)
            end

            # reset flag4item
            for j = 1:2
                flag4item[current_item[j]] = false 
            end
        end    
    end

    xFaceNodes = reshape(xFaceNodes,(2,Ti(length(xFaceNodes)/2)))
    xgrid[FaceGeometries] = VectorOfConstants{ElementGeometries,Int}(facetype_of_cellface(Edge1D, 1), face)
    xgrid[CellFaces] = xCellFaces
    xgrid[CellFaceSigns] = xCellFaceSigns
    xgrid[FaceCells] = reshape(xFaceCells,(2,face))
    xgrid[FaceNodes] = xFaceNodes
end