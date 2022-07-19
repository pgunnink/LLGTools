function get_conversion_Nd(l::LatticeDescription)
    if length(l.primitiveVectors) == 1
        # we have only x points, so add artificial y direction
        return pos -> length(pos) == 1 ? Point2f0((pos[1], 0.0)...) : Point2f0(pos...)
    elseif length(l.primitiveVectors) == 2
        return pos -> Point2f0(pos...)
    elseif length(l.primitiveVectors) == 3
        return pos -> Point3f0(pos...)
    end
end


function positionsPlottingPerAtom(unit_cells, l::LatticeDescription)
    # gets absolute positions for plotting, including an artificial y direction if needed, per atom for easy plotting
    positions, mapping_index = positionsLattice(unit_cells, l)
    points = []
    convert_Nd = get_conversion_Nd(l)
    for x in keys(sort(l.basisAtoms))
        p = [convert_Nd(positions[x][i]) for i = 1:length(positions[x])]
        push!(points, p)
    end
    return points
end


function positionsPlotting(unit_cells, l::LatticeDescription)
    # gets absolute positions for plotting, including an artificial y direction if needed, with the ordering the same as returned by mapping
    positions, mapping_index = positionsLattice(unit_cells, l)
    convert_Nd = get_conversion_Nd(l)

    res = zeros(Point{2,Float32}, length(mappingAtoms(unit_cells)))
    for x in keys(l.basisAtoms)
        for i = 1:length(positions[x])
            res[mapping_index[x][i]] = convert_Nd(positions[x][i])
        end
    end
    return res
end