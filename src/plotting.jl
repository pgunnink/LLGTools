

function positionsPlottingXY(unit_cells, l::LatticeDescription)
    # gets absolute positions for plotting, with the ordering the same as returned by mapping
    positions, mapping_index = positionsLattice(unit_cells, l)

    res = [[0.0, 0.0] for _ in 1:length(mappingAtoms(unit_cells))]
    for x in keys(l.basisAtoms)
        for i = 1:length(positions[x])
            res[mapping_index[x][i], :] = positions[x][i]
        end
    end
    return res
end
