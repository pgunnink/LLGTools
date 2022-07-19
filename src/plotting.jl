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

function plot_atoms!(unit_cells, l::LatticeDescription; markersize = 15)
    points = positionsPlottingPerAtom(unit_cells, l)
    colors = [:red :blue :green :yellow]

    for (p, c) in zip(points, colors)
        scatter!(p, markersize = markersize, color = c)
    end

end

function plot_density!(unit_cells, l::LatticeDescription, dens)
    points = positionsPlotting(unit_cells, l)
    cg = cgrad(:viridis)
    dens = dens ./ maximum(dens)
    scatter!(points, colormap = cg, color = dens, markersize = 100 / sqrt(length(points)))
end

function plot_density(unit_cells, l::LatticeDescription, dens)

    f = Figure()
    ax = Axis(f[1, 1])

    plot_density!(unit_cells, l, dens)
    current_axis().aspect = AxisAspect(1)

    return f
end
function plot_atoms(unit_cells, l::LatticeDescription; markersize = 15, force_2D = false)
    f = Figure()
    if length(l.primitiveVectors) < 3 || force_2D
        ax = Axis(f[1, 1])
        current_axis().aspect = AxisAspect(1)
    else
        lscene = LScene(f[1, 1], scenekw = (camera = cam3d!, raw = false))
    end
    plot_atoms!(unit_cells, l, markersize = markersize)


    return f
end

function plot_connections!(unit_cells, l::LatticeDescription, coupling::Matrix; linecolor = :grey, linewidth = 1, plot_arrows = false)
    mapping = mappingAtoms(unit_cells)
    inversemapping = Dict(value => key for (key, value) in mapping)

    x = []
    y = []
    z = []
    u = []
    v = []
    w = []
    directions = []
    for i = 1:size(coupling)[1]
        for j = 1:size(coupling)[1]
            if coupling[i, j] != 0
                p1, p2 = (
                    absPositionsParticle(inversemapping[i], l),
                    absPositionsParticle(inversemapping[j], l))
                if length(p1) == 1
                    push!(p1, 0, 0)
                    push!(p2, 0, 0)
                elseif length(p1) == 2
                    push!(p1, 0)
                    push!(p2, 0)
                end
                # push!(arrow_positions, p1)
                # push!(directions, p2 - p1)
                push!(x, p1[1])
                push!(y, p1[2])
                push!(z, p1[3])

                push!(u, p2[1] - p1[1])
                push!(v, p2[2] - p1[2])
                push!(w, p2[3] - p1[3])
                if coupling[i, j] > 0
                    push!(directions, '▼')
                else
                    push!(directions, '▲')
                end
                # push!(x, Point2f0(p1...))
                # push!(u, Point2f0((p2 - p1)...))
            end
        end
    end
    if length(l.primitiveVectors) < 3
        arrows!(x, y, u, v; linecolor = linecolor, linewidth = linewidth, arrowsize = 0.0)
        if plot_arrows
            for i in 1:length(x)
                arrows!([x[i]], [y[i]], 0.6 .* [u[i]], 0.6 .* [v[i]]; linecolor = linecolor, linewidth = linewidth, arrowhead = directions[i])
            end
        end
    else
        arrows!(x, y, z, u, v, w; linecolor = linecolor, linewidth = linewidth)
    end

    # current_axis().aspect = AxisAspect(1)

end

function plot_connections(unit_cells, l::LatticeDescription, coupling::Matrix; linecolor = :grey, linewidth = 1, markersize = 15, plot_arrows = false)
    fig = plot_atoms(unit_cells, l, markersize = markersize)
    plot_connections!(unit_cells, l, coupling; linecolor = linecolor, linewidth = linewidth, plot_arrows = plot_arrows)
    return fig
end