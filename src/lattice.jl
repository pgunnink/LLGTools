

struct Atom
    # an Atom is a repeatable unit and thus not unique in a lattice
    number::Integer
end

Base.isless(x::Atom, y::Atom) = x.number < y.number

struct LatticeDescription
    primitiveVectors::Vector{Vector{AbstractFloat}}
    basisAtoms::Dict{Atom,Vector{AbstractFloat}}
end

struct Particle
    # a Particle is a unique position in the lattice
    atom::Atom
    UC::Tuple
end


function generateLatticeShape(L::LatticeDescription, shape)    
    d = length(L.primitiveVectors[1]) # dimension 
    
    inside = zeros(d)
    
    res = Dict{Tuple,Vector{Atom}}()
    # you need to add the (0,0) still!
    n = 0
    # ring = collect(Iterators.flatten([x for x in Iterators.product([(-1, 0, 1) for _ in 1:d]...)]))
    
    inside_shape = true
    while inside_shape
        ring = [x for x in Iterators.product([(-n:n) for _ in 1:d]...) if (n in x || -n in x)]
        inside_shape = false
        for i ∈ ring
            pos = inside + sum(i .* L.primitiveVectors)
            uc_atoms = []
            for (at, at_pos) in L.basisAtoms
                if shape(at_pos .+ pos)
                    inside_shape = true
                    push!(uc_atoms, at)
                end
            end
            if !isempty(uc_atoms)
                res[i] = uc_atoms
            end
        end
        n += 1
    end
    res
end



function generateLatticeShape(L::LatticeDescription, shape, region)    
    d = length(L.primitiveVectors[1]) # dimension 
    
    inside = zeros(d)
    
    res = Dict{Tuple,Vector{Atom}}()
    n = 0

    # region has to be: [[x0,x1],[y0,y1]...]

    idx = collect(Iterators.product([x[1]:x[2] for x in region]...))
    for i in idx
        pos = inside + sum(i .* L.primitiveVectors)
        uc_atoms = []
        for (at, at_pos) in L.basisAtoms
            if shape(at_pos .+ pos)
                push!(uc_atoms, at)
            end
        end
        if !isempty(uc_atoms)
            res[i] = uc_atoms
        end
        
        n += 1
    end
    res
end

function absPositionsParticle(p::Particle, l::LatticeDescription)
    pos = sum(p.UC .* l.primitiveVectors)
    at = p.atom
    l.basisAtoms[at] .+ pos
end

function positionsLattice(unit_cells::Dict{Tuple,Vector{Atom}}, L::LatticeDescription)
    # takes a dictionary of unit cells and converts it to explicit positions
    # returns a vector of positions per atom in the system and mapping_index, which is useful if you need to know in which unit cell the atoms are, which is an unique mapping
    positions = Dict{Atom,Vector}()
    mapping = mappingAtoms(unit_cells)
    mapping_index = Dict{Atom,Vector}()
    for (partic, ind) in mapping
        
        pos = sum(partic.UC .* L.primitiveVectors)
        at = partic.atom
        absolute_at_pos = L.basisAtoms[at] .+ pos
        if haskey(positions, at)
            push!(positions[at], absolute_at_pos)
            push!(mapping_index[at], ind)
        else
            positions[at] =  [absolute_at_pos]
            mapping_index[at] =  [ind]
        end
        
    end
    positions, mapping_index
end




function mappingAtoms(unit_cells::Dict{Tuple,Vector{Atom}})
    # downfold the unit_cell into a single interaction matrix, and returns the mapping. We dont care so much about the specific ordering, as long as it is consistent
    mapping = Dict{Particle,Integer}()
    i = 1
    for (ind, uc_atoms) in unit_cells
        for at in uc_atoms
            mapping[Particle(at, ind)] = i
            i += 1
        end
    end
    mapping
end

countAtoms(unit_cells::Dict{Tuple,Vector{Atom}}) = length(mappingAtoms(unit_cells))


function coupleNN(unit_cells::Dict{Tuple,Vector{Atom}}, coupling)
    # coupling should be f(at1, at2, direction)
    d = length(first(keys(unit_cells)))
    nn_list = [x for x in Iterators.product([(-1:1) for _ in 1:d]...) if (1 in x || -1 in x || all(y -> y == 0., x))]

    
    mapping = mappingAtoms(unit_cells)
    N = length(mapping)
    res = zeros(N, N)

    for (ind, uc_atoms) in unit_cells
        for at1 in uc_atoms
            index1 = mapping[Particle(at1, ind)]
            for nn_pos in nn_list
                if (ind .+ nn_pos) in keys(unit_cells)
                    for at2 in unit_cells[ind .+ nn_pos]
                        index2 = mapping[Particle(at2, ind .+ nn_pos)]
                        res[index1, index2] = coupling(at1, at2, nn_pos)
                    end
                end
            end
        end
    end
    res
end

function coupleInternal(unit_cells::Dict{Tuple,Vector{Atom}}, coupling)
    # coupling should be f(at1, at2)    
    mapping = mappingAtoms(unit_cells)
    N = length(mapping)
    res = zeros(N, N)

    for (ind, uc_atoms) in unit_cells
        for at1 in uc_atoms
            index1 = mapping[Particle(at1, ind)]
            for at2 in unit_cells[ind]
                if at1 != at2
                    index2 = mapping[Particle(at2, ind)]
                    res[index1, index2] = coupling(at1, at2)
                end
            end
        end
    end
    res
end