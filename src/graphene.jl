function nn_coupling_honeycomb(a, b, p)
    if a == Atom(1) && b == Atom(2)
        if p == (1, -1) || p == (0, -1) || p == (0, 0)
            return 1
        end
    end
    if b == Atom(1) && a == Atom(2)
        if p == (-1, 1) || p == (0, 1) || p == (0, 0)
            return 1
        end
    end
    return 0
end

function nnn_coupling_honeycomb(a, b, p)
    if a == b
        if p != (0, 0)
            return 1
        end
    end
    return 0
end

function dmi_interaction_honeycomb(a, b, p)
    if p == (1, 1) || p == (-1, -1)
        return 0
    end
    if a == b
        if p != (0, 0)
            # if a == Atom(1)
            #     return 1
            # else
            #     return -1
            # end
            if p == (1, 0) || p == (0, -1) || p == (-1, 1)
                return 1
            elseif p == (-1, 0) || p == (0, 1) || p == (1, -1)
                return -1
            end
        end
    end
    return 0
end
