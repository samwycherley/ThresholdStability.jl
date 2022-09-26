export perms, automaton_constructor, isgraphautomaton
using Graphs
"""
    perms(p)

Generates a vector of all sequences of length `p` made up of elements `-1` or `1`.
"""
function perms(p)
    vals = []
    for i in 0:p
        out = unique(permutations([ones(p-i); -ones(i)]))
        for j in 1:length(out)
            push!(vals, out[j])
        end
    end
    return vals
end


"""
    automaton_constructor(Σ::AbstractVector{<:AbstractMatrix})

Generates an automaton tracking admissible transitions between states. `Σ` is a 
vector of matrices corresponding to a TAR model, assuming a default order generated 
by [`perms`](@ref).

    automaton_constructor(Σ::AbstractVector{<:AbstractMatrix}, seqlist)

Generates an automaton tracking admissible transitions between states given a custom list 
of regime sequences.
"""
function automaton_constructor(Σ::AbstractVector{<:AbstractMatrix})
    G = GraphAutomaton(length(Σ))
    p = Int(log2(length(Σ)))
    seqlist = perms(p)
    for i in 1:length(seqlist)
        for j in 1:length(seqlist)
            if seqlist[i][1:end-1] == seqlist[j][2:end]
                add_transition!(G, i, j, i)
            end
        end
    end
    return G
end

function automaton_constructor(Σ::AbstractVector{<:AbstractMatrix}, seqlist)
    G = GraphAutomaton(length(Σ))
    for i in 1:length(seqlist)
        for j in 1:length(seqlist)
            if seqlist[i][1:end-1] == seqlist[j][2:end]
                add_transition!(G, i, j, i)
            end
        end
    end
    return G
end

function HybridSystems.event(A::AbstractAutomaton) end

function HybridSystems.event(A::GraphAutomaton, q, r)
    if has_transition(A, q, r)
        if length(A.Σ[Edge(q, r)]) > 1
            @warn "Edge associated with more than one label! Reporting first label only."
        end
        return collect(values(A.Σ[Edge(q, r)]))[1]
    end
    nothing
end

function HybridSystems.event(A::OneStateAutomaton, q, r)
    @assert q == 1
    @assert r == 1
    collect(1:A.nt)
end

isgraphautomaton(A::AbstractAutomaton) = typeof(A) <: GraphAutomaton
