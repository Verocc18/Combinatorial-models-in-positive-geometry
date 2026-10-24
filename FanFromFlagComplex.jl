using Oscar,Graphs


#Builds the candidate rays for each vertex of the simplicial complex according to the exponents a 
#g is the graph complement of the 1-skeleton of the simplicial complex
function neighbor_vectors_from_graph(g;a=nothing)
    n = Graphs.nv(g)
    vectors = Dict{Int, Vector{Int}}()
    for v in Graphs.vertices(g)
        x = zeros(Int, n)
        x[v] = 1
        for u in Graphs.neighbors(g, v)
            if isnothing(a)
                x[u] = -1
            else
                x[u] = -a[v,u]
            end
        end
        vectors[v] = x
    end
    return vectors
end

#Builds the candidate fans from the graph and the exponent vectors
#Setting stop=true makes the function stop as soon as it finds a fan
function fans_from_vertices(g;a=nothing,stop=false)
    gc = Graphs.complement(g)
    vtcs = Graphs.maximal_cliques(gc)
    K = Oscar.simplicial_complex(vtcs)
    target_f = Oscar.f_vector(K)
    IM = Oscar.incidence_matrix([v for v in vtcs])
    vectors = neighbor_vectors_from_graph(g;a)
    goodFans = []
    for c in vtcs
        try 
            R = hcat([vectors[i] for i in c]...)
            cones = [positive_hull(R[collect(row(IM, i)), :]) for i in 1:number_of_rows(IM)]
            Σ = polyhedral_fan(cones)

            if Oscar.f_vector(Σ) == target_f
                if stop
                    return Σ
                end
                push!(goodFans, Σ)
            end
        catch
            @warn "Could not build a fan."
        end
    end
    return goodFans
end

#Tests Shephard's criterion for projectivity of a simplicial fan
function shephard_criterion(Σ)
    @assert is_simplicial(Σ) "The fan is not simplicial."
    @assert is_complete(Σ) "The fan is not complete." 

    F = matrix(QQ, rays(Σ))
    K = kernel(F)
    M = maximal_cones(IncidenceMatrix, Σ)
    gale_cones = []

    for i in range(1,nrows(M))
        J = setdiff(1:nrows(F), row(M,i))
        C = positive_hull(transpose(K[:, J]))
        push!(gale_cones,C)
    end
    
    int_cone = intersect(gale_cones...)
    return is_fulldimensional(int_cone)
end
    








#E6-Pezzotope
#Graph for E6
g_E6 = SimpleGraph(15)

# start with Petersen graph 

# Inner pentagon
Graphs.add_edge!(g_E6, 1, 2)
Graphs.add_edge!(g_E6, 1, 5)
Graphs.add_edge!(g_E6, 3, 5)
Graphs.add_edge!(g_E6, 3, 4)
Graphs.add_edge!(g_E6, 2, 4)

# Outer pentagon 
Graphs.add_edge!(g_E6, 7, 6)
Graphs.add_edge!(g_E6, 6, 8)
Graphs.add_edge!(g_E6, 8, 10)
Graphs.add_edge!(g_E6, 10, 9)
Graphs.add_edge!(g_E6, 9, 7)

# Spokes
Graphs.add_edge!(g_E6, 1, 7)
Graphs.add_edge!(g_E6, 3, 6)
Graphs.add_edge!(g_E6, 2, 8)
Graphs.add_edge!(g_E6, 5, 10)
Graphs.add_edge!(g_E6, 4, 9)

# --- interweaving 10-cycle ---
Graphs.add_edge!(g_E6, 1, 15)
Graphs.add_edge!(g_E6, 15, 3)
Graphs.add_edge!(g_E6, 3, 14)
Graphs.add_edge!(g_E6, 14, 2)
Graphs.add_edge!(g_E6, 2, 12)
Graphs.add_edge!(g_E6, 12, 5)
Graphs.add_edge!(g_E6, 5, 11)
Graphs.add_edge!(g_E6, 11, 4)
Graphs.add_edge!(g_E6, 4, 13)
Graphs.add_edge!(g_E6, 13, 1)

# Outer → middle connections

Graphs.add_edge!(g_E6, 7, 14)
Graphs.add_edge!(g_E6, 7, 11)
Graphs.add_edge!(g_E6, 6, 12)
Graphs.add_edge!(g_E6, 6, 13)
Graphs.add_edge!(g_E6, 8, 11)
Graphs.add_edge!(g_E6, 8, 15)
Graphs.add_edge!(g_E6, 10, 13)
Graphs.add_edge!(g_E6, 10, 14)
Graphs.add_edge!(g_E6, 9, 15)
Graphs.add_edge!(g_E6, 9, 12)

# K5 on middle pentagon

Graphs.add_edge!(g_E6, 15, 14)
Graphs.add_edge!(g_E6, 14, 12)
Graphs.add_edge!(g_E6, 12, 11)
Graphs.add_edge!(g_E6, 11, 13)
Graphs.add_edge!(g_E6, 13, 15)
Graphs.add_edge!(g_E6, 15, 12)
Graphs.add_edge!(g_E6, 15, 11)
Graphs.add_edge!(g_E6, 14, 11)
Graphs.add_edge!(g_E6, 14, 13)
Graphs.add_edge!(g_E6, 12, 13)

goodFans = fans_from_vertices(g_E6)
for Σ in goodFans 
    @assert is_regular(Σ) "The fan is not projective."
end




#E7-Pezzotope
adjacency = Dict{Int, Vector{Int}}()

# u1 + u2*u21*u22*u23*u24*u28*u3*u30*u32
adjacency[1]  = [2,21,22,23,24,28,3,30,32]

# u2 + u1*u19*u20*u25*u26*u27*u29*u31*u33*u8
adjacency[2]  = [1,19,20,25,26,27,29,31,33,8]

# u3 + u1*u11*u13*u14*u16*u18*u25*u33*u4*u6
adjacency[3]  = [1,11,13,14,16,18,25,33,4,6]

# u4 + u12*u17*u19*u21*u29*u3*u5*u7*u8
adjacency[4]  = [12,17,19,21,29,3,5,7,8]

# u5 + u13*u20*u26*u28*u30*u31*u32*u4*u6*u9
adjacency[5]  = [13,20,26,28,30,31,32,4,6,9]

# u6 + u10*u12*u14*u17*u19*u21*u22*u24*u25*u26*u27*u28*u29*u3*u34*u5*u8
adjacency[6]  = [10,12,14,17,19,21,22,24,25,26,27,28,29,3,34,5,8]

# u7 + u10*u14*u22*u24*u25*u26*u27*u28*u34*u4
adjacency[7]  = [10,14,22,24,25,26,27,28,34,4]

# u8 + u11*u13*u14*u16*u18*u2*u21*u22*u23*u24*u25*u28*u30*u32*u33*u4*u6
adjacency[8]  = [11,13,14,16,18,2,21,22,23,24,25,28,30,32,33,4,6]

# u9 + u11*u14*u16*u19*u21*u22*u23*u24*u25*u27*u29*u33*u5
adjacency[9]  = [11,14,16,19,21,22,23,24,25,27,29,33,5]

# u10 + u15*u16*u23*u29*u30*u31*u33*u6*u7
adjacency[10] = [15,16,23,29,30,31,33,6,7]

# u11 + u15*u17*u22*u26*u27*u28*u29*u3*u30*u31*u34*u8*u9
adjacency[11] = [15,17,22,26,27,28,29,3,30,31,34,8,9]

# u12 + u14*u15*u16*u22*u23*u24*u25*u26*u27*u28*u29*u30*u31*u33*u34*u4*u6
adjacency[12] = [14,15,16,22,23,24,25,26,27,28,29,30,31,33,34,4,6]

# u13 + u14*u15*u16*u17*u19*u21*u22*u23*u24*u25*u26*u27*u28*u29*u3*u30*u31*u33*u34*u5*u8
adjacency[13] = [14,15,16,17,19,21,22,23,24,25,26,27,28,29,3,30,31,33,34,5,8]

# u14 + u12*u13*u17*u19*u20*u21*u26*u28*u29*u3*u30*u31*u32*u6*u7*u8*u9
adjacency[14] = [12,13,17,19,20,21,26,28,29,3,30,31,32,6,7,8,9]

# u15 + u10*u11*u12*u13*u19*u20*u21*u24*u25*u32
adjacency[15] = [10,11,12,13,19,20,21,24,25,32]

# u16 + u10*u12*u13*u17*u19*u20*u21*u22*u24*u25*u26*u27*u28*u29*u3*u30*u31*u32*u34*u8*u9
adjacency[16] = [10,12,13,17,19,20,21,22,24,25,26,27,28,29,3,30,31,32,34,8,9]

# u17 + u11*u13*u14*u16*u19*u20*u21*u22*u23*u24*u25*u26*u27*u28*u29*u30*u31*u32*u33*u4*u6
adjacency[17] = [11,13,14,16,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,4,6]

# u18 + u19*u20*u21*u22*u23*u24*u25*u26*u27*u28*u29*u3*u30*u31*u32*u33*u8
adjacency[18] = [19,20,21,22,23,24,25,26,27,28,29,3,30,31,32,33,8]

# u19 + u13*u14*u15*u16*u17*u18*u2*u22*u26*u28*u30*u31*u32*u34*u4*u6*u9
adjacency[19] = [13,14,15,16,17,18,2,22,26,28,30,31,32,34,4,6,9]

# u20 + u14*u15*u16*u17*u18*u2*u22*u34*u5
adjacency[20] = [14,15,16,17,18,2,22,34,5]

# u21 + u1*u13*u14*u15*u16*u17*u18*u22*u25*u26*u27*u28*u29*u30*u31*u33*u34*u4*u6*u8*u9
adjacency[21] = [1,13,14,15,16,17,18,22,25,26,27,28,29,30,31,33,34,4,6,8,9]

# u22 + u1*u11*u12*u13*u16*u17*u18*u19*u20*u21*u25*u26*u29*u30*u31*u32*u33*u6*u7*u8*u9
adjacency[22] = [1,11,12,13,16,17,18,19,20,21,25,26,29,30,31,32,33,6,7,8,9]

# u23 + u1*u10*u12*u13*u17*u18*u26*u34*u8*u9
adjacency[23] = [1,10,12,13,17,18,26,34,8,9]

# u24 + u1*u12*u13*u15*u16*u17*u18*u26*u29*u30*u31*u33*u34*u6*u7*u8*u9
adjacency[24] = [1,12,13,15,16,17,18,26,29,30,31,33,34,6,7,8,9]

# u25 + u12*u13*u15*u16*u17*u18*u2*u21*u22*u26*u28*u29*u3*u30*u31*u32*u34*u6*u7*u8*u9
adjacency[25] = [12,13,15,16,17,18,2,21,22,26,28,29,3,30,31,32,34,6,7,8,9]

# u26 + u11*u12*u13*u14*u16*u17*u18*u19*u2*u21*u22*u23*u24*u25*u29*u30*u32*u33*u5*u6*u7
adjacency[26] = [11,12,13,14,16,17,18,19,2,21,22,23,24,25,29,30,32,33,5,6,7]

# u27 + u11*u12*u13*u16*u17*u18*u2*u21*u30*u32*u6*u7*u9
adjacency[27] = [11,12,13,16,17,18,2,21,30,32,6,7,9]

# u28 + u1*u11*u12*u13*u14*u16*u17*u18*u19*u21*u25*u29*u33*u5*u6*u7*u8
adjacency[28] = [1,11,12,13,14,16,17,18,19,21,25,29,33,5,6,7,8]

# u29 + u10*u11*u12*u13*u14*u16*u17*u18*u2*u21*u22*u24*u25*u26*u28*u30*u32*u34*u4*u6*u9
adjacency[29] = [10,11,12,13,14,16,17,18,2,21,22,24,25,26,28,30,32,34,4,6,9]

# u30 + u1*u10*u11*u12*u13*u14*u16*u17*u18*u19*u21*u22*u24*u25*u26*u27*u29*u33*u34*u5*u8
adjacency[30] = [1,10,11,12,13,14,16,17,18,19,21,22,24,25,26,27,29,33,34,5,8]

# u31 + u10*u11*u12*u13*u14*u16*u17*u18*u19*u2*u21*u22*u24*u25*u32*u34*u5
adjacency[31] = [10,11,12,13,14,16,17,18,19,2,21,22,24,25,32,34,5]

# u32 + u1*u14*u15*u16*u17*u18*u19*u22*u25*u26*u27*u29*u31*u33*u34*u5*u8
adjacency[32] = [1,14,15,16,17,18,19,22,25,26,27,29,31,33,34,5,8]

# u33 + u10*u12*u13*u17*u18*u2*u21*u22*u24*u26*u28*u3*u30*u32*u34*u8*u9
adjacency[33] = [10,12,13,17,18,2,21,22,24,26,28,3,30,32,34,8,9]

# u34 + u11*u12*u13*u16*u19*u20*u21*u23*u24*u25*u29*u30*u31*u32*u33*u6*u7
adjacency[34] = [11,12,13,16,19,20,21,23,24,25,29,30,31,32,33,6,7]


# Build undirected graph
g_E7 = SimpleGraph(34)
for (i, neighbors) in adjacency
    for j in neighbors
        Graphs.add_edge!(g_E7, i, j)
    end
end


Σ=fans_from_vertices(g_E7;stop=true)

