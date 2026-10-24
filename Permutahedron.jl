using PositiveChartsToricVarieties, Oscar

#Sets variables for defining polynomials
d = 3
N, t = polynomial_ring(QQ,:t=>1:d)



function set_to_polynomial(s)
    poly = 1
        for i in 2:length(s)
            ind = s[1]+1
            temp = 1
            for j in ind:s[i]
                temp *= t[j]
            end
            poly += temp
        end
    return poly
end


f = Vector{QQMPolyRingElem}()
for k in 2:d+1
    subs = subsets([i-1 for i in 1:d+1],k)
    for s in subs
        push!(f, set_to_polynomial(s))
    end
end

J = Y_variety(f)
tropical_variety(J, skip_saturation=true, skip_decomposition=true)






#2 dim example
test, I, J =test_sat_conjecture(f)
S = gens(J)[1].parent
y = gens(S)

u = [y[1]+ y[5]+ y[2]*y[3]-1,
y[2]+ (y[1]+ y[5])*y[4]-1,
y[4] + y[2]*y[6] -1,
y[6] + y[3]*y[4] -1,
y[3] + (y[1]+ y[5])*y[6] -1]
U = ideal(u)
J == U