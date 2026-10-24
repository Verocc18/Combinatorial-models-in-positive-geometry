using PositiveChartsToricVarieties, Oscar

#Sets variables for defining polynomials
d = 2
N, t = polynomial_ring(QQ,:t=>1:d)

# Polynomials for M_{0,n} where d=n-3
f = Vector{QQMPolyRingElem}()
for i in 1:d
    for j in i:d
        # Start with the constant term 1
        term = 1
        # Build sum of products t_k * ... * t_j for k = j down to i
        for k in i:j
            produ = prod(t[k:j])  # product t_k*t_{k+1}*...*t_j
            term += produ
        end
        push!(f, term)
    end
end

#Polynomials for n-gons
n = 7
f = [1 + t[1]; 1 + t[2]]
for i in 1:(n-4)
    push!(f, (1 + t[1])^i + t[2])
end


#Parameterization of positive chart
varphi = parameterization_U(f)
#Defining equations
J = U_variety(f)
#Unimodular matrix
M, F = unimod_matrix(f)



#Hexagons examples
#Example 1
f = [1 + t[1]; 1 + t[2]; 1 + t[1] + t[1]*t[2]; 1 + t[2]+t[1]*t[2]]

#Creates positive chart and tests saturation conjecture
test, I, J =test_sat_conjecture(f)

#Checks binary presentation
S = gens(J)[1].parent
u = gens(S)
u_eqs=[u[1] + u[3]*u[4]*u[5] - 1,
  u[2] + u[3]*u[4]*u[6] - 1,
  u[1]*u[2]*u[5] + u[3] - 1,
  u[1]*u[2]*u[6] + u[4] - 1,
  u[1]*u[3]*u[6]^2 + u[5] - 1,
  u[2]*u[4]*u[5]^2 + u[6] - 1]
U = ideal(u_eqs)
J == U

#Example 2
f=[1 + t[1]; 1 + t[2]; 1 + t[1] + t[1]*t[2]; 1 + t[1]*t[2]]

#Creates positive chart and tests saturation conjecture
test, I, J =test_sat_conjecture(f)

#Checks binary presentation
S = gens(I)[1].parent
u = gens(S)
u_eqs = [u[6] + u[3]*u[4]*u[5] - 1,
  u[6]*u[1]*u[4] + u[3] - 1,
  u[1] + u[3]^2*u[5]*u[2] - 1,
  u[6]^2*u[1]*u[2] + u[5] - 1,
  u[6]*u[3]*u[2] + u[4] - 1,
  u[1]*u[5]*u[4]^2 + u[2] - 1]
U = ideal(u)
J == U


#P1xP1 examples

#Choice of matrix A in SL2(Z)
a=2
b=1
c=1
d=1
a*d-b*c == 1

#Polynomials with divisors corresponding to A
f = [(1+t[1])^a*(1+t[2])^b; (1+t[1])^c*(1+t[2])^d]

#Computes variety and tests saturation conjecture
test, I, J =test_sat_conjecture(f)




#Permutahedron example
P = project_full(permutahedron(3))
f, M = unimod_nef_polynomials(P)
Y_variety(f)

#Polynomials to fit order in Example 4.3
ind = [10,3,2,5,4,1,7,6,8,9,11]
newf = [f[ind[i]] for i in 1:length(ind)]
unimod_matrix(newf)



#Different A3 example
f = [1+t[1]; 1+t[2]; 1+t[3]; 1+t[3]+t[2]*t[3]; 1+t[2]+t[1]*t[2]; 1+t[1]+t[1]*t[3]]


#D4 examples
f = [t[4]+1;
t[3]+1;
t[2]+1;
t[1]+1;
t[2]*t[4]+t[4]+1;
t[2]*t[3]+t[3]+1;
t[1]*t[2]+t[1]+1;
t[2]*t[3]*t[4]+t[3]*t[4]+t[3]+t[4]+1;
t[1]*t[2]*t[4]+t[1]*t[4]+t[1]+t[4]+1;
t[1]*t[2]*t[3]+t[1]*t[3]+t[1]+t[3]+1;
t[1]*t[2]*t[3]*t[4]+t[1]*t[3]*t[4]+t[1]*t[3]+t[1]*t[4]+t[1]+t[3]*t[4]+t[3]+t[4]+1;
t[1]*t[2]^2*t[3]*t[4]+2*t[1]*t[2]*t[3]*t[4]+t[1]*t[2]*t[3]+t[1]*t[2]*t(4)+t[1]*t(3)*t(4)+t(1)*x(3)+x(1)*x(4)+x(1)+x(2)*x(3)*x(4)+x(3)*x(4)+x(3)+x(4)+1]


f = [t[4]+1;
t[3]+1;
t[2]+1;
t[1]+1;
t[3]*t[4]+t[4]+1;
t[2]*t[3]+t[3]+1;
t[1]*t[4]+t[1]+1;
t[1]*t[2]+t[2]+1;
t[2]*t[3]*t[4]+t[3]*t[4]+t[4]+1;
t[1]*t[3]*t[4]+t[1]*t[4]+t[1]+1;
t[1]*t[2]*t[4]+t[1]*t[2]+t[2]+1;
t[1]*t[2]*t[3]+t[2]*t[3]+t[3]+1]


#C3 example
fC = [t[3]+1;
t[2]+1;
t[1]+1;
t[2]*t[3]+t[3]+1;
t[1]*t[2]+t[1]+1;
t[2]^2*t[3]+2*t[2]*t[3]+t[3]+1;
t[1]*t[2]*t[3]+t[1]*t[3]+t[1]+t[3]+1;
t[1]*t[2]^2*t[3]+2*t[1]*t[2]*t[3]+t[1]*t[3]+t[1]+t[2]*t[3]+t[3]+1;
t[1]^2*t[2]^2*t[3]+2*t[1]^2*t[2]*t[3]+t[1]^2*t[3]+t[1]^2+2*t[1]*t[2]*t[3]+2*t[1]*t[3]+2*t[1]+t[3]+1]

MC, FC = unimod_matrix(fC)

#B3 example
fB = [t[3]+1;
t[2]+1;
t[1]+1;
t[2]*t[3]+t[3]+1;
t[1]*t[2]+t[1]+1;
t[2]*t[3]^2+t[3]^2+2*t[3]+1;
t[1]*t[2]*t[3]+t[1]*t[3]+t[1]+t[3]+1;
t[1]*t[2]*t[3]^2+t[1]*t[3]^2+2*t[1]*t[3]+t[1]+t[3]^2+2*t[3]+1;
t[1]*t[2]^2*t[3]^2+2*t[1]*t[2]*t[3]^2+2*t[1]*t[2]*t[3]+t[1]*t[3]^2+2*t[1]*t[3]+t[1]+t[2]*t[3]^2+t[3]^2+2*t[3]+1]

MB, FB = unimod_matrix(fB)
