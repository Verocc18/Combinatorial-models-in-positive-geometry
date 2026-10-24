using Oscar;

b = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
A = [-1 0 0 0 0; -1 -1 0 1 0;  0  1  0  -1 -1; 0  0  0  0  1; 0  1  0  0  0; 0  -1 -1 0  1; -1 0  1  0  -1;  0  0  -1 0  0 ; -1 0  1  1  0; 0  0  0  -1 0];
#A = [1 0 0  0  0  ; 1 1 0  -1 0  ; 0 1 0  -1 -1 ; 0 0 0  0  1  ; 0 1 0  0  0  ;  0 1 1  0  -1 ; 1 0 -1 0  1  ;  0 0 1  0  0  ; 1 0 -1 -1 0  ; 0 0 0  -1 0 ];
#A = [1  0 0 0  0 ;1  1 0 -1 0 ;0  1 0 -1 -1;0  0 0 0  -1;0  1 0 0  0 ;0  1 1 0  -1;-1 0 1 0  -1;0  0 1 0  0 ;-1 0 1 1  0 ;0  0 0 1  0  ];


P=polyhedron(A,b)
f_vec=f_vector(P)
coneGens = rays(P)



#Intersecting quartic hypersurface
R, (x1, x2, x3, x4, x5) = polynomial_ring(QQ, ["x1", "x2", "x3", "x4", "x5"])


S = R[0  x1  -x1-x2+x4  -x5+x2-x4  x5;
  x1  0  x2  -x2+x5-x3  -x5-x1+x3;
  -x1-x2+x4  x2  0  x3  x1-x3-x4;
  -x5+x2-x4  -x2+x5-x3  x3  0  x4;
  x5  -x5-x1+x3  x1-x3-x4  x4  0];
  
f = det(S[1:4,1:4])

N = 10000
pos=0
for i in 1:N
  randomCoefs = [rand(0:100) for i in 1:6]
  norm = sum(randomCoefs)
  randomPosCombRays = sum([(randomCoefs[i]//norm)*coneGens[i] for i in 1:6])
  randomPointInP = [randomPosCombRays[i] for i in 1:5]
  v = evaluate(f,randomPointInP)
  if v >= 0
    pos = pos+1
  end
end
pos

#Checking the facets
Facets = faces(P,4)
coneGens = rays(Facets[2])






#Cyclic polytope
Q = cyclic_polytope(4,6)
visualize(Q)
Qd=polarize(Q)

#Setting x[1]=1
B = A[2:10,2:5]
c = -A[2:10,1]

P=polyhedron(B,c)
F=faces(P,3)
