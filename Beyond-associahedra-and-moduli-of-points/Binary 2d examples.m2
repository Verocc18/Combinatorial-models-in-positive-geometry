--Pentagon
R = QQ[u_1..u_5]

I = ideal(u_1 + u_3*u_4-1, u_2 + u_4*u_5-1,u_3 + u_1*u_5 -1, u_4 + u_1*u_2-1, u_5 + u_2*u_3-1)
dim I == 2
varsR = gens(R)

for u in varsR do(
    print (dim(ideal(u) + I) == 1);
);

subs = subsets(5,2)
edges = {{0,1},{1,2},{2,3},{3,4},{0,4}}

for s in subs do(
    d = dim(ideal(varsR_s) + I);
    if isMember(s, edges) then(
        print(d==0);
    ) else (
        print(d==-1);
    );
);



--Hexagon
R = QQ[u_1..u_6]

--Three ideals to check
I = ideal(u_1 + u_3*u_4*u_5-1, u_2 + u_4*u_5*u_6-1,u_3 + u_1*u_5*u_6 -1, u_4 + u_1*u_2*u_6-1, u_5 + u_1*u_2*u_3-1, u_6 + u_2*u_3*u_4-1)
I = ideal(u_1 + u_3*u_4*u_5-1, u_2 + u_4*u_5*u_6-1,u_3 + u_1*u_5*u_6^2 -1, u_4 + u_1*u_2*u_6-1, u_5 + u_1*u_2*u_3-1, u_6 + u_2*u_3^2*u_4-1)
I = ideal(u_1 + u_3*u_4*u_5-1, u_2 + u_4*u_5^2*u_6-1,u_3 + u_1*u_5*u_6 -1, u_4 + u_1^2*u_2*u_6-1, u_5 + u_1*u_2*u_3-1, u_6 + u_2*u_3^2*u_4-1)



--Checks to see if it is a binary geometry
dim I == 2

varsR = gens(R)

for u in varsR do(
    print (dim(ideal(u) + I) == 1);
);

subs = subsets(6,2)
edges = {{0,1},{1,2},{2,3},{3,4},{4,5},{0,5}}

for s in subs do(
    d = dim(ideal(varsR_s) + I);
    if isMember(s, edges) then(
        print(d==0);
    ) else (
        print(d==-1);
    );
);