--Computational proof of inconsistent sign patterns
--Run function allConsistentSignPatterns from DihedralSignPatternsM0n.m2 first

--Octagon
n = 8;
candsMi = {};
candsui1 = {};
candsuij_1 = {};
for m from 2 to 16 do(
    --Computes consisten sign patterns
    cp:=allConsistentSignPattern(m,2);
    --Sign pattern conditions appearing in Table 1
    test1:=select(cp, s -> s#{1,3} == 1 and s#{1,4} == 1 and s#{1,6} == 1 and s#{2,4} == 1 and s#{2,6} == -1 and s#{2,8} == 1 and s#{3,6} == 1 and s#{3,8} == 1 and s#{4,6} == 1 and s#{4,8} == -1 and s#{6,8} == 1);
    test2:=select(cp, s -> s#{1,3} == 1 and s#{1,4} == 1 and s#{1,6} == -1 and s#{2,4} == 1 and s#{2,6} == 1 and s#{2,8} == 1 and s#{3,6} == 1 and s#{3,8} == 1 and s#{4,6} == 1 and s#{4,8} == -1 and s#{6,8} == 1);
    test3:=select(cp, s -> s#{1,3} == 1 and s#{1,4} == 1 and s#{1,6} == 1 and s#{2,4} == 1 and s#{2,6} == 1 and s#{2,8} == 1 and s#{3,6} == -1 and s#{3,8} == 1 and s#{4,6} == 1 and s#{4,8} == -1 and s#{6,8} == 1);


    if #test1>0 then(
        candsMi = append(candsMi, test);
    );
    if #test2>0 then(
        candsui1 = append(candsui1, test);
    );
    if #test3>0 then(
        candsuij_1 = append(candsuij_1, test);
    );
);
--Checks nothing was added
candsMi == {}
candsui1 == {}
candsuij_1 == {}

--Heptagon
n=7;
cp=allConsistentSignPattern(2);

test1=select(cp, s -> s#{1,3} == 1 and s#{1,4} == 1 and s#{1,5} == -1 and s#{2,4} == 1 and s#{2,5} == 1 and s#{2,7} == 1 and s#{3,5} == 1 and s#{3,7} == 1 and s#{4,7} == -1 and s#{5,7} == 1);
test2=select(cp, s -> s#{1,3} == 1 and s#{1,4} == 1 and s#{1,6} == 1 and s#{2,4} == 1 and s#{2,6} == 1 and s#{2,7} == 1 and s#{3,6} == -1 and s#{3,7} == 1 and s#{4,7} == -1 and s#{4,6} == 1);
