---- Macaulay2 file for Dihedral sign patterns in M_{0,n} ----

--INSTRUCTIONS
--n and R should be defined globally before using any function from Part 1 and 2
--Sign patterns should be inputed as Hash Tables with keys being ALL 2-subsets of {1..n} and values being 1 or -1
--List of adjacent transposition p are read in order from left to right

--Global information
n = 5;
u = symbol u;
R = QQ[apply(subsets(toList (1..n), 2), S -> u_(toSequence S))];

--PART 1: Consistency of sign patterns
--Checks if the sign pattern is consistent with respect to the extended u-equation given by P (a cyclic partion of {1..n} into 4 intervals)
isConsistentForP = (s,P) ->(
    A := P_0;
    B := P_1;
    C := P_2;
    D := P_3;

    extU1 := 1; 
    extU2 := 1;

    for i in A do(
        for j in C do(
            I := sort {i, j};
            extU1 = extU1*s#I;
        );
    );

    for i in B do(
        for j in D do(
            I := sort {i, j};
            extU2 = extU2*s#I;
        );
    );

    test := extU1 + extU2;

    return not (test == -2);
);
--Creates all possible splittings of 1..n into 4 cyclic intervals
cyclicDecompositions = (n) -> (
    if n < 4 then return {};
    decomps = {};

    -- Offsets from 0 to n-1, where 0 corresponds to element 1
    for a from 1 to n - 3 do (
        for b from a + 1 to n - 2 do (
            for c from b + 1 to n - 1 do (
                for d from c+1 to n do(
                    int1 := toList(a..b-1);
                    int2 := toList(b..c-1);
                    int3 := toList(c..d-1);
                    int4 := toList(d..n);
                    if a > 1 then int4 = flatten (append(int4, toList(1..a-1)));
                    intervals = {int1, int2, int3, int4};
                    decomps = append(decomps, intervals);
                ); 
            );
        );
    );
    return decomps;
);
--Checks if a sign pattern is consistent
isConsistent = (s) -> (
    Part := cyclicDecompositions(n);
    for P in Part do(
        if isConsistentForP(s, P) == false then return false;
    );
    return true;
);
--Creates sign pattern from negative chords in n-gon
signPatternFromList = (L) -> (
   s := new MutableHashTable;
    for i in subsets(toList(1..n),2) do(  
        if isMember(i,L) then s#i = -1
        else s#i = 1;
    );
    return s;
);
--Computes length of the chord i-j in the n-gon labeled by L
chordLength = (i,j,L) -> (
    posi := position(L, k -> k == i);
    posj := position(L, k -> k == j);

    opt1 := (posi-posj)%n;
    opt2 := (posj-posi)%n;
    return min(opt1, opt2);
);
--Generates all consistent sign patterns in an n-gon, N (optional) is the number of negative chords and l the length of the shortest negative chord
allConsistentSignPattern = method();
allConsistentSignPattern(ZZ,ZZ) := (N,l) -> (
    L := toList(1..n);
    posNegs := {};
    for i from 1 to n-1 do(
        for j from i+1 to n do(
            if chordLength(i,j,L) > l-1 then(
                posNegs = append(posNegs, {i,j});
            );
        );
    );
    subsNeg := subsets(posNegs, N);
    consistentSP := {};
    for i in subsNeg do(
        s = signPatternFromList(i);
        if isConsistent(s) then(
            consistentSP = append(consistentSP, s);
        );
    );
    return consistentSP;
);
allConsistentSignPattern(ZZ) := (l) -> (
    consistentSP := {};
    for N from 1 to (binomial(n,2)-n) do(
        consistentSP = join(consistentSP, allConsistentSignPattern(N,l));
    );
    return consistentSP;
);

--PART 2: Dihedral orderings
--Given a dihedral order L computes the corresponding ring of u-variables
ringFromDihedralOrder = (L) -> (
    use R;
    edgesTarget := {};
    for i from 0 to #L-2 do(
        edgesTarget = append(edgesTarget, u_(toSequence sort {L_i, L_(i+1)})-1);
    );
    edgesTarget = append(edgesTarget, u_(toSequence sort {L_0, L_(n-1)})-1);

    I := ideal(edgesTarget);
    T := frac(R/I);
    T#"Dihedral Order"=L;

    return T;
);
--Applies the adjacent transposition based at k to the dihedral order L
newDihedralOrder = (L,k) -> (
    newL := {};
    posk := position(L, i -> i == k);
    if posk == n-1 then(
        newL = switch(0,n-1,L);
    ) else(
        newL = switch(posk, posk+1, L);
    );
    return newL;
);
--Applies the sequence of adjacent transposition p to the dihedral order L
newDihedralOrderFromP = (L,p) -> (
    tempL := L;
    for k in p do(
        tempL = newDihedralOrder(tempL, k);
    );
    return tempL;
);
--Computes the image of the monomial map D to T corresponding to the adjacent transposition at k
imAdjTranspositionMap = (D, T, k) ->(
    L := T#"Dihedral Order";
    newL := D#"Dihedral Order";

    posk := position(L, i -> i == k);
    poskp1 := -1;
    poskm1 := -1;
    if posk == n-1 then(
        poskp1 = 0;
        poskm1 = n-2;
    ) else(
        poskp1 = posk+1;
        if posk == 0 then poskm1 = n-1 else poskm1 = posk-1
    );

    images = {};
    for x in gens(D) do(
        if x == 1 then images = append(images, 1)
        else(
            use T;
            (i, j) = last baseName x;
            if set(i,j) == set(k, L_poskm1) then(
                prodL := delete(k,L);
                prodL = delete(L_poskm1, prodL);
                prodL = delete(L_poskp1, prodL);
                images = append(images, -u_(toSequence sort{L_poskm1, L_poskp1})*product(apply(prodL, l -> (u_(toSequence sort{l, k}))^(-1))));
            ) else if isMember(k,{i,j}) then(
                J := (delete(k, {i,j}))_0;
                images = append(images, u_(i,j)*u_(toSequence sort{J, L_poskp1}));
            ) else if isMember(L_poskm1, {i,j}) then(
                J := (delete(L_poskm1,{i,j}))_0;
                images = append(images, u_(i,j)*u_(toSequence sort{J, k}));
            ) else if isMember(L_poskp1, {i,j}) then(
                J := (delete(L_poskp1, {i,j}))_0;
                images = append(images, u_(toSequence sort{J, k})^(-1));
            ) else(
                images = append(images, u_(i,j));
            );
        );
    );

    return images;
);
--Computes the mononimal map corresponding to applying the sequence of adjacent transpositions p to the dihedral order L
compMapFromPerm = (L, p) -> (
    T := ringFromDihedralOrder(L);
    alpha := newDihedralOrder(L, p_0);
    D := ringFromDihedralOrder(alpha);
    im := imAdjTranspositionMap(D,T,p_0);
    f := map(T, D, im);
    for k from 1 to #p-1 do(
        T = D;
        alpha = newDihedralOrder(alpha, p_k);
        D = ringFromDihedralOrder(alpha);
        im = imAdjTranspositionMap(D,T,p_k);
        g := map(T, D, im);
        f = f*g;
    );
    return f;
);
--Computes the sign pattern after applying the monomial map of p to the sign pattern s
signPatternFromPerm = (L, p, s) ->(
    f := compMapFromPerm(L, p);
    D := source(f);
    T := target(f);
    Lists := {};
    for x in gens(R) do (
        (i,j) = last baseName x;
        I := sort {i,j};
        Lists = append(Lists, s#I);
    );

    m := map(QQ, T, Lists);

    return (f,apply(gens(D), x -> (m*f)(x)));
);
--Compute the shortes negative chord of the sign pattern s in the n-gon labeled by L
shortNegChord = (L,s) -> (
    negs := select(subsets(toList(1..n),2), i -> s#i == -1);
    if #negs == 0 then return {}
    else(
        lengths := apply(negs, j -> chordLength(j_0,j_1,L));
        minLength := minPosition(lengths);
        return negs_(minLength);
    );
);
--Computes the sequence of adjacent tranposition that flips the sign of the chord i-j
flipChord = (L,s,i,j) -> (
    l := chordLength(i,j,L);
    posi := position(L, k -> k == i);
    posj := position(L, k -> k == j);
    p := {};
    if l == ((posj-posi)%n) then (
        if posi < n-1 then p = toList(l-1:L_(posi+1))
        else p = toList(l-1:L_0);
        temp := posj;
        for k from 1 to l-2 do(
            if temp > 0 then(
                p = append(p, L_(temp-1));
                temp = temp-1;
            )
            else (
                p = append(p, L_(n-1));
                temp = n-1;
            );
        );
    )
    else (
        if posj < n-1 then p = toList(l-1:L_(posj+1))
        else p = toList(l-1:L_0);
        temp := posi;
        for k from 1 to l-2 do(
            if temp > 0 then(
                p = append(p, L_(temp-1));
                temp = temp-1;
            )
            else (
                p = append(p, L_(n-1));
                temp = n-1;
            );
        );
    );

    return p; 
);
--Computes the dihedral ordering for a consistent sign pattern
dihedralOrderFromSignPattern = (s) -> (
  L := toList(1..n);
  shortNeg := shortNegChord(L,s);
  t := s;
    while shortNeg != {} do(
        i := shortNeg_0;
        j := shortNeg_1;
        p := flipChord(L,t,i,j);
        (f,news) = signPatternFromPerm(L, p, t);
        negs := (gens(source(f)))_(positions(news, i -> i == -1));
        t = signPatternFromList(apply(negs, x -> toList(last baseName x)));
        L = newDihedralOrderFromP(L,p);
        shortNeg = shortNegChord(L,t);
    );
    return L;
);

--PART 3: Partial compactification
--Ideal of u-relation in dihedral order 1..n
uRelationsIdeal = (n) -> (
    u := symbol u;
    R := QQ[apply(subsets(toList (1..n), 2), S -> u_(toSequence S))];
    I := ideal for x in gens(R) list(

    (i, j) = last baseName x;

    prod = product flatten for k from i + 1 to j-1 list for l in select(1..n, a -> a > j or a < i) list u_(toSequence sort {k, l});
    if prod == 1_R then continue else -1 + u_(i,j) + prod
    );
    return I;
);



