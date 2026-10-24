needsPackage "WeylGroups";
needsPackage "SpechtModule";
needsPackage "Matroids";



--- Preliminary functions ---


--- Returns the matrix M without the ith and jth columns and rows.
minorMatrix = (M, i, j) -> (
    N := numRows M;
    if N != numColumns(M) then error "Matrix must be square";

    rows = sort(toList((set(toList(0..N-1)) - set({i,j}))));
    cols = sort(toList(set(toList(0..N-1)) - set({i,j})));

    return M^rows_cols;
);

--- Returns the pfaffian of a skew symmetric matrix
pfaffian = (M) -> (
    --Ensure the matrix is square and skew-symmetric
    N := numRows M;
    if N != numColumns(M) then error "Matrix must be square";
    if N % 2 != 0 then error "Matrix must have even dimensions";
    if M + transpose(M) != 0 then error "Matrix must be skew symmetric";

    -- Base cases
    if N == 0 then return 1;
    if N == 2 then return M_1_0;

    -- Recursively compute the Pfaffian
    return sum( apply( toList(1..N-1),    j -> (-1)^(j) * M_0_j * pfaffian( minorMatrix(M, 0, j) ) ) )
);

-- SO(2n) pinning in Section 2
phi = (n,i, a,b,c,d)->(

        R := ring(b);

        g := mutableMatrix(id_(R^(2*n)));


        if i < n then(
            g_(i-1,i-1) = a;
            g_(i-1,i) = b;
            g_(i,i-1) = c;
            g_(i,i) = d;


            g_(n+i-1, n+i-1) =  d;
            g_(n+i-1, n+i)   = -c;
            g_(n+i  , n+i-1) = -b;
            g_(n+i  , n+i)   =  a;
        );

        
        if i == n then(
            
            g_(n-2, n-2) = a;
            g_(n-2, 2*n-1) = b;
            g_(2*n-1,n-2) = c;
            g_(2*n-1,2*n-1) = d;


            g_(n-1, n-1) =  a;
            g_(n-1,2*n-2)   = -b;
            g_(2*n-2  , n-1) = -c;
            g_(2*n-2  , 2*n-2)   =  d;
        );
        

        return g;

);

-- Involution giving the Gale order in Section 4
invol = (n,i)->(

    if i > n then return n + 1 + (2*n - i);

    if i<= n then return i

);

-- Lowers the permutation u as in Definition 4.14
lowerUinB = (n, u,  B)->(

    tempu := new MutableList from u;

    orderingList := join( toList(1..n) , reverse(toList(n+1..2*n)) );

    itermediateMonomials := {};
    
    for i from 0 to #tempu - 1 do(

        posi := position(orderingList, x-> x== u_i);

        for j from 0 to posi do (
            tempu#i = orderingList_j;

            if (set(B)) #? (sort(  toList(tempu) ) ) then ( break; );

        );

        itermediateMonomials = append(itermediateMonomials,  apply(toList(tempu), y -> y-1 ) );
    );

    return (toList(tempu), itermediateMonomials);

);






--- Main functions ---
-- Returns the minors of A in the positivity test
positivityTest = (A) -> (
    n := numcols(A);
    rta := {};

    for k from 1 to n-1 do(
        for j from 1 to k do(
            I := join(toList(0..(n-k-2)),toList((n-k+j-1)..(n-1)));
            J := toList(0..(n-j-1));
            rta = append(rta, ((sort(apply(I,s -> s+1)),apply(J, s -> s+1)), (-1)^(j*k)*determinant(A_J^I, Strategy=>Cofactor)));
        );

    );

    return rta;

);

-- Checks if a skew-symmetric matrix A over QQ is totally positive or not
isTotallyPositive = (A) -> (
    if not (-A == transpose(A)) then return "The matrix is not skew-symmetric or it does not have rational coefficients.";
    
    PT := positivityTest(A);
    
    for m in PT do(
        if m_1 < 0 then return (false,m);
    );

    return true;

);

-- Generates a totally positive skew-symmetric matrix of size n in variables t_i
totallyPositiveSkewMatrix = (n) -> (
    N := binomial(n,2);

    R := QQ[t_1..t_N];


    Dn := rootSystemD(n);
    w0 := longestWeylGroupElement(Dn);
    L := reducedDecomposition(w0);
    minCosL := L_(toList(N..2*N-1));
    

    M := product apply(#minCosL, i -> phi(n,minCosL_i,substitute(1,R),t_(i+1),0,1));

    X := M^(toList(0..(n-1)));

    apply(toList(1..(n-1)), i -> apply(i, j -> rowAdd(X,j,-(X^{j})_(0,i),i)));

    A := matrix (X);

    return A_(toList(n..2*n-1));
);

--Returns leading terms of non-negativity test 
nonNegativityTest = (A) -> (
    n := numcols(A);
    R := ring(A);

    Dn := rootSystemD(n);
    w0 := longestWeylGroupElement(Dn);
    L := reducedDecomposition(w0);

    S := R[e];
    z := matrix(product apply(#L, i -> phi(n,L_i,substitute(1,S),e,0,1)));
    Z := transpose(z)*z;

    B := substitute((id_(ZZ^n) | A),S);
    X := B*Z;

    minors := {};
    for k from 1 to n-1 do(
        for j from 1 to k do(
            I := join(toList(0..(n-k-2)),toList((n-k+j-1)..(n-1)));
            J := toList(0..(n-j-1));
            K := join(toList((n-k-1)..(n-k+j-2)),toList(n..(2*n-j-1)));
            minors = append(minors, ((sort(apply(I,s -> s+1)),apply(J, s -> s+1)), (-1)^(j*(n-1))*determinant(X_K, Strategy=>Cofactor)));
        );
    );

    l := apply(minors, m -> (m_0,substitute(last (flatten entries ((coefficients(m_1))_1)),R)));

    return l;

);

-- Checks if a skew-symmetric matrix A over QQ is totally non-negative or not
isTotallyNonNegative = (A) -> (
    if not (A == -transpose(A)) then return "The matrix is not skew-symmetric or it does not have rational coefficients.";
    
    NNT := nonNegativityTest(A);
    
    for m in NNT do(
        if m_1 < 0 then return (false,m);
    );

    return true;
    
);

-- Given two signed permutation v,w written as words in simple reflections it generates a totally positive skew-symmetric matrix of size n in variables t_i
totallyNonNegativeSkewMatrix = (n,v,w) -> (

    N := #select(w, i-> i!=0) - #(select(v, i -> i != 0));

    R := QQ[t_1..t_N];

    M := mutableMatrix id_(R^(2*n));
    
    j:= 1;
    
    for i from 0 to #w-1 do(

        if v_i == 0 and w_i != 0 then(
            M = M * phi(n,w_i,substitute(1,R),t_j,0,1);
            j = j+1;
        );

        if v_i != 0 and w_i != 0 then (
            M = M*phi(n,v_i,substitute(0,R),substitute(1,R),-1,0);
        );
    );

    X := matrix(M^(toList(0..(n-1))));
    X = (inverse(X_(toList(0..n-1)))*X)_(toList(n..(2*n-1)));

    return X;

);

-- Identifies the permutations corresponding to the cell of a totally non-negative skew-symmetric matrix A over QQ
nonNegativeRichardsonCell = (A) -> (
    n := numcols(A);
    X := (id_(ZZ^n) | A);
    M := matroid(X);
    B := bases(M);
    B = apply(B, b-> apply(sort toList(b), x->x+1));

    Bmax := apply( max( apply(B, b -> apply(b, x -> invol(n,x)) ) ) , x->invol(n,x));

    winv := apply(sort apply(Bmax, b -> invol(n,b)), x -> invol(n,x));
    (v, intermidiate) := lowerUinB(n,winv,B);

    w := {};
    for i from 1 to n do(
        if isMember(i,winv) then(
            w = append(w,position(winv, x->x==i)+1);
        )
        else (
            temp:= position(winv, x->x==i+n)+1;
            w = append(w, n+temp);
        );
    );
    return (v,w);
    
);

-- Returns all Pfaffians of a skew-symmetric matrix A
Pfaffians = (A)->(
    n := numcols(A);
    S := ring(A);
    evenNum := 2*toList(0..n//2);
    temp := apply(evenNum,  k-> subsets(toList(0..n-1), k));
    EvenSubsets := flatten temp;

    f := symbol f;
    PF := QQ[apply(EvenSubsets, I->f_(toSequence(I)))];


    change := map(S,PF, apply(EvenSubsets, I->pfaffian(A^I_I)));
    
    changedPfaffians := {};

    VarPF := flatten entries vars(PF);

    for i from 0 to #VarPF-1 do(
        I := EvenSubsets_i;
        changedPfaffians = append(changedPfaffians, change(VarPF_i));
    );

    return changedPfaffians;
);



--- Useful functions for examples ---

-- Gives the word of a minimal coset representative in W(Dn)^[n-1] corresponding to the even subset I
wordFromEvenSubset = (n,I) ->(

    J := sort(I);

    S := {}; 

    i1 := 0;
    i2 := 0;



    for i from 0 to (#J//2 - 1) do(

        i1 =  J_(2*i);
        i2 =  J_(2*i+1);

        S = append(S, join( {n}  , (reverse(toList(i1..n-2))), i1-2*i-1:0 ,  (reverse(toList(i2..n-1))), i2-2*i-2:0  ));

    );
    
    tempR := join( join(toSequence(S)));
    
    return join( tempR, binomial(n,2) - #tempR:0);

);

-- Gives the one line notation permutation of the word in simple reflections w
permutationFromWord = (n,w) ->(

    perm := new MutableList from toList(1..n);

    for j from 0 to #w-1 do(

        l := #w - j - 1;

        if 0 < w_l and  w_l < n then(

            tempx := perm # (w_l -1);
            tempy := perm # (w_l );

            perm # (w_l -1) = tempy;
            perm # (w_l) = tempx;
        );

        if w_l == n then(

            tempx := perm # (w_l -2);
            tempy := perm # (w_l -1);

            perm # (w_l - 2) = ( (tempy + n - 1) % (2*n) ) + 1;
            perm # (w_l - 1) = ((tempx + n - 1) % (2*n) ) + 1;
        );


    );

    return toList(perm)

);

-- Given a word in simple reflection w it generates a random distinguished subexpression of w
randomDistinguishedSubExp = (n,w)->(

    Dn := rootSystemD(n);
    N := #w;


    r := 0;


    redWord := {};
    tempv := {};

    for i from 0 to  #w-1 do(
        r = random(2);
        if (r == 0) then (tempv = append(tempv, w_i));
    );


    tempv = select(tempv, x->x!=0);

    v := reduce(Dn,tempv);

    for k from 0 to N-1 do(
        
        if w_k != 0 then(
           
           
            siv := reduce(Dn,join({w_k},tempv));

            if coxeterLength(siv) < coxeterLength(v) then(
                redWord=append(redWord, w_k);
                v = siv;
                tempv = reducedDecomposition(v);
            )

            else (
                redWord=append(redWord, 0);
            );

        )
        else (
            redWord=append(redWord, 0);
        );

    );

    return redWord;
    
);













