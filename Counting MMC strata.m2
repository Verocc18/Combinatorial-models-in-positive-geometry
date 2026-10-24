--Stirling number of second kind
s := (n,k) -> (
    if n < k then(
        return 0;
    )
    else (
        rta := {};
        for i from 0 to k do(
            rta = append(rta, (-1)^(k-i)*binomial(k,i)*(i)^n);
        );
        return sum(rta)/(k!);
    );
);

N:= (r,m,n,p) -> (
    if r>2 and r<m then(
        return s(n,m)-s(p+1,m)-s(n-p+1,m);
    )
    else if r>1 and r==m then(
        return (m!)*s(p,m)*s(n-p,m);
    )
    else (
        return 0;
    );
);

n=5;
c = binomial(n,2) - n;

--Fixed p
MMCs = mutableMatrix(QQ,c,n-2);
p = 2;
for r from 2 to n-1 do(
    for d from 1 to c do(
        mList := {};
        M := 0;
        if r==2 then M=2 
        else M=floor((d+r+binomial(r,2))/(r-1));
        for m from r to M do(
            l := -d-binomial(r,2)+(m-1)*(r-1)+n-m-1;
            for i from 0 to l do(
                mList = append(mList, binomial(p,i)*binomial(n-p,l-i)*N(r,m,n-l,p-i));
            );
        );
        MMCs_(d-1,r-2) = sum(mList);
    );
)

--No fixed p
MMC = mutableMatrix(QQ,c,n-2);
for r from 2 to n-1 do(
    for d from 1 to c do(
        mList := {};
        M := 0;
        if r==2 then M=2 
        else M=floor((d+r+binomial(r,2))/(r-1));

        for m from r to M do(
            l := -d-binomial(r,2)+(m-1)*(r-1)+n-m-1;
            for p from 2 to (n-2) do(
                for i from 0 to l do(
                    mList = append(mList, binomial(n-1,p-1)*binomial(p,i)*binomial(n-p,l-i)*N(r,m,n-l,p-i)*(1/2^l));
                );
            );
        );
        MMC_(d-1,r-2) = sum(mList);
    );
)