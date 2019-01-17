function Q = cal_Q(A, d, lambda, U)
    AU = A * U;
    aa = sum(sum(U .* AU));
    UTB = U' * d;
     
    Q = aa - lambda * (UTB' * UTB);
    Q = Q / nnz(A);
end