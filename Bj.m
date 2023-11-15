function Bj = Bj(A)
    n = size(A, 1);
    D = diag(diag(A));

    Bj = eye(n) - inv(D) * A;
endfunction
