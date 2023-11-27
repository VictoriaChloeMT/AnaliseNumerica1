function Bs = Bs(A)
    n = size(A, 1);
    D = diag(diag(A));

    Bs = inv(D) * A;
endfunction
