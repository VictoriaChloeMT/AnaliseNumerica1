function Bsor = Bsor(A, w)
    n = size(A, 1);
    D = diag(diag(A));

    lower = tril(A) - D;
    upper = triu(A) - D;

    Bsor = inv(eye(n) - w * inv(D) * lower) * ((1 - w) * eye(n) + w * inv(D) * upper);
endfunction
