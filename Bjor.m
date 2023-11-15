function Bjo = Bjor(A, w)
  D = diag(diag(A));
  E = A - D;
  Bj = inv(D) * E;
  Bjo = w * Bj + (1 - w) * eye(size(A, 1));
endfunction
