function [x, iter] = jor(A, b, x0, omega, tol, max_iter)
  [n, n] = size(A);
  iter = 0;
  r = b - A * x0;
  r0 = norm(r);
  e = norm(r);
  x = x0;

  while e > tol && iter < max_iter
    iter = iter + 1;

    for i = 1:n
      s = 0;
      for j = 1:i-1
        s = s + A(i, j) * x(j);
      end

      for j = i+1:n
        s = s + A(i, j) * x(j);
      end

      x(i) = omega * (b(i) - s) / A(i, i) + (1 - omega) * x(i);
    end

    r = b - A * x;
    e = norm(r) / r0;
  end
endfunction


