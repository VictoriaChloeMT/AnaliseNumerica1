function [x, iter] = gradiente(A, b, x0, tol, max_iter)
    x = x0;
    iter = 0;
    while (iter <= max_iter)
        r = b - A * x;
        S = (r' * r) / (r' * A * r);
        x = x0 + S * r;
        if (norm(r) < tol)
            return;
        end
        x0 = x;
        iter = iter + 1;
    end
endfunction
