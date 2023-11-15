function [x, iter] = gradconj(A, b, x, tol, N)
    x0 = x;
    iter = 0;
    r = b - A * x;
    p = r;

    while (iter < N)
        alpha = (r' * r) / (p' * A * p);
        x = x + alpha * p;
        r_old = r;
        r = r - alpha * A * p;

        if (norm(r) < tol)
            return;
        end

        beta = (r' * r) / (r_old' * r_old);
        p = r + beta * p;
        iter += 1;
    end
endfunction
