function [x, iter] = sor(A, b, x0, omega, TOL, max_iter)
    n = length(b);
    x = x0;

    for iter = 1:max_iter
        for i = 1:n
            sum1 = A(i, 1:i-1) * x(1:i-1);
            sum2 = A(i, i+1:end) * x0(i+1:end);
            x(i) = (1 - omega) * x0(i) + (omega / A(i, i)) * (b(i) - sum1 - sum2);
        end

        if norm(x - x0, 'inf') < TOL
            return;
        end

        x0 = x;
    end
endfunction
