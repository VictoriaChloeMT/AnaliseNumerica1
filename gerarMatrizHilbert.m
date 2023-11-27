function H = gerarMatrizHilbert(n)
    format long
    H = zeros(n, n);
    for i = 1:n
        for j = 1:n
            H(i, j) = 1 / (i + j - 1);
        end
    end
end

