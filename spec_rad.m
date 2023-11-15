function sr = spec_rad(A)
    eign = eig(A);
    sr = max(abs(eign));
end

