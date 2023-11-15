clear
n = 5;
H = gerarMatrizHilbert(n);
b=ones(n,1);

x0 = zeros(n,1);
omega1 = 1.8;
omega2 = 0.1;  % Fator de relaxação
tolerance = 1e-7;
max_iterations = 100;

X = linsolve(H, b);

%% LU

[L, U] = lu(H);

y = L \ b;

x1 = U \ y;

%% Cholesky

L = cholesky(H);

y = L \ b;

x2 = L' \ y;

%% JOR
[x3, iter] = jor(H, b, x0, omega1, tolerance, max_iterations);

%% SOR
[x4, iter] = sor(H, b, x0, omega2, tolerance, max_iterations);

%% Max descida
[x5, iter] = gradiente(H, b, x0, tolerance, max_iterations);

%% Gradientes Conjugados
[x6, iter] = gradconj(H, b, x0, tolerance, max_iterations);

Bjo = Bjor(H, omega1);

srJO = spec_rad(Bjo)

Bsor = Bsor(H, omega2);

srSOR = spec_rad(Bsor)

Bj= Bj(H);

srJ = spec_rad(Bj)

cond_number = cond(H);

Ea_LU = abs(X - x1);
Ea_Cholesky = abs(X - x2);
Ea_JOR = abs(X - x3);
Ea_SOR = abs(X - x4);
Ea_MaxDescida = abs(X - x5);
Ea_GradConj = abs(X - x6);


Er_LU = abs(Ea_LU ./ abs(X));
Er_Cholesky = abs(Ea_Cholesky ./ abs(X));
Er_JOR = abs(Ea_JOR ./ abs(X));
Er_SOR = abs(Ea_SOR ./ abs(X));
Er_MaxDescida = abs(Ea_MaxDescida ./ abs(X));
Er_GradConj = abs(Ea_GradConj ./ abs(X));

figure;

subplot(2, 1, 1);
semilogy(1:length(Ea_LU), Ea_LU, '-o', 1:length(Ea_Cholesky), Ea_Cholesky, '-o');
title('Erros Absolutos');
legend('LU', 'Cholesky', 'JOR', 'SOR', 'Max Descida', 'Gradientes Conjugados');
ylabel('Erro Absoluto');

subplot(2, 1, 2);
semilogy(1:length(Er_LU), Er_LU, '-o', 1:length(Er_Cholesky), Er_Cholesky, '-o')
title('Erros Relativos');
legend('LU', 'Cholesky');
ylabel('Erro Relativo');

figure;
subplot(2, 1, 1);
semilogy(1:length(Ea_JOR), Ea_JOR, '-o', 1:length(Ea_SOR), Ea_SOR, '-o', 1:length(Ea_MaxDescida), Ea_MaxDescida, '-o', 1:length(Ea_GradConj), Ea_GradConj, '-o');
title('Erros Absolutos');
legend('JOR', 'SOR', 'Max Descida', 'Gradientes Conjugados');
ylabel('Erro Absoluto');

subplot(2, 1, 2);
semilogy(1:length(Er_JOR), Er_JOR, '-o', 1:length(Er_SOR), Er_SOR, '-o', 1:length(Er_MaxDescida), Er_MaxDescida, '-o', 1:length(Er_GradConj), Er_GradConj, '-o');
title('Erros Relativos');
legend('JOR', 'SOR', 'Max Descida', 'Gradientes Conjugados');
ylabel('Erro Relativo');


