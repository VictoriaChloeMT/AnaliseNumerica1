clear; n = 100; H = gerarMatrizHilbert(n); b=ones(n,1);

x0 = zeros(n,1); omega1 = 0.01; omega2 = 0.02; tolerance = 1e-7; max_iterations = 10;
X = linsolve(H, b);

%% LU
[L, U] = lu(H);

y = L \ b; x1 = U \ y;

%% Cholesky
L = cholesky(H);

y = L \ b; x2 = L' \ y;

%% JOR
[x3, iter] = jor(H, b, x0, omega1, tolerance, max_iterations);

%% SOR
[x4, iter] = sor(H, b, x0, omega2, tolerance, max_iterations);

%% Max descida
[x5, iter] = gradiente(H, b, x0, tolerance, max_iterations);

%% Gradientes Conjugados
[x6, iter] = gradconj(H, b, x0, tolerance, max_iterations);

%% Erros
Ea_LU = abs(X - x1); Ea_Cholesky = abs(X - x2); Ea_JOR = abs(X - x3); Ea_SOR = abs(X - x4);Ea_MaxDescida = abs(X - x5); Ea_GradConj = abs(X - x6);


Er_LU = abs(Ea_LU ./ abs(X)); Er_Cholesky = abs(Ea_Cholesky ./ abs(X)); Er_JOR = abs(Ea_JOR ./ abs(X)); Er_SOR = abs(Ea_SOR ./ abs(X)); Er_MaxDescida = abs(Ea_MaxDescida ./ abs(X)); Er_GradConj = abs(Ea_GradConj ./ abs(X));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%Erros relativos%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LU
figure;
plot(Er_LU);
title('Decomposição LU');
xlabel('n');
ylabel('Erro');

% Cholesky
figure;
plot(Er_Cholesky);
title('Decomposição Cholesky');
xlabel('n');
ylabel('Erro');

% JOR
figure;
plot(Er_JOR);
title('JOR');
xlabel('n');
ylabel('Erro');

% SOR
figure;
plot(Er_SOR);
title('SOR');
xlabel('n');
ylabel('Erro');

% Max Descida
figure;
plot(Er_MaxDescida);
title('Máxima Descida');
xlabel('n');
ylabel('Erro');

% Gradientes Conjugados
figure;
plot(Er_GradConj);
title('Gradientes Conjugados');
xlabel('n');
ylabel('Erro');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%Erros absolutos%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LU
figure;
plot(Ea_LU);
title('Decomposição LU');
xlabel('n');
ylabel('Erro');

% Cholesky
figure;
plot(Ea_Cholesky);
title('Decomposição Cholesky');
xlabel('n');
ylabel('Erro');

% JOR
figure;
plot(Ea_JOR);
title('JOR');
xlabel('n');
ylabel('Erro');

% SOR
figure;
plot(Ea_SOR);
title('SOR');
xlabel('n');
ylabel('Erro');

% Max Descida
figure;
plot(Ea_MaxDescida);
title('Máxima Descida');
xlabel('n');
ylabel('Erro');

% Gradientes Conjugados
figure;
plot(Ea_GradConj);
title('Gradientes Conjugados');
xlabel('n');
ylabel('Erro');


% Comparação dos resultados
figure;
plot(x1, 'DisplayName', 'LU');
hold on;
plot(x2, 'DisplayName', 'Cholesky');
title('Resultados dos métodos LU e Cholesky');
xlabel('n');
ylabel('Resultado');
legend;
hold off;

% Comparação dos resultados
figure;
plot(X, 'DisplayName', 'Exato');
hold on;
plot(x1, 'DisplayName', 'LU');
title('Resultados exatos VS LU');
xlabel('n');
ylabel('Resultado');
legend;
hold off;

% Comparação dos resultados
figure;
plot(X, 'DisplayName', 'Exato');
hold on;
plot(x2, 'DisplayName', 'Cholesky');
title('Resultados exatos VS Cholesky');
xlabel('n');
ylabel('Resultado');
legend;
hold off;

% Comparação dos resultados
figure;
plot(X, 'DisplayName', 'Exato');
hold on;
plot(x3, 'DisplayName', 'JOR');
title('Resultados exatos VS JOR');
xlabel('n');
ylabel('Resultado');
legend;
hold off;

% Comparação dos resultados
figure;
plot(X, 'DisplayName', 'Exato');
hold on;
plot(x4, 'DisplayName', 'SOR');
title('Resultados exatos VS SOR');
xlabel('n');
ylabel('Resultado');
legend;
hold off;

% Comparação dos resultados
figure;
plot(X, 'DisplayName', 'Exato');
hold on;
plot(x5, 'DisplayName', 'Máx Descida');
title('Resultados exatos VS Máx Descida');
xlabel('n');
ylabel('Resultado');
legend;
hold off;

% Comparação dos resultados
figure;
plot(X, 'DisplayName', 'Exato');
hold on;
plot(x6, 'DisplayName', 'Grad Conj');
title('Resultados exatos VS Grad Conj');
xlabel('n');
ylabel('Resultado');
legend;
hold off;

