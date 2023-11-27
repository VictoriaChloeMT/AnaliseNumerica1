clear
n = 100;
H = gerarMatrizHilbert(n);
b=ones(n,1);

x0 = zeros(n,1);
omega1 = 0.02 ;
omega2 = 1.25;
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

Bsor_func = Bsor(H, omega2);

srSOR = spec_rad(Bsor_func)

Bj = Bj(H);

srJ = spec_rad(Bj)

Bs = Bs(H);

srS = spec_rad(Bs)
interv = 2/srS;

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


graphics_toolkit('gnuplot');
omega_1 = 0:0.01:interv;
radios_espectrais_JOR = zeros(size(omega_1));

for i = 1:length(omega_1)
    Bjo = Bjor(H, omega_1(i));
    radios_espectrais_JOR(i) = spec_rad(Bjo);
end

figure; % Primeira figura
plot(omega_1, radios_espectrais_JOR, 'LineWidth', 2);
xlabel('Omega');
ylabel('Raio Espectral');
legend('JOR');
grid on;

omega_2 = 0:0.01:2;
radios_espectrais_SOR = zeros(size(omega_2));

for i = 1:length(omega_2)
    Bsor_func = Bsor(H, omega_2(i));
    radios_espectrais_SOR(i) = spec_rad(Bsor_func);
end

figure; % Segunda figura
plot(omega_2, radios_espectrais_SOR, 'LineWidth', 2);
xlabel('Omega');
ylabel('Raio Espectral');
legend('SOR');
grid on;

% Gráfico de Convergência
figure;
plot(1:length(Ea_JOR), Ea_JOR, 'o-', 1:length(Ea_SOR), Ea_SOR, 's-', 1:length(Ea_MaxDescida), Ea_MaxDescida, '^-', 1:length(Ea_GradConj), Ea_GradConj, 'd-');
legend('JOR', 'SOR', 'Max Descida', 'Gradientes Conjugados');
xlabel('Número de Iterações');
ylabel('Erro Absoluto');
title('Gráfico de Convergência');
grid on;

% Gráfico de Dispersão de Resíduos
figure;
scatter(1:length(Ea_JOR), Ea_JOR, 'o', 'filled');
hold on;
scatter(1:length(Ea_SOR), Ea_SOR, 's', 'filled');
scatter(1:length(Ea_MaxDescida), Ea_MaxDescida, '^', 'filled');
scatter(1:length(Ea_GradConj), Ea_GradConj, 'd', 'filled');
legend('JOR', 'SOR', 'Max Descida', 'Gradientes Conjugados');
xlabel('Iteração');
ylabel('Resíduos');
title('Gráfico de Dispersão de Resíduos');
grid on;

% Histograma de Erros
figure;
hist([Ea_JOR, Ea_SOR, Ea_MaxDescida, Ea_GradConj], 20);
xlabel('Erro Absoluto');
ylabel('Frequência');
title('Histograma de Erros');
legend('JOR', 'SOR', 'Max Descida', 'Gradientes Conjugados');
grid on;

% Erro Absoluto vs Erro Relativo
figure;
scatter(Er_JOR, Ea_JOR, 'o', 'filled');
hold on;
scatter(Er_SOR, Ea_SOR, 's', 'filled');
scatter(Er_MaxDescida, Ea_MaxDescida, '^', 'filled');
scatter(Er_GradConj, Ea_GradConj, 'd', 'filled');
legend('JOR', 'SOR', 'Max Descida', 'Gradientes Conjugados');
xlabel('Erro Relativo');
ylabel('Erro Absoluto');
title('Erro Absoluto vs Erro Relativo');
grid on;

