clear
n = 100;
H = gerarMatrizHilbert(n);
b=ones(n,1);

x0 = zeros(n,1);
omega1 = 0.02 ;
omega2 = 1.25;
tolerance = 1e-7;
max_iterations = 10;
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


graphics_toolkit('gnuplot');
omega_1 = 0.01:0.01:interv;
radios_espectrais_JOR = zeros(size(omega_1));

for i = 1:length(omega_1)
    Bjo = Bjor(H, omega_1(i));
    radios_espectrais_JOR(i) = spec_rad(Bjo);
end

figure; 
plot(omega_1, radios_espectrais_JOR, 'LineWidth', 2);
xlabel('Omega');
ylabel('Raio Espectral');
legend('JOR');
grid on;

omega_2 = 0.01:0.01:1.99;
radios_espectrais_SOR = zeros(size(omega_2));

for i = 1:length(omega_2)
    Bsor_func = Bsor(H, omega_2(i));
    radios_espectrais_SOR(i) = spec_rad(Bsor_func);
end

figure;
plot(omega_2, radios_espectrais_SOR, 'LineWidth', 2);
xlabel('Omega');
ylabel('Raio Espectral');
legend('SOR');
grid on;


