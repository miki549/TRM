% Test script a trapezoid_tau függvény ellenőrzésére
clear all;
close all;

% Teszt paraméterek
Nc = 3;
cap = [100; 100; 100];
X = [50; 30; 20];
tolerance = 1e-4;

% Reakció mátrix (1->2, 2->3, 3->1)
reaction_matrix = [1, 2; 2, 3; 3, 1];

% 1. Teszt: Konstans propensity függvény
% Ebben az esetben az analitikus megoldás: tau = (S-T)/a
fprintf('1. Teszt: Konstans propensity függvény\n');

% Konstans rate függvény
k_const = @(t, from, to) 0.001;

S = 1.5;  % Példa S érték
T = 0;    % Kezdeti T érték
t = 0;    % Kezdeti időpont

% Analitikus megoldás kiszámítása
initial_prop = compute_actual_propensity(X, 1, t, k_const, reaction_matrix, cap, Nc);
tau_analytic = (S-T)/initial_prop;

% Numerikus megoldás
[tau_numeric, integral] = trapezoid_tau(S, T, t, X, 1, k_const, reaction_matrix, cap, Nc, tolerance);

fprintf('Analitikus tau: %.6f\n', tau_analytic);
fprintf('Numerikus tau: %.6f\n', tau_numeric);
fprintf('Relatív hiba: %.2e\n\n', abs(tau_analytic-tau_numeric)/tau_analytic);

% 2. Teszt: Lineáris propensity függvény
% a(t) = a0*(1 + kt)
fprintf('2. Teszt: Lineáris propensity függvény\n');

k_linear = @(t, from, to) 0.001*(1 + 0.1*t);

% Numerikus megoldás a lineáris esetben
[tau_numeric_lin, integral_lin] = trapezoid_tau(S, T, t, X, 1, k_linear, reaction_matrix, cap, Nc, tolerance);

% Ellenőrzés: az integrál értékének egyeznie kell S-T-vel
error_lin = abs(integral_lin - (S-T));
fprintf('Integrál érték: %.6f\n', integral_lin);
fprintf('Elvárt érték (S-T): %.6f\n', S-T);
fprintf('Abszolút hiba: %.2e\n\n', error_lin);

% 3. Teszt: Szinuszos propensity függvény
fprintf('3. Teszt: Szinuszos propensity függvény\n');

k_sin = @(t, from, to) 0.001*(1 + 0.5*sin(t));

% Különböző tolerancia értékek tesztelése
tolerances = [1e-2, 1e-4, 1e-6];
for tol = tolerances
    [tau_sin, integral_sin] = trapezoid_tau(S, T, t, X, 1, k_sin, reaction_matrix, cap, Nc, tol);
    fprintf('Tolerancia: %.0e\n', tol);
    fprintf('Tau érték: %.6f\n', tau_sin);
    fprintf('Integrál érték: %.6f\n', integral_sin);
    fprintf('Integrál hiba: %.2e\n\n', abs(integral_sin - (S-T)));
end

% 4. Teszt: Konvergencia vizsgálat
fprintf('4. Teszt: Konvergencia vizsgálat\n');

iterations = 100;
tau_values = zeros(iterations, 1);
integral_values = zeros(iterations, 1);

for i = 1:iterations
    tol = 10^(-i/10);  % Fokozatosan csökkenő tolerancia
    [tau_values(i), integral_values(i)] = trapezoid_tau(S, T, t, X, 1, k_sin, reaction_matrix, cap, Nc, tol);
end

% Konvergencia plot
figure;
semilogx(10.^(-((1:iterations)/10)), tau_values, 'b.-');
xlabel('Tolerancia');
ylabel('Tau érték');
title('Tau konvergencia vizsgálat');
grid on;

% Integrál érték hiba plot
figure;
loglog(10.^(-((1:iterations)/10)), abs(integral_values - (S-T)), 'r.-');
xlabel('Tolerancia');
ylabel('Integrál hiba');
title('Integrál hiba vs. Tolerancia');
grid on;

% 5. Teszt: Propensity függvény vizualizáció
fprintf('5. Teszt: Propensity függvény vizualizáció\n');

t_range = linspace(0, max(tau_values)*1.5, 1000);
props = zeros(size(t_range));

for i = 1:length(t_range)
    props(i) = compute_actual_propensity(X, 1, t_range(i), k_sin, reaction_matrix, cap, Nc);
end

figure;
plot(t_range, props, 'b-', 'LineWidth', 1.5);
hold on;
plot([0, tau_numeric], [props(1), compute_actual_propensity(X, 1, tau_numeric, k_sin, reaction_matrix, cap, Nc)], 'r--');
xlabel('Idő');
ylabel('Propensity érték');
title('Propensity függvény és a számított tau');
legend('Propensity függvény', 'Számított tau intervallum');
grid on; 