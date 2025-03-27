close all

% Fizikai paraméterek
V_max = 5; % maximum sebesség [m/s]
l = 5; % átlagos járműhossz minimális távolsággal [m]
rho_max = 1/l; % maximum járműsűrűség [jármű / m]
omega = V_max/rho_max; % [m^2 / (jármű * s)]
delta_x = 200; % cellméret [m]
k = omega/delta_x; % átmeneti ráta cellák között [m / (jármű * s)]
cm = delta_x/l; % cellák maximális kapacitása [jármű]

% Reakció (átmeneti) ráta
k_stoch = k / delta_x;

% 10 csomópontos úthálózat
Nc = 5; % útszakaszok száma
cap = cm*ones(Nc, 1); % útszakaszok kapacitásai

% Kezdeti állapot (járművek száma szakaszonként)
X = [40; 2; 0; 10; 5];

% Szimuláció paraméterei
tfinal = 400; % Szimuláció vége

% Reakció mátrix - komplex útkapcsolatokkal
reaction_matrix = [
    1, 2, -1, 1;  % 1-es csomópontból 2-es csomópontba 
    1, 3, -1, 1;  % 1-es csomópontból 3-as csomópontba
    2, 4, -1, 1;  % 2-es csomópontból 4-es csomópontba
    3, 4, -1, 1   % 3-as csomópontból 4-es csomópontba
    1, 4, -1, 1   % 1-es csomópontból 4-es csomópontba (átlós él)
];

% Gillespie szimuláció futtatása
[t_history, X_history] = gillespie_simulation(Nc, X, cap, k_stoch, reaction_matrix, tfinal);

% Plotolás
figure;
stairs(t_history, X_history, 'LineWidth', 2);

% Színek generálása a csomópontok számának megfelelően
newcolors = jet(Nc);
       
colororder(newcolors);

% Jelmagyarázat dinamikus generálása
legend_labels = arrayfun(@(x) sprintf('n_%d', x), 1:Nc, 'UniformOutput', false);
legend(legend_labels);

title('TRM Stochasztikus Szimuláció - Gillespie Algoritmus (10 csomópont)');
ylabel('Járművek száma');
xlabel('Idő');
grid on;

% Kapcsolati mátrix létrehozása
K = zeros(Nc, Nc);
kij = k;
for i = 1:size(reaction_matrix, 1)
    from = reaction_matrix(i, 1);
    to = reaction_matrix(i, 2);
    K(from, to) = kij;
end
K = K/delta_x;

% Determinisztikus szimuláció
y0 = X;
tspan = [0 tfinal];
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
[t_det, y_det] = ode45(@(t,y) TRM_red_ode_mod(t,y,K,cap), tspan, y0, options);

% Két szimuláció összehasonlítása
figure;
subplot(2,1,1);
stairs(t_history, X_history, 'LineWidth', 1.5);
title('Stochasztikus Szimuláció (Gillespie Algoritmus)');
ylabel('Járművek száma');
xlabel('Idő');
grid on;
legend(legend_labels);
colororder(newcolors);

subplot(2,1,2);
plot(t_det, y_det, 'LineWidth', 1.5);
title('Determinisztikus Szimuláció (ODE)');
ylabel('Járművek száma');
xlabel('Idő');
grid on;
legend(legend_labels);
colororder(newcolors);