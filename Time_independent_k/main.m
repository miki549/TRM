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
Nc = 4; % útszakaszok száma
cap = cm*ones(Nc, 1); % útszakaszok kapacitásai

% Kezdeti állapot (járművek száma szakaszonként)
X = [40; 2; 0; 10];

% Szimuláció paraméterei
tfinal = 400; % Szimuláció vége

% Reakció mátrix - komplex útkapcsolatokkal
reaction_matrix = [
    1, 2, -1, 1; 
    2, 3, -1, 1; 
    3, 4, -1, 1;
    4, 1, -1, 1;
    3, 2, -1, 1
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

stairs(t_history, X_history, 'LineWidth', 1.5);
title('Stochasztikus Szimuláció (Gillespie Algoritmus)');
ylabel('Járművek száma');
xlabel('Idő');
grid on;
legend(legend_labels);
colororder(newcolors);
plot(2,1,1);
