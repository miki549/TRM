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
tfinal = 100; % Szimuláció vége

% Reakció mátrix - komplex útkapcsolatokkal
reaction_matrix = [
    1, 2, -1, 1; 
    2, 3, -1, 1; 
    3, 1, -1, 1;
];


% Plotolás

% Színek generálása a csomópontok számának megfelelően
newcolors = jet(Nc);     
colororder(newcolors);

% Jelmagyarázat dinamikus generálása
legend_labels = arrayfun(@(x) sprintf('n_%d', x), 1:Nc, 'UniformOutput', false);
legend(legend_labels);

% Időfüggő reakcióráta definiálása (pl. szinuszos)
k_stoch_time_f = @(t) 0.1+ 0.05 * sin(2*pi*t/10);
[t_history, X_history] = time_dependent_gillespie(Nc, X, cap, k_stoch_time_f, reaction_matrix, tfinal, 100,20);
stairs(t_history, X_history, 'LineWidth', 1.5);
title('Időfüggő Stochasztikus Szimuláció)');
ylabel('Járművek száma');
xlabel('Idő');
grid on;
legend(legend_labels);
colororder(newcolors);
plot(2,1,1);
