clear all
close all

% Folytonos idejű modell fizikai paraméterei (járműsűrűség [jármű/m]):
V_max = 5;           % Maximális sebesség [m/s]
l = 5;               % Átlagos járműhossz a minimális követési távolsággal [m]
rho_max = 1/l;       % Maximális járműsűrűség [jármű/m]
omega = V_max/rho_max; % [m^2/(jármű*s)]
delta_x = 200;       % Cella méret [m]
k = omega/delta_x;   % Átmeneti ráta cellák között [m/(jármű*s)]
cm = delta_x/l;      % Cellák maximális kapacitása [jármű]

% Reakció (átmeneti) ráta a modellben, ahol az állapotváltozó a járművek száma
k_stoch = k / delta_x;

Nc = 3;              % Cellák száma
cc = (delta_x / l) * ones(Nc,1); % Cellák kapacitása

% Átmeneti mátrix létrehozása
K = zeros(Nc,Nc);
kij = k;
K(1,2) = kij;
K(2,3) = kij;
K(3,1) = kij;

K = K/delta_x;       % Állapotváltozók járműszámra konvertálása

cap = 40*ones(Nc,1); % Cellák kapacitása

% Kezdeti állapot
n0 = [30; 5; 15];
y0 = [n0];

% Szimuláció időtartama
tfinal = 200;
tspan = [0 tfinal];

% ODE megoldó beállításai
options = odeset('RelTol',1e-8,'AbsTol',1e-8);

% Differenciálegyenlet megoldása
[t_det, y_det] = ode45(@(t,y) TRM_red_ode_mod(t,y,K,cap), tspan, y0, options);

% Eredmények ábrázolása
figure
plot(t_det, y_det, 'LineWidth', 2)

% Színek beállítása
newcolors = [0.83 0.14 0.14; ...
             1.00 0.54 0.00; ...
             0.2 0.1 1];
       
colororder(newcolors)

legend('n_1','n_2','n_3')
title('Determinisztikus háromszög hálózat')
ylabel('Állapotok')
xlabel('Idő')
colororder(newcolors)
grid





