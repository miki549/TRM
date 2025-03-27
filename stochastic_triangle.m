clear all
close all

%physical parameters of continuous time model where the state variable is vehicle density [vehicle / m]:
V_max = 5; %maximum speed [m/s]
l = 5; %average vehicle length including minimum space between vehicles [m]
rho_max = 1/l; %maximum vehicle density [vehicle / m]
omega = V_max/rho_max; %[m^2 / (vehicle * s)]
delta_x=200; %cell size [m]
k = omega/delta_x; %transition rate between cells [m / (vehicle * s)]
cm = delta_x/l; %maximum capacity of cells [vehicle]

%reaction (transition) rate of models where the state variable is the number of vehicles (spaces) in each cell:
%this is simply obtained by transforming the state variables in the kinetic ODEs from densities to vehicle numbers
%needed to compute directly in vehicle / space units (otherwise not important)
k_stoch = k / delta_x;

% Háromszög úthálózat: 3 útszakasz
Nc = 3; %number of compartments
cap = 40*ones(Nc, 1); %capacities of compartments

% Kezdeti állapot (járművek száma szakaszonként)
X = [30; 5; 15];

% Gillespie algoritmus paraméterei
tfinal = 400; % Szimuláció vége
t = 0; % Szimuláció kezdete
t_history = [0]; % Időpontok tárolása
X_history = X'; % Állapotok tárolása

% Reakció mátrix - minden lehetséges reakciót definiál
% Formátum: [honnan, hová, változás a honnan-ban, változás a hová-ban]
reaction_matrix = [
    1, 2, -1, 1;   % 1-es útszakaszról 2-esre
    2, 3, -1, 1;   % 2-es útszakaszról 3-asra 
    3, 1, -1, 1;   % 3-as útszakaszról 1-esre
];

j = 0; % Reakció számláló

% Gillespie algoritmus
while t < tfinal
    j = j + 1;
    
    % 1. lépés: Propensity/intenzitás függvények kiszámítása minden reakcióra
    propensities = zeros(size(reaction_matrix, 1), 1);
    
    for r = 1:size(reaction_matrix, 1)
        from = reaction_matrix(r, 1);
        to = reaction_matrix(r, 2);
        
        % Propensity: k_stoch * n_i * (c_j - n_j)
        propensities(r) = k_stoch * X(from) * (cap(to) - X(to));
    end
    %Külső flow hozzáadása
    flows = TRM_external_flows(t,Nc);
    inflow_rates = flows(:,1);
    outflow_rates = flows(:,2);
    
    % Külső flow propensity-jei
    inflow_propensities = inflow_rates .* (cap - X);
    outflow_propensities = outflow_rates .* X;
    
    % Összes reakció propensity-jének összegyűjtése
    all_propensities = [propensities; inflow_propensities; outflow_propensities];
    total_propensity = sum(all_propensities);
    
    % Ha nincs lehetséges reakció, kilép
    if total_propensity == 0
        break;
    end
    
    % 2-3. lépés: Két véletlenszám generálása
    r1 = rand();
    r2 = rand();
    
    % 4. lépés: Várakozási idő kiszámítása
    delta_t = log(1/r1) / total_propensity;
    
    % Idő frissítése
    t = t + delta_t;
    if t > tfinal
        break;
    end
    
    % 5. lépés: Melyik reakció következik be?
    cum_prop = cumsum(all_propensities) / total_propensity;
    mu = find(cum_prop >= r2, 1);
    
    % Állapot frissítése a választott reakció alapján
    if mu <= size(reaction_matrix, 1)
        % Belső átmenet két útszakasz között
        from = reaction_matrix(mu, 1);
        to = reaction_matrix(mu, 2);
        X(from) = X(from) - 1;
        X(to) = X(to) + 1;
    elseif mu <= size(reaction_matrix, 1) + Nc
        % Külső inflow
        inflow_idx = mu - size(reaction_matrix, 1);
        X(inflow_idx) = X(inflow_idx) + 1;
    else
        % Külső outflow
        outflow_idx = mu - size(reaction_matrix, 1) - Nc;
        X(outflow_idx) = X(outflow_idx) - 1;
    end
    
    % Adatok tárolása
    t_history = [t_history, t];
    X_history = [X_history; X'];
end

% Plotolás
figure;
stairs(t_history, X_history, 'LineWidth', 2);

newcolors = [0.83 0.14 0.14; ...
             1.00 0.54 0.00; ...
             0.2 0.1 1];
       
colororder(newcolors);

legend('n_1', 'n_2', 'n_3');
title('TRM Stochastic Simulation - Gillespie Algorithm');
ylabel('Number of vehicles');
xlabel('Time');
grid on;

% Determinisztikus szimuláció futtatása összehasonlításhoz
% Kapcsolati mátrix létrehozása
K = zeros(Nc, Nc);
kij = k;
K(1,2) = kij; % 1-es útszakaszról 2-esre
K(2,3) = kij; % 2-es útszakaszról 3-asra
K(3,1) = kij; % 3-as útszakaszról 1-esre
K = K/delta_x;

% Determinisztikus szimuláció
y0 = [30; 5; 15]; % Ugyanaz a kezdeti állapot
tspan = [0 tfinal];
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
[t_det, y_det] = ode45(@(t,y) TRM_red_ode_mod(t,y,K,cap), tspan, y0, options);

% Két szimuláció összehasonlítása
figure;
subplot(2,1,1);
stairs(t_history, X_history, 'LineWidth', 1.5);
title('Stochastic Simulation (Gillespie Algorithm)');
ylabel('Number of vehicles');
xlabel('Time');
grid on;
legend('n_1', 'n_2', 'n_3');
colororder(newcolors);

subplot(2,1,2);
plot(t_det, y_det, 'LineWidth', 1.5);
title('Deterministic Simulation (ODE)');
ylabel('Number of vehicles');
xlabel('Time');
grid on;
legend('n_1', 'n_2', 'n_3');
colororder(newcolors);