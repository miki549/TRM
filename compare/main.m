% Főprogram a Modified Next Reaction Method használatához
close all

% Fizikai paraméterek
V_max = 5;                % Maximum sebesség [m/s]
l = 5;                    % Átlagos járműhossz minimális távolsággal [m]
rho_max = 1/l;            % Maximum járműsűrűség [jármű / m]
omega = V_max/rho_max;    % [m^2 / (jármű * s)]
delta_x = 200;            % Cellméret [m]
k_base = omega/delta_x;   % Alapértelmezett átmeneti ráta [m / (jármű * s)]
cm = delta_x/l;           % Cellák maximális kapacitása [jármű]

% Kezdeti állapot (járművek száma szakaszonként)
X = [5;0;40];

Nc = length(X);                   % Útszakaszok száma
cap = cm*ones(Nc, 1);     % Útszakaszok kapacitásai

% Szimuláció paraméterei
tfinal = 400;             % Szimuláció vége

% Reakció mátrix - komplex útkapcsolatokkal
reaction_matrix = [
    1, 2; 
    2, 3;
    3, 1
];

% Időfüggő átmeneti ráta függvény definiálása
% Ez a függvény visszaadja a "from" és "to" útszakaszok közötti 
% átmeneti rátát a "t" időpillanatban
function k = time_dependent_rate(t, from, to)
    % Az alap átmeneti ráta
    k_base = 0.00025;
    %{
    %torlódás a 2-es útszakaszon bizonyos időszakban
    if from == 1 && to == 2
        k = k_base * (1 + 0.5*sin(t/30));  % Ingadozó forgalom
    elseif from == 2 && to == 3
        if t > 100 && t < 200
            k = k_base * 0.5;  % Lassabb áthaladás (forgalmi dugó)
        else
            k = k_base;
        end
    elseif from == 3 && to == 4
        k = k_base * (1 + 0.2*cos(t/50));  % Enyhe ingadozás
    elseif from == 4 && to == 1
        k = k_base * (1 + 0.3*sin(t/40));  % Közepes ingadozás
    else 
    %}
    if from == 1 && to == 2
        k = k_base * (1 + 0.5*sin(t/30));  % Ingadozó forgalom
    else
        k = k_base;  % Alapértelmezett ráta
    end
end

% Modified Next Reaction Method futtatása
[t_history1, X_history1] = MNRM(Nc, X, cap, @time_dependent_rate, reaction_matrix, tfinal);
[t_history2, X_history2] = tRSSA(Nc, X, cap, @time_dependent_rate, reaction_matrix, tfinal);

% Színek generálása a csomópontok számának megfelelően
newcolors = jet(Nc);

% Eredmények ábrázolása külön ablakban két különböző ploton
figure('Position', [100, 100, 900, 600]); % Nagyobb ablak a két subplot számára

% MNRM eredmények ábrázolása felső subplot-ban
subplot(2, 1, 1);
stairs(t_history1, X_history1, 'LineWidth', 2);
colororder(newcolors);
legend_labels = arrayfun(@(x) sprintf('n_%d', x), 1:Nc, 'UniformOutput', false);
legend(legend_labels);
title('Sztochasztikus Forgalom Szimuláció (MNRM Algoritmus)');
ylabel('Járművek száma');
xlabel('Idő [s]');
grid on;
ylim([0, cap(1)]);  % Y-tengely korlátozása a maximális kapacitásra

% tRSSA eredmények ábrázolása alsó subplot-ban
subplot(2, 1, 2);
stairs(t_history2, X_history2, 'LineWidth', 2);
colororder(newcolors);
legend(legend_labels);
title('Sztochasztikus Forgalom Szimuláció (tRSSA Algoritmus)');
ylabel('Járművek száma');
xlabel('Idő [s]');
grid on;
ylim([0, cap(1)]);  % Y-tengely korlátozása a maximális kapacitásra

% Külön ábra az átmeneti ráták időbeli változásának ábrázolásához
figure;
t_range = linspace(0, tfinal, 1000);
k_values = zeros(length(t_range), size(reaction_matrix, 1));

for i = 1:length(t_range)
    for r = 1:size(reaction_matrix, 1)
        from = reaction_matrix(r, 1);
        to = reaction_matrix(r, 2);
        k_values(i, r) = time_dependent_rate(t_range(i), from, to);
    end
end

plot(t_range, k_values, 'LineWidth', 1.5);
title('Időfüggő Átmeneti Ráták');
xlabel('Idő [s]');
ylabel('k_{stoch} érték');
legend_labels = arrayfun(@(x) sprintf('k_{%d→%d}', reaction_matrix(x,1), reaction_matrix(x,2)), 1:size(reaction_matrix, 1), 'UniformOutput', false);
legend(legend_labels);
grid on;

figure;

% Az összes jármű számának kiszámítása minden időpontban
total_cars = zeros(size(t_history1));  % Előre lefoglaljuk a memóriát
for i = 1:length(t_history1)
    row_sum = 0;
    for j = 1:Nc
        row_sum = row_sum + X_history1(i,j);
    end
    total_cars(i) = row_sum;
end

% Ábrázolás
plot(t_history1, total_cars, 'LineWidth', 2, 'Color', 'k');
title('Összes jármű száma az idő függvényében');
xlabel('Idő [s]');
ylabel('Összes jármű száma');
grid on;

% Egyszerű teszteset az analitikus ellenőrzéshez
function test_analytical()
    % Egyszerű eset: konstans propensity függvény
    % Ekkor az integrál egyszerűen: a * τ = S - T
    % Tehát τ = (S - T) / a
    
    % Teszt paraméterek
    S = 1.0;
    T = 0.0;
    t = 0;
    X = [1,2];
    mu = 1;
    
    % Egyszerű konstans átmeneti ráta függvény
    k_stoch_func = @(t, from, to) 0.1;  % Konstans átmeneti ráta
    
    % Reakciómátrix és kapacitás beállítása
    reaction_matrix = [1, 2];
    cap = [20,10];
    Nc = 2;
    tolerance = 1e-8;
    
    % Debug: Propensity értékek kiírása
    from = reaction_matrix(mu, 1);
    to = reaction_matrix(mu, 2);
    k = k_stoch_func(t, from, to);
    actual_propensity = compute_actual_propensity(X, mu, t, k_stoch_func, reaction_matrix, cap, Nc);
    
    fprintf('Debug információk:\n');
    fprintf('k = %f\n', k);
    fprintf('X(from) = %d\n', X(from));
    fprintf('cap(to) = %d\n', cap(to));
    fprintf('X(to) = %d\n', X(to));
    fprintf('Számított propensity = %f\n', actual_propensity);
    
    % Numerikus megoldás
    [tau_numerical, steps] = adaptive_rk4_solver(S, T, t, X, mu, k_stoch_func, reaction_matrix, cap, Nc, tolerance);
    
    % Analitikus megoldás
    % A propensity = k * X(from) * (cap(to) - X(to))
    a = actual_propensity;  % Használjuk a számított propensity értéket
    tau_analytical = (S - T) / a;
    
    % Eredmények kiírása
    fprintf('\nEredmények:\n');
    fprintf('Numerikus megoldás: τ = %f (lépések száma: %d)\n', tau_numerical, steps);
    fprintf('Analitikus megoldás: τ = %f\n', tau_analytical);
    fprintf('Különbség: %e\n', abs(tau_numerical - tau_analytical));
end
    
test_analytical();