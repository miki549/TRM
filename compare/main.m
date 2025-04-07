clear all
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

% Reakció mátrix
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
    %torlódás a 2-es útszakaszon bizonyos időszakban
    if from == 1 && to == 2
        k = k_base * (1 + 0.5*sin(t/30));  % Ingadozó forgalom
    elseif from == 3 && to == 1
        k = k_base * (1 + 0.2*cos(t/50));  % Enyhe ingadozás
    elseif from == 2 && to == 3
        if t>100 && t<200
            k = k_base * 0.5;
        else
            k = k_base;
        end
    else 
        k = k_base;  % Alapértelmezett ráta
    end
end

% Modified Next Reaction Method futtatása
[t_history1, X_history1, integral_check] = MNRM(Nc, X, cap, @time_dependent_rate, reaction_matrix, tfinal);
%[t_history2, X_history2] = tRSSA(Nc, X, cap, @time_dependent_rate, reaction_matrix, tfinal);

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
%{
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
%}
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

% Ábra készítése az integrálok összehasonlításához
figure('Position', [100, 100, 1200, 800]);

% Subplot 1: Időbeli lefutás
subplot(2,1,1);
hold on;
plot(integral_check.time, integral_check.S_minus_T, 'k-', 'LineWidth', 2, 'DisplayName', 'S-T');
plot(integral_check.time, integral_check.trap_integral, 'b--', 'LineWidth', 1.5, 'DisplayName', 'Trapéz integrál');
xlabel('Idő');
ylabel('Integrál érték');
title('Integrálok időbeli összehasonlítása');
legend('Location', 'best');
grid on;
hold off;