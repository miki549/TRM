% Főprogram a time-dependent Rejection-based SSA használatához
close all

% Fizikai paraméterek
V_max = 5;                % Maximum sebesség [m/s]
l = 5;                    % Átlagos járműhossz minimális távolsággal [m]
rho_max = 1/l;            % Maximum járműsűrűség [jármű / m]
omega = V_max/rho_max;    % [m^2 / (jármű * s)]
delta_x = 200;            % Cellméret [m]
k_base = omega/delta_x;   % Alapértelmezett átmeneti ráta [m / (jármű * s)]
cm = delta_x/l;           % Cellák maximális kapacitása [jármű]

% 4 csomópontos úthálózat
Nc = 4;                   % Útszakaszok száma
cap = cm*ones(Nc, 1);     % Útszakaszok kapacitásai

% Kezdeti állapot (járművek száma szakaszonként)
X = [40; 2; 0; 0];

% Szimuláció paraméterei
tfinal = 400;             % Szimuláció vége

% Reakció mátrix - komplex útkapcsolatokkal
reaction_matrix = [
    1, 2; 
    2, 3; 
    3, 4;
    4, 1;
];

% Időfüggő átmeneti ráta függvény definiálása
% Ez a függvény visszaadja a "from" és "to" útszakaszok közötti 
% átmeneti rátát a "t" időpillanatban
function k = time_dependent_rate(t, from, to)
    % Az alap átmeneti ráta
    k_base = 0.00025;
    
    % Alap forgalmi mintázatok időbeli változása
    % Például torlódás a 2-es útszakaszon bizonyos időszakban
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
        k = k_base;  % Alapértelmezett ráta
    end
end

% tRSSA algoritmus futtatása
[t_history, X_history] = time_dependent_rejection_SSA(Nc, X, cap, @time_dependent_rate, reaction_matrix, tfinal);

% Eredmények ábrázolása
figure;
stairs(t_history, X_history, 'LineWidth', 2);

% Színek generálása a csomópontok számának megfelelően
newcolors = jet(Nc);
colororder(newcolors);

% Jelmagyarázat dinamikus generálása
legend_labels = arrayfun(@(x) sprintf('n_%d', x), 1:Nc, 'UniformOutput', false);
legend(legend_labels);

title('Sztochasztikus Forgalom Szimuláció (tRSSA Algoritmus)');
ylabel('Járművek száma');
xlabel('Idő [s]');
grid on;
ylim([0, cap(1)]);  % Y-tengely korlátozása a maximális kapacitásra

% Átmeneti ráták időbeli változásának ábrázolása
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