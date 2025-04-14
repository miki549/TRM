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
X_init = [5;1;30;13;40];

Nc = length(X_init);      % Útszakaszok száma
cap = cm*ones(Nc, 1);     % Útszakaszok kapacitásai

% Szimuláció paraméterei
tfinal = 400;             % Szimuláció vége

% Reakció mátrix
reaction_matrix = [
    1, 2; 
    2, 3;
    3, 4;
    4, 1;
    4, 5;
    5, 3;
];

% Időfüggő átmeneti ráta függvény definiálása
% Ez a függvény visszaadja a "from" és "to" útszakaszok közötti 
% átmeneti rátát a "t" időpillanatban
function k = time_dependent_rate(t, from, to)
    % Az alap átmeneti ráta
    k_base = 0.00025;
    
    if from == 1 && to == 2
        k = k_base * (1 + 0.5*sin(t/30));  % Ingadozó forgalom
    elseif from == 3 && to == 4
        k = k_base * (1 + 0.2*cos(t/50));  % Enyhe ingadozás
    elseif from == 2 && to == 3
        % Tanh alapú átmenet a 2->3 szakaszon
        center = 150; % Az átmenet középpontja
        width = 20;   % Az átmenet szélessége
        k = k_base * (0.75 + 0.25*tanh((t-center)/width));
    elseif from == 4 && to == 1
        % Periodikus változás nagyobb amplitúdóval a 4->1 szakaszon
        k = k_base * (1 + 0.7*sin(t/40 + pi/4));  % Fáziseltolt, erősebb ingadozás
    else 
        k = k_base;  % Alapértelmezett ráta
    end
end

% Átmeneti ráta mátrix a determinisztikus modellhez
K = zeros(Nc, Nc);
for i = 1:size(reaction_matrix, 1)
    from = reaction_matrix(i, 1);
    to = reaction_matrix(i, 2);
    K(from, to) = k_base;
end
K = K / delta_x;

% 50 futtatás mindkét algoritmussal
num_runs = 200;
common_time_points = linspace(0, tfinal, 1000);

% Inicializáljuk a tároló mátrixokat a közös időpontokra interpolált értékeknek
all_MNRM_results = zeros(num_runs, length(common_time_points), Nc);
all_tRSSA_results = zeros(num_runs, length(common_time_points), Nc);

fprintf('Szimulációk futtatása (%d futtatás):\n', num_runs);

for i = 1:num_runs
    fprintf('Futtatás %d/%d...\n', i, num_runs);
    
    % Inicializáljuk az X vektort minden futtatás előtt
    X = X_init;
    
    % MNRM futtatása
    [t_history_MNRM, X_history_MNRM] = MNRM(Nc, X, cap, @time_dependent_rate, reaction_matrix, tfinal);
    
    % tRSSA futtatása
    [t_history_tRSSA, X_history_tRSSA] = tRSSA(Nc, X, cap, @time_dependent_rate, reaction_matrix, tfinal);
    
    % Adatok interpolálása a közös időpontokra
    for j = 1:Nc
        all_MNRM_results(i,:,j) = interp1(t_history_MNRM, X_history_MNRM(:,j), common_time_points, 'previous');
        all_tRSSA_results(i,:,j) = interp1(t_history_tRSSA, X_history_tRSSA(:,j), common_time_points, 'previous');
    end
end

% Átlagok számítása
mean_MNRM = squeeze(mean(all_MNRM_results, 1));
mean_tRSSA = squeeze(mean(all_tRSSA_results, 1));

% Determinisztikus modell futtatása időfüggő rátákkal
y0 = X_init;
tspan = [0 tfinal];
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);

% Átadjuk a time_dependent_rate függvényt a TRM_red_ode_mod-nak
[t_det, y_det] = ode45(@(t,y) TRM_red_ode_mod(t, y, K, cap, @time_dependent_rate), tspan, y0, options);

% Interpoláljuk a determinisztikus eredményeket is a közös időpontokra
det_interp = zeros(length(common_time_points), Nc);
for j = 1:Nc
    det_interp(:,j) = interp1(t_det, y_det(:,j), common_time_points);
end

% Egyetlen futtatás az MNRM és tRSSA algoritmusokkal (az eredeti plotokhoz)
X = X_init;
[t_history1, X_history1] = MNRM(Nc, X, cap, @time_dependent_rate, reaction_matrix, tfinal);
[t_history2, X_history2] = tRSSA(Nc, X, cap, @time_dependent_rate, reaction_matrix, tfinal);

% Színek generálása a csomópontok számának megfelelően
newcolors = jet(Nc);

% EREDETI PLOTOK MEGJELENÍTÉSE 
% ============================

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

% Összes jármű számának kiszámítása minden időpontban
figure;
total_cars = zeros(size(t_history1));
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

% Propensity függvények vizualizálása
t_start = 0;
t_end = tfinal;
num_points = 1000;  % Ennyi pontban számoljuk ki az értékeket
plot_propensities(t_history1, X_history1, reaction_matrix, cap, @time_dependent_rate, Nc, t_start, t_end, num_points);

% ÚJ ÖSSZEHASONLÍTÓ PLOT AZ ÁTLAGOKKAL ÉS DETERMINISZTIKUS MODELLEL
% =================================================================

% Új ábra az átlagolt eredmények és determinisztikus modell összehasonlításához
figure('Position', [100, 100, 1000, 600]);

% Vonaltípusok és markerek meghatározása a kért módon
line_styles = {'--', '-', '-'};
markers = {'none', 'o', 'none'};
marker_size = [0, 6, 0];
marker_interval = [1, 40, 1];  % tRSSA minden 40. pontban marker

hold on;
% Első útszakasz mindhárom módszerrel - ezeket használjuk a legend-hez
h1 = plot(common_time_points, mean_MNRM(:,1), line_styles{1}, 'LineWidth', 2, 'Color', newcolors(1,:), 'Marker', markers{1});
h2 = plot(common_time_points, mean_tRSSA(:,1), line_styles{2}, 'LineWidth', 2, 'Color', newcolors(1,:), 'Marker', markers{2}, 'MarkerSize', marker_size(2), 'MarkerIndices', 1:marker_interval(2):length(common_time_points), 'MarkerFaceColor', newcolors(1,:));
h3 = plot(common_time_points, det_interp(:,1), line_styles{3}, 'LineWidth', 2.5, 'Color', newcolors(1,:), 'Marker', markers{3});

% Többi útszakasz hozzáadása
for j = 2:Nc
    plot(common_time_points, mean_MNRM(:,j), line_styles{1}, 'LineWidth', 2, 'Color', newcolors(j,:), 'Marker', markers{1});
    plot(common_time_points, mean_tRSSA(:,j), line_styles{2}, 'LineWidth', 2, 'Color', newcolors(j,:), 'Marker', markers{2}, 'MarkerSize', marker_size(2), 'MarkerIndices', 1:marker_interval(2):length(common_time_points), 'MarkerFaceColor', newcolors(j,:));
    plot(common_time_points, det_interp(:,j), line_styles{3}, 'LineWidth', 2.5, 'Color', newcolors(j,:), 'Marker', markers{3});
end
hold off;

title('Algoritmusok összehasonlítása (50 futtatás átlaga)');
ylabel('Járművek száma');
xlabel('Idő [s]');
grid on;
ylim([0, cap(1)]);

leg1 = legend([h1, h2, h3], 'MNRM (átlag)', 'tRSSA (átlag)', 'Determinisztikus', 'Location', 'northwest');

% Determinisztikus modell eredményeinek külön ábrázolása
figure('Position', [100, 100, 800, 500]);

% Determinisztikus eredmények ábrázolása
plot(t_det, y_det, 'LineWidth', 2.5);
colororder(newcolors);
title('Determinisztikus Forgalom Szimuláció');
ylabel('Járművek száma');
xlabel('Idő [s]');
grid on;
ylim([0, cap(1)]);

legend(legend_labels, 'Location', 'best');