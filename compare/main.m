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

%% Időfüggő átmeneti ráta függvény definiálása
function k = time_dependent_rate(t, from, to)
    % Alapértelmezett átmeneti ráta
    k_base = 0.00025;
    
    if from == 1 && to == 2
        % Szinuszos ingadozás az 1->2 szakaszon
        k = k_base * (1 + 0.5*sin(t/30));
    elseif from == 3 && to == 4
        % Enyhe koszinuszos ingadozás a 3->4 szakaszon
        k = k_base * (1 + 0.2*cos(t/50));
    elseif from == 2 && to == 3
        % Tanh alapú átmenet a 2->3 szakaszon
        center = 150;  % Átmenet középpontja
        width = 20;    % Átmenet szélessége
        k = k_base * (0.75 + 0.25*tanh((t-center)/width));
    elseif from == 4 && to == 1
        % Erősebb periodikus változás a 4->1 szakaszon
        k = k_base * (1 + 0.7*sin(t/40 + pi/4));
    else 
        % Alapértelmezett konstans ráta
        k = k_base;
    end
end

%% Determinisztikus modell előkészítése
% Átmeneti ráta mátrix létrehozása
K = zeros(Nc, Nc);
for i = 1:size(reaction_matrix, 1)
    from = reaction_matrix(i, 1);
    to = reaction_matrix(i, 2);
    K(from, to) = k_base;
end
K = K / delta_x;

%% Szimulációk futtatása
% Paraméterek
num_runs = 200;  % Futtatások száma
common_time_points = linspace(0, tfinal, 1000);  % Közös időpontok az összehasonlításhoz

% Eredmények tárolása
all_MNRM_results = zeros(num_runs, length(common_time_points), Nc);
all_tRSSA_results = zeros(num_runs, length(common_time_points), Nc);

fprintf('Szimulációk futtatása (%d futtatás):\n', num_runs);

% Szimulációk végrehajtása
for i = 1:num_runs
    fprintf('Futtatás %d/%d...\n', i, num_runs);
    
    % Állapotvektor inicializálása minden futtatáshoz
    X = X_init;
    
    % MNRM algoritmus futtatása
    [t_history_MNRM, X_history_MNRM] = MNRM(Nc, X, cap, @time_dependent_rate, reaction_matrix, tfinal);
    
    % tRSSA algoritmus futtatása
    [t_history_tRSSA, X_history_tRSSA] = tRSSA(Nc, X, cap, @time_dependent_rate, reaction_matrix, tfinal);
    
    % Eredmények interpolálása a közös időpontokra
    for j = 1:Nc
        all_MNRM_results(i,:,j) = interp1(t_history_MNRM, X_history_MNRM(:,j), common_time_points, 'previous');
        all_tRSSA_results(i,:,j) = interp1(t_history_tRSSA, X_history_tRSSA(:,j), common_time_points, 'previous');
    end
end

%% Átlagok számítása
% MNRM átlagok számítása
mean_MNRM = zeros(length(common_time_points), Nc);
for j = 1:Nc
    mean_MNRM(:,j) = mean(all_MNRM_results(:,:,j));
end

% tRSSA átlagok számítása
mean_tRSSA = zeros(length(common_time_points), Nc);
for j = 1:Nc
    mean_tRSSA(:,j) = mean(all_tRSSA_results(:,:,j));
end

%% Determinisztikus modell futtatása
% ODE megoldó beállításai és futtatása
y0 = X_init;
tspan = [0 tfinal];
options = odeset('RelTol', 1e-8, 'AbsTol', 1e-8);
[t_det, y_det] = ode45(@(t,y) TRM_red_ode_mod(t, y, K, cap, @time_dependent_rate), tspan, y0, options);

% Determinisztikus eredmények interpolálása
det_interp = zeros(length(common_time_points), Nc);
for j = 1:Nc
    det_interp(:,j) = interp1(t_det, y_det(:,j), common_time_points);
end

%% Egyedi futtatások vizualizációhoz
X = X_init;
[t_history1, X_history1] = MNRM(Nc, X, cap, @time_dependent_rate, reaction_matrix, tfinal);
[t_history2, X_history2] = tRSSA(Nc, X, cap, @time_dependent_rate, reaction_matrix, tfinal);

%% Eredmények megjelenítése
% Színek beállítása
newcolors = jet(Nc);

% MNRM és tRSSA összehasonlító ábrák
figure('Position', [100, 100, 900, 600])

% MNRM eredmények
subplot(2, 1, 1)
stairs(t_history1, X_history1, 'LineWidth', 2)
colororder(newcolors)
legend_labels = arrayfun(@(x) sprintf('n_%d', x), 1:Nc, 'UniformOutput', false);
legend(legend_labels)
title('Sztochasztikus Forgalom Szimuláció (MNRM Algoritmus)')
ylabel('Járművek száma')
xlabel('Idő [s]')
grid on
ylim([0, cap(1)])

% tRSSA eredmények
subplot(2, 1, 2)
stairs(t_history2, X_history2, 'LineWidth', 2)
colororder(newcolors)
legend(legend_labels)
title('Sztochasztikus Forgalom Szimuláció (tRSSA Algoritmus)')
ylabel('Járművek száma')
xlabel('Idő [s]')
grid on
ylim([0, cap(1)])

%% Összes jármű számának ábrázolása
figure
total_cars = sum(X_history1, 2);
plot(t_history1, total_cars, 'LineWidth', 2, 'Color', 'k')
title('Összes jármű száma az idő függvényében')
xlabel('Idő [s]')
ylabel('Összes jármű száma')
grid on

%% Propensity függvények vizualizációja
plot_propensities(t_history1, X_history1, reaction_matrix, cap, @time_dependent_rate, Nc, 0, tfinal, 1000);

%% Átlagolt eredmények és determinisztikus modell összehasonlítása
figure('Position', [100, 100, 1000, 600])

% Vonaltípusok és markerek beállítása
line_styles = {'--', '-', '-'};
markers = {'none', 'o', 'none'};
marker_size = [0, 6, 0];
marker_interval = [1, 40, 1];

hold on
% Első útszakasz mindhárom módszerrel
h1 = plot(common_time_points, mean_MNRM(:,1), line_styles{1}, 'LineWidth', 2, 'Color', newcolors(1,:), 'Marker', markers{1});
h2 = plot(common_time_points, mean_tRSSA(:,1), line_styles{2}, 'LineWidth', 2, 'Color', newcolors(1,:), ...
    'Marker', markers{2}, 'MarkerSize', marker_size(2), ...
    'MarkerIndices', 1:marker_interval(2):length(common_time_points), ...
    'MarkerFaceColor', newcolors(1,:));
h3 = plot(common_time_points, det_interp(:,1), line_styles{3}, 'LineWidth', 2.5, 'Color', newcolors(1,:), 'Marker', markers{3});

% Többi útszakasz ábrázolása
for j = 2:Nc
    plot(common_time_points, mean_MNRM(:,j), line_styles{1}, 'LineWidth', 2, 'Color', newcolors(j,:), 'Marker', markers{1});
    plot(common_time_points, mean_tRSSA(:,j), line_styles{2}, 'LineWidth', 2, 'Color', newcolors(j,:), ...
        'Marker', markers{2}, 'MarkerSize', marker_size(2), ...
        'MarkerIndices', 1:marker_interval(2):length(common_time_points), ...
        'MarkerFaceColor', newcolors(j,:));
    plot(common_time_points, det_interp(:,j), line_styles{3}, 'LineWidth', 2.5, 'Color', newcolors(j,:), 'Marker', markers{3});
end
hold off

title('Algoritmusok összehasonlítása (200 futtatás átlaga)')
ylabel('Járművek száma')
xlabel('Idő [s]')
grid on
ylim([0, cap(1)])
legend([h1, h2, h3], 'MNRM (átlag)', 'tRSSA (átlag)', 'Determinisztikus', 'Location', 'northwest')

%% Determinisztikus modell külön ábrázolása
figure('Position', [100, 100, 800, 500])
plot(t_det, y_det, 'LineWidth', 2.5)
colororder(newcolors)
title('Determinisztikus Forgalom Szimuláció')
ylabel('Járművek száma')
xlabel('Idő [s]')
grid on
ylim([0, cap(1)])
legend(legend_labels, 'Location', 'best')