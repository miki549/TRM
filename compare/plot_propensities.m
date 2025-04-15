function plot_propensities(t_history, X_history, reaction_matrix, cap, k_stoch_func, Nc, t_start, t_end, num_points)
    % Propensity függvények ábrázolása az idő függvényében
    % Bemenetek:
    % - t_history: szimuláció időpontjai
    % - X_history: szimuláció állapotai
    % - reaction_matrix: reakció mátrix
    % - cap: kapacitások
    % - k_stoch_func: átmeneti ráta függvény
    % - t_start, t_end: időintervallum
    % - num_points: mintavételi pontok száma
    
    % Időpontok generálása
    t = linspace(t_start, t_end, num_points);
    num_reactions = size(reaction_matrix, 1);
    propensities = zeros(num_reactions, num_points);
    
    % Propensity értékek kiszámítása
    for i = 1:num_points
        [~, closest_idx] = min(abs(t_history - t(i)));
        current_X = X_history(closest_idx, :);
        
        for r = 1:num_reactions
            from = reaction_matrix(r, 1);
            to = reaction_matrix(r, 2);
            k = k_stoch_func(t(i), from, to);
            propensities(r,i) = k * current_X(from) * (cap(to) - current_X(to));
        end
    end
    
    % Ábrázolás
    figure('Position', [100, 100, 1200, 600]);
    
    % Propensity függvények
    subplot(2,1,1);
    hold on;
    for r = 1:num_reactions
        plot(t, propensities(r,:), 'LineWidth', 2);
    end
    hold off;
    grid on;
    title('Propensity függvények');
    xlabel('Idő [s]');
    ylabel('Propensity érték');
    legend_labels = arrayfun(@(x) sprintf('a_{%d→%d}', ...
        reaction_matrix(x,1), reaction_matrix(x,2)), ...
        1:num_reactions, 'UniformOutput', false);
    legend(legend_labels, 'Location', 'best');
    
    % Átmeneti ráták ábrázolása
    subplot(2,1,2);
    hold on;
    k_values = zeros(num_reactions, num_points);
    for i = 1:num_points
        for r = 1:num_reactions
            from = reaction_matrix(r, 1);
            to = reaction_matrix(r, 2);
            k_values(r,i) = k_stoch_func(t(i), from, to);
        end
    end
    for r = 1:num_reactions
        plot(t, k_values(r,:), 'LineWidth', 2);
    end
    hold off;
    grid on;
    title('Átmeneti ráták (k_{stoch})');
    xlabel('Idő [s]');
    ylabel('k_{stoch} érték');
    legend(legend_labels, 'Location', 'best');
end