function propensity = compute_actual_propensity(X, mu, t, k_stoch_func, reaction_matrix, cap, Nc)
    % Aktuális propensity kiszámítása a kiválasztott reakcióra
    
    % Áramlások kiszámítása és cache-elése
    flows = TRM_external_flows(t, Nc);
    inflow_rates = flows(:,1);
    outflow_rates = flows(:,2);
    
    % Reakciók számának meghatározása
    num_internal_reactions = size(reaction_matrix, 1);
    num_inflows = sum(inflow_rates > 0);
    
    % Propensity számítása a reakció típusa alapján
    if mu <= num_internal_reactions
        % Belső átmenetek
        from = reaction_matrix(mu, 1);
        to = reaction_matrix(mu, 2);
        propensity = k_stoch_func(t, from, to) * X(from) * (cap(to) - X(to));
    else
        if mu <= num_internal_reactions + num_inflows
            % Beáramlások
            inflow_idx = mu - num_internal_reactions;
            target_idx = find(inflow_rates > 0, inflow_idx, 'last');
            target_idx = target_idx(1);
            propensity = inflow_rates(target_idx) * (cap(target_idx) - X(target_idx));
        else
            % Kiáramlások
            outflow_idx = mu - (num_internal_reactions + num_inflows);
            target_idx = find(outflow_rates > 0, outflow_idx, 'last');
            target_idx = target_idx(1);
            propensity = outflow_rates(target_idx) * X(target_idx);
        end
    end
end