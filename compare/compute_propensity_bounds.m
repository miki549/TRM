function [aj_lower, aj_upper] = compute_propensity_bounds(X, X_min, X_max, t_start, t_end, k_stoch_func, reaction_matrix, cap, Nc)
    % Propensity alsó és felső korlátainak kiszámítása egy időintervallumra
    
    % Reakciók számának meghatározása
    num_internal_reactions = size(reaction_matrix, 1);
    
    % Külső áramlások meghatározása az intervallum elején és végén
    flows_start = TRM_external_flows(t_start, Nc);
    flows_end = TRM_external_flows(t_end, Nc);
    
    inflow_rates_start = flows_start(:,1);
    outflow_rates_start = flows_start(:,2);
    inflow_rates_end = flows_end(:,1);
    outflow_rates_end = flows_end(:,2);
    
    num_inflows = sum(max(inflow_rates_start, inflow_rates_end) > 0);
    num_outflows = sum(max(outflow_rates_start, outflow_rates_end) > 0);
    
    total_reactions = num_internal_reactions + num_inflows + num_outflows;
    
    % Korlátok inicializálása
    aj_lower = zeros(total_reactions, 1);
    aj_upper = zeros(total_reactions, 1);
    
    % Időpontok mintavételezése az intervallumon belül
    num_samples = 5;
    t_samples = linspace(t_start, t_end, num_samples);
    
    % Belső átmenetek propensity korlátainak számítása
    for r = 1:num_internal_reactions
        from = reaction_matrix(r, 1);
        to = reaction_matrix(r, 2);
        
        prop_samples = zeros(num_samples, 1);
        
        for s = 1:num_samples
            t = t_samples(s);
            k = k_stoch_func(t, from, to);
            
            % Alsó korlát számítása
            min_n_from = max(X_min(from), 0);
            max_n_to = min(X_max(to), cap(to));
            min_available = max(cap(to) - max_n_to, 0);
            
            % Felső korlát számítása
            max_n_from = min(X_max(from), cap(from));
            min_n_to = max(X_min(to), 0);
            max_available = max(cap(to) - min_n_to, 0);
            
            % Aktuális propensity
            current_prop = k * X(from) * (cap(to) - X(to));
            prop_samples(s) = current_prop;
        end
        
        aj_lower(r) = min(prop_samples);
        aj_upper(r) = max(prop_samples);
    end
    
    % Beáramlások propensity korlátainak számítása
    inflow_idx = 0;
    for i = 1:Nc
        if (max(inflow_rates_start(i), inflow_rates_end(i)) > 0)
            inflow_idx = inflow_idx + 1;
            
            prop_samples = zeros(num_samples, 1);
            for s = 1:num_samples
                t = t_samples(s);
                flows = TRM_external_flows(t, Nc);
                inflow_rate = flows(i, 1);
                current_prop = inflow_rate * (cap(i) - X(i));
                prop_samples(s) = current_prop;
            end
            
            aj_lower(num_internal_reactions + inflow_idx) = min(prop_samples);
            aj_upper(num_internal_reactions + inflow_idx) = max(prop_samples);
        end
    end
    
    % Kiáramlások propensity korlátainak számítása
    outflow_idx = 0;
    for i = 1:Nc
        if (max(outflow_rates_start(i), outflow_rates_end(i)) > 0)
            outflow_idx = outflow_idx + 1;
            
            prop_samples = zeros(num_samples, 1);
            for s = 1:num_samples
                t = t_samples(s);
                flows = TRM_external_flows(t, Nc);
                outflow_rate = flows(i, 2);
                current_prop = outflow_rate * X(i);
                prop_samples(s) = current_prop;
            end
            
            aj_lower(num_internal_reactions + num_inflows + outflow_idx) = min(prop_samples);
            aj_upper(num_internal_reactions + num_inflows + outflow_idx) = max(prop_samples);
        end
    end
    
    % Numerikus stabilitás biztosítása
    aj_lower(aj_lower < 1e-10) = 0;
    aj_upper(aj_upper < 1e-10) = 1e-10;
end