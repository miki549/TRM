function [aj_lower, aj_upper] = compute_propensity_bounds(X, X_min, X_max, t_start, t_end, k_stoch_func, reaction_matrix, cap, Nc)
    % Propensity alsó és felső korlátainak kiszámítása egy időintervallumra
    
    % Reakciók száma
    num_internal_reactions = size(reaction_matrix, 1);
    
    % Külső flow-k számának meghatározása
    flows_start = TRM_external_flows(t_start, Nc);
    flows_end = TRM_external_flows(t_end, Nc);
    
    inflow_rates_start = flows_start(:,1);
    outflow_rates_start = flows_start(:,2);
    inflow_rates_end = flows_end(:,1);
    outflow_rates_end = flows_end(:,2);
    
    num_inflows = sum(max(inflow_rates_start, inflow_rates_end) > 0);
    num_outflows = sum(max(outflow_rates_start, outflow_rates_end) > 0);
    
    total_reactions = num_internal_reactions + num_inflows + num_outflows;
    
    % Propensity korlátok inicializálása
    aj_lower = zeros(total_reactions, 1);
    aj_upper = zeros(total_reactions, 1);
    
    % Időpontok az intervallumon belül mintavételezéshez
    num_samples = 5;
    t_samples = linspace(t_start, t_end, num_samples);
    
    % 1. Belső átmenetek propensity korlátai
    for r = 1:num_internal_reactions
        from = reaction_matrix(r, 1);
        to = reaction_matrix(r, 2);
        
        prop_samples = zeros(num_samples, 1);
        
        % Minden időpontra kiszámítjuk a propensity-t
        for s = 1:num_samples
            t = t_samples(s);
            k = k_stoch_func(t, from, to);
            
            % Alsó korlát: min n_i * min (c_j - n_j) * min k
            min_n_from = max(X_min(from), 0);
            max_n_to = min(X_max(to), cap(to));
            min_available = max(cap(to) - max_n_to, 0);
            
            % Felső korlát: max n_i * max (c_j - n_j) * max k
            max_n_from = min(X_max(from), cap(from));
            min_n_to = max(X_min(to), 0);
            max_available = max(cap(to) - min_n_to, 0);
            
            % Aktuális állapotra vonatkozó propensity
            current_prop = k * X(from) * (cap(to) - X(to));
            
            % Mintavétel tárolása
            prop_samples(s) = current_prop;
        end
        
        % Alsó és felső korlátok az intervallumon
        aj_lower(r) = min(prop_samples);
        aj_upper(r) = max(prop_samples);
    end
    
    % 2. Külső flow-k propensity korlátai
    % Beáramlások
    inflow_idx = 0;
    for i = 1:Nc
        if (max(inflow_rates_start(i), inflow_rates_end(i)) > 0)
            inflow_idx = inflow_idx + 1;
            
            % Minden időpontra kiszámítjuk a propensity-t
            prop_samples = zeros(num_samples, 1);
            for s = 1:num_samples
                t = t_samples(s);
                flows = TRM_external_flows(t, Nc);
                inflow_rate = flows(i, 1);
                
                % Aktuális állapotra vonatkozó propensity
                current_prop = inflow_rate * (cap(i) - X(i));
                
                % Mintavétel tárolása
                prop_samples(s) = current_prop;
            end
            
            % Alsó és felső korlátok az intervallumon
            aj_lower(num_internal_reactions + inflow_idx) = min(prop_samples);
            aj_upper(num_internal_reactions + inflow_idx) = max(prop_samples);
        end
    end
    
    % Kiáramlások
    outflow_idx = 0;
    for i = 1:Nc
        if (max(outflow_rates_start(i), outflow_rates_end(i)) > 0)
            outflow_idx = outflow_idx + 1;
            
            % Minden időpontra kiszámítjuk a propensity-t
            prop_samples = zeros(num_samples, 1);
            for s = 1:num_samples
                t = t_samples(s);
                flows = TRM_external_flows(t, Nc);
                outflow_rate = flows(i, 2);
                
                % Aktuális állapotra vonatkozó propensity
                current_prop = outflow_rate * X(i);
                
                % Mintavétel tárolása
                prop_samples(s) = current_prop;
            end
            
            % Alsó és felső korlátok az intervallumon
            aj_lower(num_internal_reactions + num_inflows + outflow_idx) = min(prop_samples);
            aj_upper(num_internal_reactions + num_inflows + outflow_idx) = max(prop_samples);
        end
    end
    
    % Nullára állítjuk a nagyon kis értékeket az alsó korlátnál
    aj_lower(aj_lower < 1e-10) = 0;
    
    % Kis pozitív értéket adunk a nulla felső korlátoknak a numerikus stabilitás érdekében
    aj_upper(aj_upper < 1e-10) = 1e-10;
end