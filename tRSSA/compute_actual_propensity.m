function propensity = compute_actual_propensity(X, mu, t, k_stoch_func, reaction_matrix, cap, Nc)
    % Aktuális propensity kiszámítása a kiválasztott reakcióra
    
    % Reakciók száma
    num_internal_reactions = size(reaction_matrix, 1);
    
    % Külső flow-k számának meghatározása
    flows = TRM_external_flows(t, Nc);
    inflow_rates = flows(:,1);
    outflow_rates = flows(:,2);
    
    num_inflows = sum(inflow_rates > 0);
    
    % Propensity kiszámítása a reakció típusa alapján
    if (mu <= num_internal_reactions)
        % Belső átmenet
        from = reaction_matrix(mu, 1);
        to = reaction_matrix(mu, 2);
        
        % Időfüggő átmeneti ráta használata
        k = k_stoch_func(t, from, to);
        
        % Propensity: k_stoch(t) * n_i * (c_j - n_j)
        propensity = k * X(from) * (cap(to) - X(to));
    else
        % Külső flow-k kezelése
        if (mu <= num_internal_reactions + num_inflows)
            % Beáramlás
            inflow_idx = mu - num_internal_reactions;
            actual_idx = 0;
            for i = 1:Nc
                if (inflow_rates(i) > 0)
                    actual_idx = actual_idx + 1;
                    if (actual_idx == inflow_idx)
                        propensity = inflow_rates(i) * (cap(i) - X(i));
                        break;
                    end
                end
            end
        else
            % Kiáramlás
            outflow_idx = mu - num_internal_reactions - num_inflows;
            actual_idx = 0;
            for i = 1:Nc
                if (outflow_rates(i) > 0)
                    actual_idx = actual_idx + 1;
                    if (actual_idx == outflow_idx)
                        propensity = outflow_rates(i) * X(i);
                        break;
                    end
                end
            end
        end
    end
end