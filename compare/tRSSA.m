function [t_history, X_history] = tRSSA(Nc, X, cap, k_stoch_func, reaction_matrix, tfinal)
    % Time-dependent Rejection-based SSA (tRSSA) algoritmus implementációja
    % Bemenetek:
    % - Nc: útszakaszok száma
    % - X: kezdeti állapot (járművek száma)
    % - cap: útszakaszok kapacitása
    % - k_stoch_func: időfüggő sztochasztikus átmeneti ráta függvény
    % - reaction_matrix: reakciók mátrixa
    % - tfinal: szimuláció vége
    
    % Inicializálás
    t = 0;
    t_history = [0];
    X_history = [X'];
    
    % Reakciók számának meghatározása
    num_internal_reactions = size(reaction_matrix, 1);
    
    % Külső áramlások meghatározása
    flows_info = TRM_external_flows(t, Nc);
    inflow_rates = flows_info(:,1);
    outflow_rates = flows_info(:,2);
    num_inflows = sum(inflow_rates > 0);
    
    % Állapottér határok definiálása
    X_min = zeros(Nc, 1);
    X_max = cap;
    
    % Időintervallumok diszkretizálása
    num_intervals = 10;
    time_points = linspace(0, tfinal, num_intervals + 1);
    
    % Intervallum index inicializálása
    i = 1;
    
    % Propensity korlátok számítása az első intervallumra
    [aj_lower, aj_upper] = compute_propensity_bounds(X, X_min, X_max, time_points(i), time_points(i+1), k_stoch_func, reaction_matrix, cap, Nc);
    aj_upper_sum = sum(aj_upper);
    
    % Fő ciklus
    while (t < tfinal)
        % Következő reakció idejének kiszámítása
        r1 = rand();
        tau = -log(r1) / aj_upper_sum;
        t = t + tau;
        
        % Időintervallum ellenőrzése
        if (t > time_points(i+1))
            t = time_points(i+1);
            i = i + 1;
            
            if (i <= num_intervals)
                [aj_lower, aj_upper] = compute_propensity_bounds(X, X_min, X_max, time_points(i), time_points(i+1), k_stoch_func, reaction_matrix, cap, Nc);
                aj_upper_sum = sum(aj_upper);
            end
            continue;
        end
        
        if (t > tfinal)
            break;
        end
        
        % Reakció kiválasztása és elfogadás-elutasítás teszt
        r2 = rand();
        r3 = rand();
        
        cum_prop = cumsum(aj_upper) / aj_upper_sum;
        mu = find(cum_prop >= r2, 1);
        
        accepted = false;
        actual_propensity = compute_actual_propensity(X, mu, t, k_stoch_func, reaction_matrix, cap, Nc);
        
        if (r3 <= actual_propensity / aj_upper(mu))
            accepted = true;
        end
        
        % Állapot frissítése elfogadott reakció esetén
        if (accepted)
            if (mu <= num_internal_reactions)
                from = reaction_matrix(mu, 1);
                to = reaction_matrix(mu, 2);
                X(from) = X(from) - 1;
                X(to) = X(to) + 1;
            else
                if (mu <= num_internal_reactions + num_inflows)
                    inflow_idx = mu - num_internal_reactions;
                    actual_idx = 0;
                    for idx = 1:Nc
                        if (inflow_rates(idx) > 0)
                            actual_idx = actual_idx + 1;
                            if (actual_idx == inflow_idx)
                                X(idx) = X(idx) + 1;
                                break;
                            end
                        end
                    end
                else
                    outflow_idx = mu - num_internal_reactions - num_inflows;
                    actual_idx = 0;
                    for idx = 1:Nc
                        if (outflow_rates(idx) > 0)
                            actual_idx = actual_idx + 1;
                            if (actual_idx == outflow_idx)
                                X(idx) = X(idx) - 1;
                                break;
                            end
                        end
                    end
                end
            end
            
            % Állapottér határok frissítése ha szükséges
            if (any(X < X_min) || any(X > X_max))
                X_min = min(X_min, X);
                X_max = max(X_max, X);
                [aj_lower, aj_upper] = compute_propensity_bounds(X, X_min, X_max, time_points(i), time_points(i+1), k_stoch_func, reaction_matrix, cap, Nc);
                aj_upper_sum = sum(aj_upper);
            end
            
            % Eredmények tárolása
            t_history = [t_history; t];
            X_history = [X_history; X'];
        end
    end
end