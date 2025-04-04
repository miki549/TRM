function [t_history, X_history] = tRSSA(Nc, X, cap, k_stoch_func, reaction_matrix, tfinal)
    % Time-dependent Rejection-based SSA (tRSSA) algoritmus implementációja
    % Bemenetek:
    % - Nc: útszakaszok száma
    % - X: kezdeti állapot (járművek száma)
    % - cap: útszakaszok kapacitása
    % - k_stoch_func: időfüggő sztochasztikus átmeneti ráta függvény
    % - reaction_matrix: reakciók mátrixa
    % - tfinal: szimuláció vége
    
    % 1. Inicializálás
    t = 0;                  % Kezdeti idő
    t_history = [0];        % Időpontok tárolása
    X_history = [X'];       % Állapotok tárolása
    
    % Reakciók száma
    num_internal_reactions = size(reaction_matrix, 1);
    
    % Külső flow-k számának meghatározása
    flows_info = TRM_external_flows(t, Nc);
    inflow_rates = flows_info(:,1);
    outflow_rates = flows_info(:,2);
    
    num_inflows = sum(inflow_rates > 0);
    num_outflows = sum(outflow_rates > 0);
    
    total_reactions = num_internal_reactions + num_inflows + num_outflows;
    
    % 2. Állapottér határok definiálása
    X_min = zeros(Nc, 1);
    X_max = cap;
    
    % 3. Időintervallumok diszkretizálása
    num_intervals = 10;  % k érték a pszeudokódban
    time_points = linspace(0, tfinal, num_intervals + 1);
    
    % 4. Intervallum index inicializálása
    i = 1;
    
    % 5-7. Időfüggő propensity korlátok kiszámítása az első intervallumra
    [aj_lower, aj_upper] = compute_propensity_bounds(X, X_min, X_max, time_points(i), time_points(i+1), k_stoch_func, reaction_matrix, cap, Nc);
    
    % 8. Az összes felső korlát összege
    aj_upper_sum = sum(aj_upper);
    
    % 9. Fő ciklus
    while (t < tfinal)
        % 10. Véletlenszám generálása
        r1 = rand();
        
        % 11-12. Következő reakció idejének kiszámítása és idő frissítése
        tau = -log(r1) / aj_upper_sum;
        t = t + tau;
        
        % 13. Ellenőrizzük, túlléptük-e az aktuális időintervallumot
        if (t > time_points(i+1))
            % 14-15. Idő és intervallum frissítése
            t = time_points(i+1);
            i = i + 1;
            
            % 16-17. Propensity korlátok újraszámítása az új intervallumra
            if (i <= num_intervals)
                [aj_lower, aj_upper] = compute_propensity_bounds(X, X_min, X_max, time_points(i), time_points(i+1), k_stoch_func, reaction_matrix, cap, Nc);
                aj_upper_sum = sum(aj_upper);
            end
            
            % 18. Ugrás vissza a 9. lépésre (folytatjuk a ciklust)
            continue;
        end
        
        % Ha túlléptük a szimuláció végét, kilépünk
        if (t > tfinal)
            break;
        end
        
        % 19-20. Két új véletlenszám generálása
        r2 = rand();
        r3 = rand();
        
        % 21. Minimum index kiválasztása
        cum_prop = cumsum(aj_upper) / aj_upper_sum;
        mu = find(cum_prop >= r2, 1);
        
        % 22. Elfogadás jelző inicializálása
        accepted = false;
        
        % 23-29. Elfogadás-elutasítás teszt
        % Aktuális propensity kiszámítása a kiválasztott reakcióra
        actual_propensity = compute_actual_propensity(X, mu, t, k_stoch_func, reaction_matrix, cap, Nc);
        
        % 23-24. Ellenőrzés a lower bound-dal
        if (r3 <= actual_propensity / aj_upper(mu))
            accepted = true;
        % 25-28. Ellenőrzés a tényleges propensity-vel
        else
            accepted = (r3 <= actual_propensity / aj_upper(mu));
        end
        
        % 31-37. Ha a reakció elfogadásra került
        if (accepted)
            % 32. Állapot frissítése
            if (mu <= num_internal_reactions)
                % Belső átmenet
                from = reaction_matrix(mu, 1);
                to = reaction_matrix(mu, 2);
                X(from) = X(from) - 1;
                X(to) = X(to) + 1;
            else
                % Külső flow-k kezelése
                if (mu <= num_internal_reactions + num_inflows)
                    % Beáramlás
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
                    % Kiáramlás
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
            
            % 33-36. Ellenőrizzük, hogy az új állapot az állapottéren belül van-e
            if (any(X < X_min) || any(X > X_max))
                % 34. Új állapottér határok definiálása
                X_min = min(X_min, X);
                X_max = max(X_max, X);
                
                % 35-36. Propensity korlátok újraszámítása
                [aj_lower, aj_upper] = compute_propensity_bounds(X, X_min, X_max, time_points(i), time_points(i+1), k_stoch_func, reaction_matrix, cap, Nc);
                aj_upper_sum = sum(aj_upper);
            end
            
            % Adatok tárolása
            t_history = [t_history; t];
            X_history = [X_history; X'];
        end
    end
end