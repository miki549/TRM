function [t_history, X_history] = MNRM(Nc, X, cap, k_stoch_func, reaction_matrix, tfinal)
    % Modified Next Reaction Method (MNRM) algoritmus implementációja
    % Bemenetek:
    % - Nc: útszakaszok száma
    % - X: kezdeti állapot (járművek száma)
    % - cap: útszakaszok kapacitása
    % - k_stoch_func: időfüggő sztochasztikus átmeneti ráta függvény (function handle)
    % - reaction_matrix: reakciók mátrixa
    % - tfinal: szimuláció vége
    
    % Inicializálás
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
    
    % MNRM változók inicializálása
    T = zeros(total_reactions, 1);    % Az utolsó reakció óta eltelt idő
    S = -log(rand(total_reactions, 1)); % Következő reakció időpontja
    
    % Fő ciklus
    while t < tfinal
        % Propensity függvények kiszámítása minden reakcióra az aktuális állapotban
        propensities = zeros(total_reactions, 1);
        
        % 1. Belső átmenetek propensity értékei
        for r = 1:num_internal_reactions
            from = reaction_matrix(r, 1);
            to = reaction_matrix(r, 2);
            
            % Időfüggő átmeneti ráta használata
            current_k_stoch = k_stoch_func(t, from, to);
            
            % Propensity: k_stoch(t) * n_i * (c_j - n_j)
            propensities(r) = current_k_stoch * X(from) * (cap(to) - X(to));
        end
        
        % 2. Külső flow-k propensity értékei
        flows_info = TRM_external_flows(t, Nc);
        inflow_rates = flows_info(:,1);
        outflow_rates = flows_info(:,2);
        
        % Beáramlások propensity-jei
        inflow_idx = 0;
        for i = 1:Nc
            if inflow_rates(i) > 0
                inflow_idx = inflow_idx + 1;
                propensities(num_internal_reactions + inflow_idx) = inflow_rates(i) * (cap(i) - X(i));
            end
        end
        
        % Kiáramlások propensity-jei
        outflow_idx = 0;
        for i = 1:Nc
            if outflow_rates(i) > 0
                outflow_idx = outflow_idx + 1;
                propensities(num_internal_reactions + num_inflows + outflow_idx) = outflow_rates(i) * X(i);
            end
        end
        
        % Ha minden propensity nulla, kilépünk
        if all(propensities == 0)
            break;
        end
        
        % Várakozási idők kiszámítása
        tau = zeros(total_reactions, 1);
        integral = zeros(total_reactions,1);
        for j = 1:total_reactions
            if propensities(j) > 0
                % Az integrál egyenlet megoldása a trapezoid_tau függvénnyel
                [tau(j), integral(j)] = trapezoid_tau(S(j), T(j), t, X, j, k_stoch_func, reaction_matrix, cap, Nc, 1e-6);
            else
                tau(j) = Inf;
            end
        end
        
        % Legkorábbi reakció kiválasztása
        [delta_t, mu] = min(tau);
        
        % Idő frissítése
        t = t + delta_t;
        if t > tfinal
            break;
        end
        % EZT A RÉSZT KELL MEGNÉZNI_____________________________________
        % T értékek frissítése minden reakcióra
        for j = 1:total_reactions
            % T frissítése az integrál közelítésével
            if j <= num_internal_reactions
                from = reaction_matrix(j, 1);
                to = reaction_matrix(j, 2);
                start_k = k_stoch_func(t - delta_t, from, to);
                end_k = k_stoch_func(t, from, to);
                avg_k = (start_k + end_k) / 2;
                avg_propensity = avg_k * X(from) * (cap(to) - X(to));
            else
                % Külső flow-k becslése
                avg_propensity = propensities(j);
            end
            
            T(j) = T(j) + avg_propensity * delta_t;
        end
       % IDÁIG______________________________________________________ 
        % Állapot frissítése a választott reakció alapján
        if mu <= num_internal_reactions
            % Belső átmenet két útszakasz között
            from = reaction_matrix(mu, 1);
            to = reaction_matrix(mu, 2);
            X(from) = X(from) - 1;
            X(to) = X(to) + 1;
        else
            % Külső flow-k kezelése
            if mu <= num_internal_reactions + num_inflows
                % Beáramlás
                inflow_idx = mu - num_internal_reactions;
                actual_idx = 0;
                for i = 1:Nc
                    if inflow_rates(i) > 0
                        actual_idx = actual_idx + 1;
                        if actual_idx == inflow_idx
                            X(i) = X(i) + 1;
                            break;
                        end
                    end
                end
            else
                % Kiáramlás
                outflow_idx = mu - num_internal_reactions - num_inflows;
                actual_idx = 0;
                for i = 1:Nc
                    if outflow_rates(i) > 0
                        actual_idx = actual_idx + 1;
                        if actual_idx == outflow_idx
                            X(i) = X(i) - 1;
                            break;
                        end
                    end
                end
            end
        end
        
        % Kiválasztott reakció S értékének frissítése
        S(mu) = S(mu) - log(rand());
        
        % Adatok tárolása
        t_history = [t_history; t];
        X_history = [X_history; X'];
    end
end