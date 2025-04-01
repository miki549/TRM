function [t_history, X_history] = gillespie_simulation(Nc, X, cap, k_stoch, reaction_matrix, tfinal)
    % Gillespie algoritmus szimulációja
    % Bemenetek:
    % - Nc: útszakaszok száma
    % - X: kezdeti állapot (járművek száma)
    % - cap: útszakaszok kapacitása
    % - k_stoch: sztochasztikus átmeneti ráta
    % - reaction_matrix: reakciók mátrixa
    % - tfinal: szimuláció vége
    
    % Inicializálás
    t = 0; % Szimuláció kezdete
    t_history = 0; % Időpontok tárolása
    X_history = X'; % Állapotok tárolása

    j = 0; % Reakció számláló

    % Gillespie algoritmus fő ciklusa
    while t < tfinal
        j = j + 1;
        
        % 1. lépés: Propensity/intenzitás függvények kiszámítása minden reakcióra
        propensities = zeros(size(reaction_matrix, 1), 1);
        
        for r = 1:size(reaction_matrix, 1)
            from = reaction_matrix(r, 1);
            to = reaction_matrix(r, 2);
            
            % Propensity: k_stoch * n_i * (c_j - n_j)
            propensities(r) = k_stoch * X(from) * (cap(to) - X(to));
        end
        
        % Külső flow hozzáadása
        flows = TRM_external_flows(t, Nc);
        inflow_rates = flows(:,1);
        outflow_rates = flows(:,2);
        
        % Külső flow propensity-jei
        inflow_propensities = inflow_rates .* (cap - X);
        outflow_propensities = outflow_rates .* X;
        
        % Összes reakció propensity-jének összegyűjtése
        all_propensities = [propensities; inflow_propensities; outflow_propensities];
        total_propensity = sum(all_propensities);
        
        % Ha nincs lehetséges reakció, kilép
        if total_propensity == 0
            break;
        end
        
        % 2-3. lépés: Két véletlenszám generálása
        r1 = rand();
        r2 = rand();
        
        % 4. lépés: Várakozási idő kiszámítása
        delta_t = log(1/r1) / total_propensity;
        
        % Idő frissítése
        t = t + delta_t;
        if t > tfinal
            break;
        end
        
        % 5. lépés: Melyik reakció következik be?
        cum_prop = cumsum(all_propensities) / total_propensity;
        mu = find(cum_prop >= r2, 1);
        
        % Állapot frissítése a választott reakció alapján
        if mu <= size(reaction_matrix, 1)
            % Belső átmenet két útszakasz között
            from = reaction_matrix(mu, 1);
            to = reaction_matrix(mu, 2);
            X(from) = X(from) - 1;
            X(to) = X(to) + 1;
        elseif mu <= size(reaction_matrix, 1) + Nc
            % Külső inflow
            inflow_idx = mu - size(reaction_matrix, 1);
            X(inflow_idx) = X(inflow_idx) + 1;
        else
            % Külső outflow
            outflow_idx = mu - size(reaction_matrix, 1) - Nc;
            X(outflow_idx) = X(outflow_idx) - 1;
        end
        
        % Adatok tárolása
        t_history = [t_history, t];
        X_history = [X_history; X'];
    end
end