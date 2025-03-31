function [t_history, X_history] = gillespie_time_dependent(Nc, X, cap, k_stoch_func, reaction_matrix, tfinal)
    % Időfüggő Gillespie algoritmus szimulációja
    % Bemenetek:
    % - Nc: útszakaszok száma
    % - X: kezdeti állapot (járművek száma)
    % - cap: útszakaszok kapacitása
    % - k_stoch_func: időfüggő sztochasztikus átmeneti ráta függvény
    % - reaction_matrix: reakciók mátrixa
    % - tfinal: szimuláció vége
    
    % Inicializálás
    t = 0; % Szimuláció kezdete
    t_history = 0; % Időpontok tárolása
    X_history = X'; % Állapotok tárolása

    j = 0; % Reakció számláló
    
    % Időintervallumok diszkretizációja a tRSSA alapján
    % A szimulációs időt [0, tfinal] felosztjuk kisebb intervallumokra
    time_intervals = 0:tfinal/20:tfinal;  % 20 intervallumra osztjuk
    current_interval = 1;

    % Gillespie algoritmus fő ciklusa
    while t < tfinal
        j = j + 1;
        
        % Aktuális időfüggő k_stoch érték meghatározása
        current_k_stoch = k_stoch_func(t);
        
        % 1. lépés: Propensity/intenzitás függvények kiszámítása minden reakcióra
        propensities = zeros(size(reaction_matrix, 1), 1);
        
        for r = 1:size(reaction_matrix, 1)
            from = reaction_matrix(r, 1);
            to = reaction_matrix(r, 2);
            
            % Időfüggő propensity: current_k_stoch(t) * n_i * (c_j - n_j)
            propensities(r) = current_k_stoch * X(from) * (cap(to) - X(to));
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
        t_new = t + delta_t;
        
        % Ellenőrizzük, hogy átléptünk-e egy új időintervallumba
        if current_interval < length(time_intervals) && t_new > time_intervals(current_interval + 1)
            % Ha a következő reakció ideje átlépné az aktuális időintervallumot
            % az időt az időintervallum végére állítjuk
            t = time_intervals(current_interval + 1);
            current_interval = current_interval + 1;
            
            % Nem történik reakció, csak az időt frissítjük
            % és a következő iterációban új propensityk számolódnak
            t_history = [t_history, t];
            X_history = [X_history; X'];
            continue;
        else
            t = t_new;
            if t > tfinal
                break;
            end
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