function [t_history, X_history] = time_dependent_gillespie(Nc, X, cap, k_stoch_time_fn, reaction_matrix, tfinal, time_intervals, interval_split_num)
% Időfüggő elutasítás-alapú SSA (tRSSA) algoritmus
% Bemenetek:
% - Nc: útszakaszok száma
% - X: kezdeti állapot (járművek száma)
% - cap: útszakaszok kapacitása
% - k_stoch_time_fn: időfüggő sztochasztikus átmeneti ráta függvény, formátum: @(t) k(t)
% - reaction_matrix: reakciók mátrixa
% - tfinal: szimuláció vége
% - time_intervals: az idő diszkretizáció intervallumainak száma (opcionális)

% Alapértelmezett értékek beállítása
if nargin < 7
    time_intervals = 100; % Intervallumok száma
end
if nargin < 8
    interval_split_num=20; % Egy intervallum felosztása
end

% Inicializálás
t = 0;
t_history = 0;
X_history = X';

% Állapothatárok definiálása
X_lower = max(0, X - 10); % Alsó határ minden útszakaszra
X_upper = min(cap, X + 10); % Felső határ minden útszakaszra

% Idő diszkretizálása
time_points = linspace(0, tfinal, time_intervals + 1);
i = 1; % Kezdő időintervallum index

% 4. Reakcióráta határok kiszámítása az első időintervallumra
[c_lower, c_upper] = calculate_rate_bounds(k_stoch_time_fn, time_points(i), time_points(i+1),interval_split_num);

% 5. Propensity határok számítása minden reakcióra
[a_lower, a_upper] = calculate_propensity_bounds(reaction_matrix, c_lower, c_upper, X_lower, X_upper, cap);

% 6. Összesített propensity számítása
a_upper_sum = sum(a_upper);

% Fő ciklus
while t < tfinal
  
    % 7. Várakozási idő generálása
    r1 = rand();
    tau = (-1/a_upper_sum) * log(r1);
    
    % 8. Idő frissítése
    t_next = t + tau;
    
    % 9. Ellenőrizzük, hogy átléptünk-e az időintervallum határán
    if t_next > time_points(i+1)
        % Időintervallum határa, frissíteni kell a rátákat
        t = time_points(i+1);
        i = i + 1;
        
        % Ha elértük az utolsó időintervallumot, kilépünk
        if i > time_intervals
            break;
        end
        
        % Új rátahatárok számítása
        [c_lower, c_upper] = calculate_rate_bounds(k_stoch_time_fn, time_points(i), time_points(i+1),interval_split_num);
        
        % Propensity határok frissítése
        [a_lower, a_upper] = calculate_propensity_bounds(reaction_matrix, c_lower, c_upper, X_lower, X_upper, cap);
        
        % Összesített propensity frissítése
        a_upper_sum = sum(a_upper);
        
        % Adatok tárolása
        t_history = [t_history, t];
        X_history = [X_history; X'];
        
        continue;
    else
        t = t_next;
    end
    
    % Reakció kiválasztása és elutasítás-alapú elfogadása
    r2 = rand();
    r3 = rand();
    
    % Reakció kiválasztása
    cum_prop = cumsum(a_upper) / a_upper_sum;
    mu = find(cum_prop >= r2, 1);
    
    % Elutasítás-alapú elfogadás
    accepted = false;
    
    % Gyors út: ha a reakció biztosan elfogadható az alsó határral
    if r3 <= (a_lower(mu) / a_upper(mu))
        accepted = true;
    else
        % Pontos propensity számítása a jelenlegi állapotban
        a_current = calculate_exact_propensity(mu, reaction_matrix, X, cap, k_stoch_time_fn(t));
        
        % Elfogadás a pontos propensity alapján
        if r3 <= (a_current / a_upper(mu))
            accepted = true;
        end
    end
    
    % Ha a reakció elfogadott, frissítjük az állapotot
    if accepted
        % Reakció végrehajtása
        if mu <= size(reaction_matrix, 1)
            % Belső átmenet két útszakasz között
            from = reaction_matrix(mu, 1);
            to = reaction_matrix(mu, 2);
            X(from) = X(from) - 1;
            X(to) = X(to) + 1;
        end
        
        % Ellenőrizzük, hogy az új állapot a határokon belül van-e
        if any(X < X_lower) || any(X > X_upper)
            % Ha kilépett a határokból, újraszámoljuk a határokat
            X_lower = max(0, X - 10);
            X_upper = min(cap, X + 10);
            
            % Propensity határok újraszámítása
            [a_lower, a_upper] = calculate_propensity_bounds(reaction_matrix, c_lower, c_upper, X_lower, X_upper, cap);
            
            % Összesített propensity frissítése
            a_upper_sum = sum(a_upper);
        end
    end
    
    % Külső flow kezelése (ha szükséges)
    flows = TRM_external_flows(t, Nc);
    inflow_rates = flows(:,1);
    outflow_rates = flows(:,2);
    
    % Külső flow-k végrehajtása (determinisztikus módon vagy propensity alapon)
    for node = 1:Nc
        % Inflow kezelése
        if inflow_rates(node) > 0 && X(node) < cap(node)
            p_in = inflow_rates(node) * tau * (cap(node) - X(node));
            if rand() < p_in
                X(node) = X(node) + 1;
            end
        end
        
        % Outflow kezelése
        if outflow_rates(node) > 0 && X(node) > 0
            p_out = outflow_rates(node) * tau * X(node);
            if rand() < p_out
                X(node) = X(node) - 1;
            end
        end
    end
    
    % Adatok tárolása
    t_history = [t_history, t];
    X_history = [X_history; X'];
    
    % Ha elértük az időlimitet, kilépünk
    if t >= tfinal
        break;
    end
end
end
