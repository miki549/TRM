function tau = adaptive_rk4_solver(S, T, t, X, mu, k_stoch_func, reaction_matrix, cap, Nc, tolerance)
    % Az adaptív 4. rendű Runge-Kutta módszer implementációja az MNRM algoritmushoz
    % 
    % Bemeneti paraméterek:
    % S - A következő reakció időpontja (az MNRM algoritmusból)
    % T - Az utolsó reakció óta eltelt idő
    % t - Aktuális szimulációs idő
    % X - A rendszer állapotvektora (járművek száma útszakaszonként)
    % mu - Az aktuális reakció indexe
    % k_stoch_func - Időfüggő átmeneti ráta függvény
    % reaction_matrix - Reakciómátrix (útszakaszok közötti kapcsolatok)
    % cap - Útszakaszok kapacitásvektora
    % Nc - Útszakaszok száma
    % tolerance - Megengedett numerikus hiba
    
    % Kezdeti propensity kiszámítása a kezdeti lépésköz becsléséhez
    initial_propensity = compute_actual_propensity(X, mu, t, k_stoch_func, reaction_matrix, cap, Nc);
    % Kezdeti lépésköz becslése - kisebb kezdeti lépésköz a pontosabb integráláshoz
    initial_h = (S - T) / initial_propensity / 10;
    
    % Lépésköz korlátok beállítása
    h_min = initial_h * 1e-6;  % Minimális lépésköz a túl kis lépések elkerülésére
    h_max = initial_h * 10;    % Maximális lépésköz a pontosság megőrzésére
    
    % Kezdeti értékek inicializálása
    tau = 0;              % A teljes integrálási idő
    integral_sum = 0;     % Az integrál összege
    h = initial_h;        % Aktuális lépésköz
    
    % Fő integrálási ciklus
    while integral_sum < (S - T)  % Addig folytatjuk, amíg el nem érjük a kívánt integrálértéket
        % RK4 lépés végrehajtása és hibabecslés
        [increment1, error_est] = rk4_step(t, tau, X, mu, h, k_stoch_func, reaction_matrix, cap, Nc);
        
        % Adaptív lépésköz szabályozás a hibabecslés alapján
        if error_est > tolerance
            % Ha a hiba túl nagy, csökkentjük a lépésközt
            % A 0.9 biztonsági faktor, a 0.2 empirikus konstans
            h = max(h_min, h * 0.9 * (tolerance/error_est)^0.2);
            continue;  % Újrapróbáljuk az aktuális lépést kisebb lépésközzel
        elseif error_est < tolerance/10
            % Ha a hiba sokkal kisebb a megengedettnél, növeljük a lépésközt
            h = min(h_max, h * 1.1);  % Maximum 10%-os növelés
        end
        
        % Az aktuális lépés hozzáadása az integrálhoz
        integral_sum = integral_sum + increment1;
        tau = tau + h;
        
        % Túllövés ellenőrzése és korrekciója
        if integral_sum > (S - T)
            % Ha túlléptük a célt, lineáris interpolációval korrigálunk
            excess = integral_sum - (S - T);
            tau = tau - h * (excess/increment1);  % Arányos visszalépés
            break;
        end
    end
end

function [increment, error_est] = rk4_step(t, tau, X, mu, h, k_stoch_func, reaction_matrix, cap, Nc)
    % A klasszikus 4. rendű Runge-Kutta lépés implementációja
    % 
    % Bemeneti paraméterek megegyeznek a fő függvénnyel
    % h - Az aktuális lépésköz
    %
    % Kimeneti paraméterek:
    % increment - Az integrál növekménye ebben a lépésben
    % error_est - A numerikus hiba becslése
    
    % A négy RK4 együttható kiszámítása különböző időpontokban
    k1 = compute_actual_propensity(X, mu, t + tau, k_stoch_func, reaction_matrix, cap, Nc);
    k2 = compute_actual_propensity(X, mu, t + tau + h/2, k_stoch_func, reaction_matrix, cap, Nc);
    k3 = compute_actual_propensity(X, mu, t + tau + h/2, k_stoch_func, reaction_matrix, cap, Nc);
    k4 = compute_actual_propensity(X, mu, t + tau + h, k_stoch_func, reaction_matrix, cap, Nc);
    
    % Az RK4 formula alkalmazása: (k1 + 2k2 + 2k3 + k4)/6
    increment = (h/6) * (k1 + 2*k2 + 2*k3 + k4);
    
    % Hibabecslés egy alacsonyabb rendű módszerrel való összehasonlítással
    % Itt a trapéz szabályt használjuk az összehasonlításhoz
    k_lower = (h/2) * (k1 + k4);  % Trapéz szabály
    error_est = abs(increment - k_lower);  % A két módszer különbsége adja a hibabecslést
end