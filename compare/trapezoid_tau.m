function [tau, final_integral] = trapezoid_tau(S, T, t, X, mu, k_stoch_func, reaction_matrix, cap, Nc, tolerance)
    % Trapéz szabállyal történő tau számítás az MNRM algoritmushoz
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
    
    % Kezdeti propensity kiszámítása
    initial_propensity = compute_actual_propensity(X, mu, t, k_stoch_func, reaction_matrix, cap, Nc);
    
    if initial_propensity == 0
        tau = Inf;
        final_integral = 0;
        return;
    end
    
    % Kezdeti becslés a tau-ra
    target_area = S - T;
    tau_guess = target_area / initial_propensity;
    
    % Iterációs paraméterek
    max_iterations = 200;
    iteration = 0;
    tau = tau_guess;
    
    while iteration < max_iterations
        % Trapéz szabály alkalmazása
        start_prop = compute_actual_propensity(X, mu, t, k_stoch_func, reaction_matrix, cap, Nc);
        end_prop = compute_actual_propensity(X, mu, t + tau, k_stoch_func, reaction_matrix, cap, Nc);
        
        integral_value = (tau/2) * (start_prop + end_prop);
        error = abs(integral_value - target_area);
        
        if error < tolerance
            final_integral = integral_value;
            return;
        end
        
        % Tau frissítése
        if integral_value > 0
            tau = max(0.0, tau * target_area / integral_value);
        else
            tau = Inf;
            final_integral = 0;
            return;
        end
        
        iteration = iteration + 1;
    end
    %fprintf('S-T: %f Computed: %f',target_area,integral_value);
    %fprintf('\n');
    final_integral = integral_value;
end 