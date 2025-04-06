% Egyszerű teszteset az analitikus ellenőrzéshez
function test_analytical()
    % Egyszerű eset: konstans propensity függvény
    % Ekkor az integrál egyszerűen: a * τ = S - T
    % Tehát τ = (S - T) / a
    
    % Teszt paraméterek
    S = 1.0;
    T = 0.1;
    t = 0;
    X = [1];
    mu = 1;
    
    % Egyszerű konstans átmeneti ráta függvény
    k_stoch_func = @(t, from, to) 0.5;  % Konstans átmeneti ráta
    
    % Reakcómátrix és kapacitás beállítása
    reaction_matrix = [1, 1];
    cap = [30];
    Nc = 1;
    tolerance = 1e-8;
    
    % Numerikus megoldás
    [tau_numerical, steps] = adaptive_rk4_solver(S, T, t, X, mu, k_stoch_func, reaction_matrix, cap, Nc, tolerance);
    
    % Analitikus megoldás
    % A propensity = k * X(from) * (cap(to) - X(to))
    % Ahol k = 0.1, X(from) = 1, cap(to) = 10, X(to) = 1
    a = 0.1 * 1 * (10 - 1);  % Konstans propensity
    tau_analytical = (S - T) / a;
    
    % Eredmények kiírása
    fprintf('Numerikus megoldás: τ = %f (lépések száma: %d)\n', tau_numerical, steps);
    fprintf('Analitikus megoldás: τ = %f\n', tau_analytical);
    fprintf('Különbség: %e\n', abs(tau_numerical - tau_analytical));
end

% Konvergencia teszt
function test_convergence()
    % Teszt paraméterek
    S = 1.0;
    T = 0.0;
    t = 0;
    X = [1];
    mu = 1;
    k_stoch_func = @(t, from, to) 0.1 * (1 + sin(t));  % Időfüggő átmeneti ráta
    reaction_matrix = [1, 1];
    cap = [10];
    Nc = 1;
    
    % Különböző toleranciák tesztelése
    tolerances = [1e-4, 1e-6, 1e-8, 1e-10];
    results = zeros(size(tolerances));
    
    fprintf('Konvergencia teszt:\n');
    fprintf('Tolerancia\tτ érték\t\tLépések száma\n');
    
    for i = 1:length(tolerances)
        [tau, steps] = adaptive_rk4_solver(S, T, t, X, mu, k_stoch_func, reaction_matrix, cap, Nc, tolerances(i));
        results(i) = tau;
        fprintf('%e\t%f\t%d\n', tolerances(i), tau, steps);
    end
    
    % Konvergencia ellenőrzése
    differences = abs(diff(results));
    fprintf('\nKülönbségek a szomszédos eredmények között:\n');
    for i = 1:length(differences)
        fprintf('|τ%d - τ%d| = %e\n', i, i+1, differences(i));
    end
end

test_analytical();
test_convergence();