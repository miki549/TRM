% Segédfüggvény: Pontos propensity számítása
function a_j = calculate_exact_propensity(reaction_idx, reaction_matrix, X, cap, k_stoch_current)
% Egy adott reakcióra kiszámítjuk a pontos propensity-t
if reaction_idx <= size(reaction_matrix, 1)
    from = reaction_matrix(reaction_idx, 1);
    to = reaction_matrix(reaction_idx, 2);
    a_j = k_stoch_current * X(from) * (cap(to) - X(to));
else
    a_j = 0; % Ha nem belső átmenet, akkor 0
end
end