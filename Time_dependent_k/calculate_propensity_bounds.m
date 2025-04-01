% Segédfüggvény: Propensity határok számítása
function [a_lower, a_upper] = calculate_propensity_bounds(reaction_matrix, c_lower, c_upper, X_lower, X_upper, cap)
% Minden reakcióra kiszámítjuk a propensity határokat
num_reactions = size(reaction_matrix, 1);
a_lower = zeros(num_reactions, 1);
a_upper = zeros(num_reactions, 1);
h_lower = zeros(num_reactions, 1);
h_upper = zeros(num_reactions, 1);

for r = 1:num_reactions
    from = reaction_matrix(r, 1);
    to = reaction_matrix(r, 2);
    
    % Állapotfüggvény határok (h_j)
    h_lower(r) = X_lower(from) * (cap(to) - X_upper(to));
    h_upper(r) = X_upper(from) * (cap(to) - X_lower(to));
    
    % Propensity határok (a_j)
    a_lower(r) = c_lower * h_lower(r);
    a_upper(r) = c_upper * h_upper(r);
    
    % Biztosítjuk, hogy a határok pozitívak legyenek
    a_lower(r) = max(0, a_lower(r));
    a_upper(r) = max(0, a_upper(r));
end
end