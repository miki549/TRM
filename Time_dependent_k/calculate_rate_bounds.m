% Segédfüggvény: Időfüggő reakcióráta határok számítása
function [c_lower, c_upper] = calculate_rate_bounds(k_stoch_time_fn, t_start, t_end, split_num)
% Időintervallumban a minimum és maximum k_stoch értékek keresése
% Egyszerű módszer: diszkretizáljuk az intervallumot és keressük a min/max-ot
time_points = linspace(t_start, t_end, split_num);
k_values = zeros(size(time_points));

for i = 1:length(time_points)
    k_values(i) = k_stoch_time_fn(time_points(i));
end

c_lower = min(k_values);
c_upper = max(k_values);
end