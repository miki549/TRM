function dy = TRM_red_ode_mod(t, y, K, cap, k_stoch_func)
    % Általánosított riboszóma áramlási modell ODE implementációja
    % Bemenetek:
    % - t: idő
    % - y: állapotvektor [n(1) ... n(m)]
    % - K: alap kapcsolódási mátrix (konstans)
    % - cap: kapacitásvektor
    % - k_stoch_func: időfüggő átmeneti ráta függvény

    m = size(K,1);    % Cellák száma
    n = y(1:m);       % Járművek száma

    % Monomális mátrix összeállítása
    NS = zeros(m,m);
    for i = 1:m
        for j = 1:m
            NS(i,j) = n(i)*(cap(j)-n(j));
        end
    end

    tmp_dy = zeros(m, 1);

    % Időfüggő átmeneti rátákkal történő egyenletrendszer felépítése
    for i = 1:m
        outflow = 0;
        inflow = 0;
        
        % Kiáramlások számítása az i-edik cellából
        for j = 1:m
            if K(i,j) > 0
                rate = k_stoch_func(t, i, j);
                outflow = outflow + rate * n(i) * (cap(j) - n(j));
            end
        end
        
        % Beáramlások számítása az i-edik cellába
        for j = 1:m
            if K(j,i) > 0
                rate = k_stoch_func(t, j, i);
                inflow = inflow + rate * n(j) * (cap(i) - n(i));
            end
        end
        
        tmp_dy(i) = inflow - outflow;
    end

    % Külső áramlások hozzáadása
    ext_flows = TRM_external_flows(t,m);
    ext_inflow = ext_flows(:,1).*(cap-y);
    ext_outflow = ext_flows(:,2).*y;

    dy = tmp_dy + ext_inflow - ext_outflow;
end



