function out = TRM_external_flows(t, N)
    % Külső be- és kiáramlások definiálása a hálózati TRM modellhez
    % Bemenetek:
    % - t: idő
    % - N: útszakaszok száma
    %
    % Kimenet:
    % - out first column: beáramlások, second column: kiáramlások
    % - out(i,1): külső beáramlási ráta az i-edik útszakaszra
    % - out(i,2): kiáramlási ráta az i-edik útszakaszról

    flow = zeros(N, 2);

    % Beáramlási ráták (időfüggő)
    flow(1, 1) = 0.05 * (sin(t/50) + 1.5);  % Első útszakaszra változó beáramlás
    
    % Reggeli és délutáni csúcsidőszak szimulálása
    %morning_peak = exp(-((t-100)^2)/1000) * 0.1;  % Reggeli csúcs kb. t=100 körül
    %evening_peak = exp(-((t-300)^2)/1000) * 0.08; % Délutáni csúcs kb. t=300 körül
    %flow(2, 1) = morning_peak + evening_peak;     % Második útszakaszra csúcsidőszaki beáramlás
    
    % Kiáramlási ráták (időfüggő)
    flow(3, 2) = 0.03 * (cos(t/70) + 1.2);  % Harmadik útszakaszról változó kiáramlás
    
    % Speciális esemény: útlezárás/baleset szimulálása
    %if t > 200 && t < 250
        % Harmadik útszakasz kiáramlási rátájának növelése (pl. terelés miatt)
    %    flow(3, 2) = 0.06;
    %else
    %    flow(3, 2) = 0.02;
    %end

    out = flow;
end