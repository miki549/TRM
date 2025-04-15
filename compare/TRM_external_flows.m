function out = TRM_external_flows(t, N)
    % Külső be- és kiáramlások definiálása a hálózati TRM modellhez
    % Bemenetek:
    % - t: idő
    % - N: útszakaszok száma
    %
    % Kimenet:
    % - out első oszlop: beáramlások, második oszlop: kiáramlások
    % - out(i,1): külső beáramlási ráta az i-edik útszakaszra
    % - out(i,2): kiáramlási ráta az i-edik útszakaszról

    flow = zeros(N, 2);

    % Itt definiálhatók a különböző időfüggő áramlások
    % Példák kommentben:
    flow(1, 1) = 0.05 * (sin(t/50) + 1.5);  % Első útszakaszra változó beáramlás
    %flow(2, 1) = morning_peak + evening_peak; % Második útszakaszra csúcsidőszaki beáramlás
    %flow(3, 2) = 0.03 * (cos(t/70) + 1.2);   % Harmadik útszakaszról változó kiáramlás

    out = flow;
end