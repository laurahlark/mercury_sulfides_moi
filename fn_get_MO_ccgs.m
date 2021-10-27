function ccgs = fn_get_MO_ccgs(X_Si)
    % Calculates carbon content at graphite saturation (CCGS) in wt. frac
    % under magma ocean conditions (2200 K, 5 GPa) using equation of 
    % Steenstra et al. (2021) https://doi.org/10.1016/j.icarus.2019.113391
    a = -23.6;
    b = -1156;
    c = -36;
    d = 1.12;
    e = 11.39;
    T = 2200;
    P = 5e9;
    X_Ni = 0;
    ccgs = 10.^(...
        a + ...
        b/T + ...
        c*P/1e9/T + ...
        d*log10(100-X_Ni*100) + ...
        e*log10(100-X_Si*100))/100;  % /100 to get fraction from wt.%
end