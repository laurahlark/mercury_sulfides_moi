function core_params=fn_core_params(I, I_oss, rho_m, rho_S, rho_cr, R_S, R_m)
    R = 2439.36e3;
    M = 3.30111e23;
    G = 6.67408e-11;

    M_cr = 4/3*pi*(R^3-R_S^3) * rho_cr;
    M_S = 4/3*pi*(R_S^3-R_m^3) * rho_S;
    
    I_m = I_oss - 8*pi/15*rho_cr * (R^5-R_S^5) - 8*pi/15*rho_S * (R_S^5-R_m^5);
    R_c = (R_m^5 - 15/(8*pi)*I_m./rho_m)^(1/5);
    M_m = 4/3*pi*(R_m^3-R_c^3) * rho_m;

    M_silicate =  M_cr + M_S + M_m;
    M_c = M - M_silicate;

    I_c = I - I_oss;
    g_c = G/R_c^2 * M_c;
    
    P_S = G * rho_cr * ((M - 4*pi/3*rho_cr*R^3) * (1/R_S - 1/R) ...
        + 4*pi/3 * rho_cr * (R^2 - R_S^2)/2);
    P_m = P_S + G * rho_S * ((M - M_cr - 4*pi/3*rho_S*R_S^3) * (1/R_m - 1/R_S) ...
        + 4*pi/3 * rho_S * (R_S^2 - R_m^2)/2);
    P_c = P_m + G * rho_m * (...
        (M - M_cr - M_S - 4*pi/3*rho_m*R_m^3) * (1/R_c - 1/R_m) ...
        + 4*pi/3 * rho_m * (R_m^2 - R_c^2)/2);
    
    core_params.M_c = M_c;
    core_params.R_c = R_c;
    core_params.I_c = I_c;
    core_params.g_c = g_c;
    core_params.P_c = P_c;
end