% Generates data used for figure S7. How sensitive are the results to the
% choices of temperature profile, core radius, and CMB pressure?

T_of_r = fn_get_default_T_profile();

[EOS_params_l, EOS_l, EOS_params_s, EOS_s] = get_FeSi_Knibbe_2021(.07);
EOS_liq = @(T, P) EOS_l(EOS_params_l, T, P);
EOS_sol = @(T, P) EOS_s(EOS_params_s, T, P);

R_cmb = 1990e3;
P_cmb = 5.5e9;
T_cmb = interp1(T_of_r.r, T_of_r.T, R_cmb);

r_ic = 0;  % Change to inspect results for other inner core radii.

% Baseline.
[I, rho]  = fn_get_core_profile(...
    R_cmb, r_ic, EOS_liq, EOS_sol, P_cmb, T_of_r);

% Part a: high temp vs low temp.
T_lo_T_cmb = T_of_r;
T_lo_T_cmb.T = T_lo_T_cmb.T - 200;
[I_lo, rho_lo] = fn_get_core_profile(...
    R_cmb, r_ic, EOS_liq, EOS_sol, P_cmb, T_lo_T_cmb);
T_hi_T_cmb = T_of_r;
T_hi_T_cmb.T = T_hi_T_cmb.T + 200;
[I_hi, rho_hi] = fn_get_core_profile(...
    R_cmb, r_ic, EOS_liq, EOS_sol, P_cmb, T_hi_T_cmb);

fprintf("+- 200 K\t\t rho: %g %g %g, I: %g %g %g\n", ...
    rho_hi, rho, rho_lo, I_hi, I, I_lo);

% Part b: high vs low CMB heat flux.
T_lo_q = get_T_of_r(T_cmb, R_cmb, 3e-3);
[I_lo, rho_lo] = fn_get_core_profile(...
    R_cmb, r_ic, EOS_liq, EOS_sol, P_cmb, T_lo_q);
T_hi_q = get_T_of_r(T_cmb, R_cmb, 20e-3);
[I_hi, rho_hi] = fn_get_core_profile(...
    R_cmb, r_ic, EOS_liq, EOS_sol, P_cmb, T_hi_q);

fprintf("Low and high q\t\t rho: %g %g %g, I: %g %g %g\n", ...
    rho_lo, rho, rho_hi, I_lo, I, I_hi);

% Part c: extremes of R_c and P_c
R_cmb_lo = 1949e3;
P_cmb_hi = 5.9e9;
[I_lo, rho_lo] = fn_get_core_profile(...
    R_cmb_lo, r_ic, EOS_liq, EOS_sol, P_cmb_hi, T_of_r);
R_cmb_hi = 2020e3;
P_cmb_lo = 5.1e9;
[I_hi, rho_hi] = fn_get_core_profile(...
    R_cmb_hi, r_ic, EOS_liq, EOS_sol, P_cmb_lo, T_of_r);

fprintf("Extreme R_c and P_c\t\t rho: %g %g %g, I: %g %g %g\n", ...
    rho_lo, rho, rho_hi, I_lo, I, I_hi);


function T_of_r=get_T_of_r(T_cmb, R_cmb, q_cmb)
    % Following methods of Knibbe et al. (2021) 
    % https://doi.org/10.1029/2020JE006651
    F = -q_cmb;
    k = 41;
    r = flip(linspace(0, R_cmb));

    a = -.0308e-6;
    b = -.0154e-12;
    c = -.062e-18;
    d = .0045e-24;
    dTad_dr = T_cmb * (a + b*r + c*r.^2 + d*r.^3);
    r_cond = interp1(dTad_dr, r, F/k);
    dTdr = dTad_dr;
    dTdr(r > r_cond) = F/k;

    T_of_r.T = flip(T_cmb + cumtrapz(r, dTdr));
    T_of_r.r = flip(r);
end