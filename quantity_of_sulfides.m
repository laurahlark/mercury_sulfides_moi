% Produces data for figure 1. What is the expected quantity of sulfide?
lib_amp_Bertone_2021 = 39.03;
dlib_amp_Bertone_2021 = 1.1;
[C_m, dC_m] = fn_get_C_m(lib_amp_Bertone_2021, dlib_amp_Bertone_2021);
C_m = (C_m + [-dC_m 0 dC_m]);

rho_sulfide = 2600;
rho_silicate = 3200;
mm_sulfide = 60;

X_S = 0:.005:.17;

% Part a: simple mass fraction of sulfide.
X_sulfide = mm_sulfide/molar_mass("S") * X_S;

% Part b: sulfide layer thickness.
[~, R_c1, P_c1, d_sulfide, ~] = ...
    fn_get_C_m_implications(3200, 2600, X_sulfide, false, C_m);

% Part c: mantle density.
[~, R_c2, P_c2, ~, rho_m] = ...
    fn_get_C_m_implications(3200, 2600, X_sulfide, true, C_m);
