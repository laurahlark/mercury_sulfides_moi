% Produces data for figure 1. What is the expected quantity of sulfide?
R = 2439.36e3;
M = 3.30111e23;
lib_amp_Bertone_2021 = 39.03 + [-1.1 0 1.1];
I_m = flip(fn_get_I_m(lib_amp_Bertone_2021)) * M*R^2;

rho_sulfide = 2600;
rho_silicate = 3200;
mm_sulfide = 60;

X_S = 0:.005:.17;

% Part a: simple mass fraction of sulfide.
X_sulfide = mm_sulfide/molar_mass("S") * X_S;

% Part b: sulfide layer thickness.
[~, R_c1, P_c1, d_sulfide, ~] = ...
    fn_get_I_m_implications(3200, 2600, X_sulfide, false, I_m);

% Part c: mantle density.
[~, R_c2, P_c2, ~, rho_m] = ...
    fn_get_I_m_implications(3200, 2600, X_sulfide, true, I_m);
