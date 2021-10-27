% Produces data for figure 2. What are the implications of sulfide for the
% density and radius of the core?
R = 2439.36e3;
M = 3.30111e23;
lib_amp_Bertone_2021 = 39.03 + [-1.1 0 1.1];
I_m = flip(fn_get_I_m(lib_amp_Bertone_2021)) * M*R^2;

rho_sulfide = 2600;
rho_silicate = 3200;
mm_sulfide = 60;

[rho_c_extracted, R_c_extracted, rho_c_dispersed, R_c_dispersed, X_S] = ...
    fn_get_I_m_implications_by_S(rho_sulfide, mm_sulfide, rho_silicate, I_m);
