% Produces data for figure 2. What are the implications of sulfide for the
% density and radius of the core?
lib_amp_Bertone_2021 = 39.03;
dlib_amp_Bertone_2021 = 1.1;
[C_m, dC_m] = fn_get_C_m(lib_amp_Bertone_2021, dlib_amp_Bertone_2021);
C_m = (C_m + [-dC_m 0 dC_m]);

rho_sulfide = 2600;
rho_silicate = 3200;
mm_sulfide = 60;

[rho_c_extracted, R_c_extracted, P_c_extracted, ...
    rho_c_dispersed, R_c_dispersed, P_c_dispersed, X_S] = ...
    fn_get_C_m_implications_by_S(rho_sulfide, mm_sulfide, rho_silicate, C_m);
