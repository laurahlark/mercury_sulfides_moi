% Produces data for figures S3.

% How does the geodetically constrained core density and moment of inertia 
% compare to the geochemically predicted?
lib_amp_Bertone_2021 = 39.03;
dlib_amp_Bertone_2021 = 1.1;
[C_m, dC_m] = fn_get_C_m(lib_amp_Bertone_2021, dlib_amp_Bertone_2021);
C_m = (C_m + [-dC_m 0 dC_m]);

rho_sulfide = 2300;  % To explore other sulfide densities, change this.
mm_sulfide = 75;  % To explore other sulfide molar masses, change this.
rho_silicate = 3200;  % To explore other silicate densities, change this.

[rho_c_extracted, R_c_extracted, P_c_extracted, ...
    rho_c_dispersed, R_c_dispersed, P_c_dispersed, X_S] = ...
    fn_get_C_m_implications_by_S(rho_sulfide, mm_sulfide, rho_silicate, C_m);

M_c_extracted = rho_c_extracted .* 4/3*pi.*R_c_extracted.^3;
M_c_dispersed = rho_c_dispersed .* 4/3*pi.*R_c_dispersed.^3;

R = 2439.36e3;
M = 3.30111e23;

C_lo = .333;
geodetic_core_C_extracted_lo = (C_lo*M*R^2 - C_m)'./(M_c_extracted.*R_c_extracted.^2);
geodetic_core_C_dispersed_lo = (C_lo*M*R^2 - C_m)'./(M_c_dispersed.*R_c_dispersed.^2);

C_hi = .343;
geodetic_core_C_extracted_hi = (C_hi*M*R^2 - C_m)'./(M_c_extracted.*R_c_extracted.^2);
geodetic_core_C_dispersed_hi = (C_hi*M*R^2 - C_m)'./(M_c_dispersed.*R_c_dispersed.^2);
