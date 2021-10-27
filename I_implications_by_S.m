% Produces data for figures 3 and 4. 

% How does the geodetically constrained core density and moment of inertia 
% compare to the geochemically predicted?

% Figure 3, part a: geodetic constraints.

R = 2439.36e3;
M = 3.30111e23;
lib_amp_Bertone_2021 = 39.03 + [-1.1 0 1.1];
I_m = flip(fn_get_I_m(lib_amp_Bertone_2021)) * M*R^2;

rho_sulfide = 2600;  % To explore other sulfide densities, change this.
mm_sulfide = 60;  % To explore other sulfide molar masses, change this.
rho_silicate = 3200;  % To explore other silicate densities, change this.

[rho_c_extracted, R_c_extracted, rho_c_dispersed, R_c_dispersed, X_S] = ...
    fn_get_I_m_implications_by_S(rho_sulfide, mm_sulfide, rho_silicate, I_m);

M_c_extracted = rho_c_extracted .* 4/3*pi.*R_c_extracted.^3;
M_c_dispersed = rho_c_dispersed .* 4/3*pi.*R_c_dispersed.^3;


I_lo = .333;
geodetic_core_I_extracted_lo = (I_lo*M*R^2 - I_m)'./(M_c_extracted.*R_c_extracted.^2);
geodetic_core_I_dispersed_lo = (I_lo*M*R^2 - I_m)'./(M_c_dispersed.*R_c_dispersed.^2);

I_hi = .343;
geodetic_core_I_extracted_hi = (I_hi*M*R^2 - I_m)'./(M_c_extracted.*R_c_extracted.^2);
geodetic_core_I_dispersed_hi = (I_hi*M*R^2 - I_m)'./(M_c_dispersed.*R_c_dispersed.^2);

% Figure 3, part b: geochemical predictions for FeSi.

T_of_r = fn_get_default_T_profile();
r_ic__ = [0 600e3 1000e3 1200e3];
X_Si__ = 0:.02:.2;

[geochemical_core_rho, geochemical_core_I] = fn_get_cores_rho_I(...
    X_Si__, r_ic__, T_of_r, @get_FeSi_Knibbe_2021);


% Figure 4.

% What are the geodetic constraints on core silicon as a function of mantle
% sulfur, and what are the implications for planetary moment of inertia?

% To make these calculations for alternate alloys, use values for those
% instead of geochemical_core_rho, geochemical_core_I, and X_Si__.

X_Si_of_X_S_extracted = nan(length(r_ic__), length(I_m), length(X_S));
X_Si_of_X_S_dispersed = nan(length(r_ic__), length(I_m), length(X_S));
I_of_X_S_extracted = nan(length(r_ic__), length(I_m), length(X_S));
I_of_X_S_dispersed = nan(length(r_ic__), length(I_m), length(X_S));
for i=1:length(r_ic__)
    for j=1:length(I_m)
        fprintf("%g %g\n", i, j);
        X_Si_of_X_S_extracted(i, j, :) = interp1(...
            geochemical_core_rho(:, i), X_Si__, rho_c_extracted(j, :));
        core_I_of_X_S_extracted = squeeze(interp1(...
            X_Si__, geochemical_core_I(:, i), ...
            X_Si_of_X_S_extracted(i, j, :), 'linear', 'extrap'))';
        I_of_X_S_extracted(i, j, :) = core_I_of_X_S_extracted .* ...
            (M_c_extracted(j, :) .* R_c_extracted(j, :).^2) / (M*R^2) + ...
            I_m(j)/(M*R^2);
        
        X_Si_of_X_S_dispersed(i, j, :) = interp1(...
            geochemical_core_rho(:, i), X_Si__, rho_c_dispersed(j, :));
        core_I_of_X_S_dispersed = squeeze(interp1(...
            X_Si__, geochemical_core_I(:, i), ...
            X_Si_of_X_S_dispersed(i, j, :), 'linear', 'extrap'))';
        I_of_X_S_dispersed(i, j, :) = core_I_of_X_S_dispersed .* ...
            (M_c_dispersed(j, :) .* R_c_dispersed(j, :).^2) / (M*R^2) + ...
            I_m(j)/(M*R^2);
    end
end

