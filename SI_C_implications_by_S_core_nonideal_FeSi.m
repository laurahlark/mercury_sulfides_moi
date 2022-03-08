% Produces data for figure S4

% How does the geodetically constrained core density and moment of inertia 
% compare to the geochemically predicted?

% Part a

lib_amp_Bertone_2021 = 39.03;
dlib_amp_Bertone_2021 = 1.1;
[C_m, dC_m] = fn_get_C_m(lib_amp_Bertone_2021, dlib_amp_Bertone_2021);
C_m = (C_m + [-dC_m 0 dC_m]);

rho_sulfide = 2600;  % To explore other sulfide densities, change this.
mm_sulfide = 60;  % To explore other sulfide molar masses, change this.
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

T_of_r = fn_get_default_T_profile();
r_ic__ = [0 1200e3];
X_Si__ = 0:.01:.2;

default_R_c = 1990e3;
default_P_c = 5.5e9;
[geochemical_core_rho, geochemical_core_C] = fn_get_cores_rho_C(...
    X_Si__, r_ic__, T_of_r, default_R_c, default_P_c, ...
    @(X_Si) get_nonideal_FeSi_EOS(X_Si));


% Part b

% What are the geodetic constraints on core silicon as a function of mantle
% sulfur, and what are the implications for planetary moment of inertia?

% To make these calculations for alternate alloys, use values for those
% instead of geochemical_core_rho, geochemical_core_C, and X_Si__.
X_Si_of_X_S_extracted = nan(length(r_ic__), length(C_m), length(X_S));
X_Si_of_X_S_dispersed = nan(length(r_ic__), length(C_m), length(X_S));
C_of_X_S_extracted = nan(length(r_ic__), length(C_m), length(X_S));
C_of_X_S_dispersed = nan(length(r_ic__), length(C_m), length(X_S));
error_ex = nan(length(r_ic__), length(C_m), length(X_S));
error_d = nan(length(r_ic__), length(C_m), length(X_S));
for i=1:length(r_ic__)
    for j=1:length(C_m)
        for k=1:length(X_S)
            fprintf("%g %g %g\n", i, j, k);
            % Two cases.
            
            % 1. Extracted.
            preliminary_core_Si_ex = interp1(...
                geochemical_core_rho(:, i), X_Si__, rho_c_extracted(j, k), 'linear', 'extrap');
            X_Si_ex = preliminary_core_Si_ex-.005:.005:preliminary_core_Si_ex+.005;
            [final_rho_ex, final_C_ex] = fn_get_cores_rho_C(...
                X_Si_ex, r_ic__(i), T_of_r, ...
                R_c_extracted(j, k), P_c_extracted(j, k), ...
                @(X_Si) get_nonideal_FeSi_EOS(X_Si));
            final_core_Si_ex = interp1(...
                final_rho_ex, X_Si_ex, rho_c_extracted(j, k));
            error_ex(i, j, k) = final_core_Si_ex - preliminary_core_Si_ex;
            fprintf("Extracted case: default R/P: %g, exact R/P: %g\n", preliminary_core_Si_ex, final_core_Si_ex);
            
            X_Si_of_X_S_extracted(i, j, k) = final_core_Si_ex;
            core_C_of_X_S_extracted = interp1(X_Si_ex, final_C_ex, final_core_Si_ex);
            C_of_X_S_extracted(i, j, k) = core_C_of_X_S_extracted .* ...
                (M_c_extracted(j, k) .* R_c_extracted(j, k).^2) / (M*R^2) + ...
                C_m(j)/(M*R^2);
            
            % 2. Dispersed.
            preliminary_core_Si_d = interp1(...
                geochemical_core_rho(:, i), X_Si__, rho_c_dispersed(j, k), 'linear', 'extrap');
            X_Si_d = preliminary_core_Si_d-.005:.005:preliminary_core_Si_d+.005;
            [final_rho_d, final_C_d] = fn_get_cores_rho_C(...
                X_Si_d, r_ic__(i), T_of_r, ...
                R_c_dispersed(j, k), P_c_dispersed(j, k), ...
                @(X_Si) get_nonideal_FeSi_EOS(X_Si));
            final_core_Si_d = interp1(final_rho_d, X_Si_d, rho_c_dispersed(j, k));
            error_d(i, j, k) = final_core_Si_d - preliminary_core_Si_d;
            fprintf("Dispersed case: default R/P: %g, exact R/P: %g\n", preliminary_core_Si_d, final_core_Si_d);
            
            X_Si_of_X_S_dispersed(i, j, k) = final_core_Si_d;
            core_C_of_X_S_dispersed = interp1(X_Si_d, final_C_d, final_core_Si_d);
            C_of_X_S_dispersed(i, j, k) = core_C_of_X_S_dispersed .* ...
                (M_c_dispersed(j, k) .* R_c_dispersed(j, k).^2) / (M*R^2) + ...
                C_m(j)/(M*R^2);
        end
    end
end

function [EOS_params_l, EOS_l, EOS_params_s, EOS_s] = get_nonideal_FeSi_EOS(X_Si)
    [EOS_params_l, EOS_l, ~, ~] = get_FeSi_Steinbrugge_2021(X_Si);
    [~, ~, EOS_params_s, EOS_s] = get_FeSiC_Knibbe_2021(zeros(size(X_Si)), X_Si);
end