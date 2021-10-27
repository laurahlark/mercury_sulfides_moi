function I_m = fn_get_I_m(libration_amplitude)
    to_radians = 2*pi/(360*60*60);
    f_of_e_c = @(e_c) 1 - 11*e_c^2 + 959/48 * e_c^4;
    eccentricity = 0.20563;  % Baland 2017
    C_22 = .80389e-5;  % pm .0006 (uncertainty dwarfed by , from Margot 2018
    I_m = 6 * f_of_e_c(eccentricity) * C_22 ./ (libration_amplitude*to_radians);
end