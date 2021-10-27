% Generates data for alternate core alloys:
%       - non-ideal FeSi
%       - FeS
%       - C-saturated FeSi

T_of_r = fn_get_default_T_profile();
r_ic__ = [0 600e3 1000e3 1200e3];
X_Si__ = 0:.02:.2;
X_S__ = 0:.005:.2;

[rho_avg_c_Si_nonideal, I_c_Si_nonideal] = fn_get_cores_rho_I(...
    X_Si__, r_ic__, T_of_r, @get_nonideal_FeSi_EOS);

[rho_avg_c_S, I_c_S] = fn_get_cores_rho_I(...
    X_S__, r_ic__, T_of_r, @get_FeS_Steinbrugge_2021);

[rho_avg_c_SiCsat, I_c_SiCsat] = fn_get_cores_rho_I(X_Si__, r_ic__, T_of_r, ...
    @(X_Si) get_FeSiC_Knibbe_2021(fn_get_MO_ccgs(X_Si), X_Si));


function [EOS_params_l, EOS_l, EOS_params_s, EOS_s] = get_nonideal_FeSi_EOS(X_Si)
    [EOS_params_l, EOS_l, ~, ~] = get_FeSi_Steinbrugge_2021(X_Si);
    [~, ~, EOS_params_s, EOS_s] = get_FeSi_Knibbe_2021(X_Si);
end