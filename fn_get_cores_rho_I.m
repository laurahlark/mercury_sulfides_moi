function [rho_avg_c, I_c]=fn_get_cores_rho_I(X_le, r_ic, T_of_r, EOS_generator)
    I_c = nan(length(X_le), length(r_ic));
    rho_avg_c = nan(length(X_le), length(r_ic));

    [EOS_params_l, EOS_l, EOS_params_s, EOS_s] = EOS_generator(X_le);
    default_R_c = 1990e3;
    default_P_c = 5.5e9;

    for i=1:length(X_le)
        for j=1:length(r_ic)
            [I_c(i, j), rho_avg_c(i, j)] = fn_get_core_profile(...
                default_R_c, r_ic(j), ...
                @(T, P) EOS_l(EOS_params_l(i), T, P), ...
                @(T, P) EOS_s(EOS_params_s(i), T, P), ...
                default_P_c, T_of_r);
            fprintf("X_le=%g,\t\t r_ic=%g,\t\t rho_c = %g\n", X_le(i), r_ic(j)/1e3, rho_avg_c(i, j));
        end
    end
end