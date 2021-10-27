function [I_c, rho_avg_c, r_in_core, P_in_core, rho_in_core] = fn_get_core_profile(...
    R_oc, R_ic, EOS_oc, EOS_ic, P_c, T_of_r)
    M = 3.30111e23;
    M_c_last = 0;
    M_c = .75 * M;
    while abs(M_c_last - M_c) > .0001 * M_c
        M_c_last = M_c;
        [r_in_core, P_in_core, rho_in_core] = eval_core(...
            R_oc, R_ic, EOS_oc, EOS_ic, P_c, T_of_r, M_c);
        M_c = -trapz(r_in_core, 4*pi*r_in_core.^2.*rho_in_core);
    end
    
    rho_avg_c = M_c/(4/3*pi*R_oc^3);
    I_c = -trapz(r_in_core, 8/3*pi*r_in_core.^4.*rho_in_core)/(M_c*R_oc^2);
end

function [r_in_core, P_in_core, rho_in_core] = eval_core(...
    R_oc, R_ic, EOS_oc, EOS_ic, P_c, T_of_r, M_c_guess) 
    r = flip(sqrt(linspace(0, R_oc^2, 1000)));
    [r_in_core, rho_in_core, P_in_core, M_interior_to_r_oc] = integrate_profile(...
        [r(r>R_ic) R_ic], P_c, M_c_guess, EOS_oc, T_of_r);
    if R_ic > 0
        [r_in_ic, rho_of_r_ic, P_of_r_ic, ~] = integrate_profile(...
            [R_ic r(r<R_ic)], P_in_core(r_in_core==R_ic), ...
            M_interior_to_r_oc(r_in_core==R_ic), EOS_ic, T_of_r);
        r_in_core = [r_in_core r_in_ic];
        P_in_core = [P_in_core P_of_r_ic];
        rho_in_core = [rho_in_core rho_of_r_ic];
    end
end


function [r, rho_of_r, P_of_r, M_interior_to_r]=integrate_profile(r, P_t, M_t, EOS, T_of_r)
    G = 6.67408e-11;
    P_of_r = nan(size(r));
    P_of_r(1) = P_t;
    M_interior_to_r = nan(size(r));
    M_interior_to_r(1) = M_t;
    rho_of_r = nan(size(r));
    for i=1:length(r)
        rho_of_r(i) = EOS(interp1(T_of_r.r, T_of_r.T, r(i), 'spline', 'extrap'), P_of_r(i));
        if i == length(r)
            break;
        end
        
        V_shell = 4/3*pi * (r(i)^3 - r(i+1)^3);  % R descending
        M_interior_to_r(i+1) = max(M_interior_to_r(i) - rho_of_r(i) * V_shell, 0);
        
        % Don't allow ridiculously high M resulting from overestimating the
        % mass of the core. Bound it at a reasonable max density (11K).
        M_for_g = min(M_interior_to_r(i), 4/3*pi*r(i)^3 * 11000);
        g = G * M_for_g / r(i)^2;
        P_of_r(i+1) = P_of_r(i) + rho_of_r(i)*g*(r(i) - r(i+1));
    end
end