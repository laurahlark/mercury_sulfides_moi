function [M_c, R_c, P_c, d_S, rho_mantle] = fn_get_I_m_implications(...
    rho_silicate, rho_sulfide, m_frac_sulfide, dispersed, I_m)
    if length(rho_silicate) == 1
        rho_silicate = rho_silicate * ones(size(m_frac_sulfide));
    end
    if length(rho_sulfide) == 1
        rho_sulfide = rho_sulfide * ones(size(m_frac_sulfide));
    end

    R = 2439.36e3;
    M = 3.30111e23;

    R_S = R - 38e3;  % Top of the sulfide layer.  38 km thick crust (PRMM).
    rho_cr = 2900;
    d_S_guess = m_frac_sulfide * 450e3;

    M_c = nan(length(I_m), length(m_frac_sulfide));
    R_c = nan(length(I_m), length(m_frac_sulfide));
    P_c = nan(length(I_m), length(m_frac_sulfide));
    rho_mantle = nan(length(I_m), length(m_frac_sulfide));
    d_S = nan(length(I_m), length(m_frac_sulfide));
    for i=1:length(I_m)
        for j=1:length(m_frac_sulfide)
            if dispersed
                R_m = R_S;
                X_S_mantle = m_frac_sulfide(j);  % Start with a guess as to volume fraction
                for k=1:5
                    rho_mantle(i, j) = (X_S_mantle * rho_sulfide(j) + ...
                        (1 - X_S_mantle) * rho_silicate(j));
                    core_params = fn_core_params(...
                        0, I_m(i), rho_mantle(i, j), rho_sulfide(j), rho_cr, ...
                        R_S, R_m);
                    if m_frac_sulfide(j)==0
                        break
                    end
                    M_sil = M - core_params.M_c;
                    V_m = 4/3*pi*(R_m^3 - core_params.R_c^3);
                    X_S_mantle = m_frac_sulfide(j) * M_sil / (rho_sulfide(j) * V_m);
                end
            else
                R_m = R_S-d_S_guess(j);
                for k=1:5
                    core_params = fn_core_params(...
                        0, I_m(i), rho_silicate(j), rho_sulfide(j), rho_cr, ...
                        R_S, R_m);
                    if m_frac_sulfide(j)==0
                        break
                    end
                    M_sil = M - core_params.M_c;
                    M_S = m_frac_sulfide(j)*M_sil;
                    V_S = M_S/rho_sulfide(j);
                    R_m = nthroot(R_S^3 - 3/(4*pi)*V_S, 3);
                end
                d_S(i, j) = R_S - R_m;
            end

            M_c(i, j) = core_params.M_c;
            R_c(i, j) = core_params.R_c;
            P_c(i, j) = core_params.P_c;
        end
    end
end