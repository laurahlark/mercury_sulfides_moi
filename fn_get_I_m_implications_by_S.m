function [rho_c_extracted, R_c_extracted, P_c_extracted, ...
    rho_c_dispersed, R_c_dispersed, P_c_dispersed, m_frac_S] = ...
    fn_get_I_m_implications_by_S(rho_sulfide, mm_sulfide, rho_silicate, I_m)
    % What are the implications for M_c (rho_c), R_c of I_m?
    m_frac_S = 0:.01:.3;
    m_frac_sulfide = mm_sulfide/molar_mass("S")*m_frac_S;

    [M_c_extracted, R_c_extracted, P_c_extracted] = fn_get_I_m_implications(...
        rho_silicate, rho_sulfide, m_frac_sulfide, false, I_m);
    rho_c_extracted = M_c_extracted./(4/3*pi*R_c_extracted.^3);

    [M_c_dispersed, R_c_dispersed, P_c_dispersed] = fn_get_I_m_implications(...
        rho_silicate, rho_sulfide, m_frac_sulfide, true, I_m);
    rho_c_dispersed = M_c_dispersed./(4/3*pi*R_c_dispersed.^3);
end
