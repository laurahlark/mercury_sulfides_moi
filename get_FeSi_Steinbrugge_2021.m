function [EOS_params_l, EOS_l, EOS_params_s, EOS_s]=get_FeSi_Steinbrugge_2021(X_Si)
    mm = [55.845 83.9250] * 1e-3;

    Fe.molar_V = 6.88 * 1e-6;
    Fe.K = 148e9;
    Fe.K_d = 5.8;
    Fe.alpha_0 = 9e-5;
    Fe.T_0 = 298;
    Fe.ag_param = 5.1;
    Fe.kappa = .56;
    Fe.rho_0 = mm(1) / Fe.molar_V;

    FeSi.molar_V = 16.5839 * 1e-6;
    FeSi.K = 69.0074e9;
    FeSi.K_d = 7.76;
    FeSi.alpha_0 = 17.6525e-5;
    FeSi.T_0 = 1723;
    FeSi.ag_param = 4.07505;
    FeSi.kappa = .56;
    FeSi.rho_0 = mm(2) / FeSi.molar_V;

    EOS_params_l = [];
    for i=1:length(X_Si)
        comp.mm = mm;
        comp.X_le = X_Si(i);
        comp.params1 = Fe;
        comp.params2 = FeSi;

        % From https://github.com/gregorsteinbruegge/MercuryInterior
        W1=-2.3199284685783192;
        W2=-1.2489264297620897;
        a = [W1 W2];
        b = [0 0];
        comp.fn_W = @(P) a + b * log(3/2 + P);
        EOS_params_l = [EOS_params_l, comp];
    end

    EOS_l = @EOS_liq;

    Fe.molar_V = 6.82 * 1e-6;
    Fe.K = 163.4e9;
    Fe.K_d = 5.38;
    Fe.alpha_0 = 7e-5;
    Fe.T_0 = 300;
    Fe.ag_param = 5.5;
    Fe.kappa = 1.4;
    Fe.rho_0 = mm(1) / Fe.molar_V;

    EOS_params_s = [];
    for i=1:length(X_Si)
        comp.Fe_params = Fe;
        EOS_params_s = [EOS_params_s, comp];
    end

    EOS_s = @EOS_sol;
end

function P=EOS_Vinet(rho_norm, K, K_d)
    P = 3 * K .* rho_norm.^(2/3) .* (1-rho_norm.^(-1/3)) .* ...
            exp(3/2 * (K_d - 1) * (1 - rho_norm.^(-1/3)));
end

function rho=EOS_Komabayashi(EOS_params, T, P)
    rho_norm = .8:.001:1.8;
    P_all = EOS_Vinet(rho_norm, EOS_params.K, EOS_params.K_d);
    
    rho_P_0_over_rho_0 = interp1(P_all, rho_norm, P);
    alpha_P = EOS_params.alpha_0 * exp(...
        - EOS_params.ag_param/EOS_params.kappa * ...
        (1 - rho_P_0_over_rho_0.^-EOS_params.kappa));

    rho_P_T_over_rho_0 = rho_P_0_over_rho_0 ./ exp(alpha_P * (T - EOS_params.T_0));
    rho = rho_P_T_over_rho_0 * EOS_params.rho_0;
end

function V_ex=V_excess(mol_frac, W)
    x2 = mol_frac(2);
    V_ex = x2 * (1 - x2) * (x2 * W(1) + (1-x2) * W(2)) * 1e-6;
end

function rho=EOS_liq(EOS_params, T, P)
    rho1 = EOS_Komabayashi(EOS_params.params1, T, P);
    rho2 = EOS_Komabayashi(EOS_params.params2, T, P);
    V1 = EOS_params.mm(1) ./ rho1;
    V2 = EOS_params.mm(2) ./ rho2;
    
    X_Si = EOS_params.X_le;
    mol_X_Si = wt_Si_to_mol_Si(X_Si);
    chi(2) = mol_X_Si / (1 - mol_X_Si);
    chi(1) = 1 - chi(2);
    
    rho = nan(size(P));
    for i=1:length(P)
        W = EOS_params.fn_W(P(i));
        V = [V1(i) V2(i)];
        V_ideal = sum(V .* chi);
        V_mix = V_ideal + V_excess(chi, W);
        rho(i) = sum(EOS_params.mm .* chi) / V_mix;
    end
end

function rho=EOS_sol(EOS_params, T, P)
    rho1 = EOS_Komabayashi(EOS_params.Fe_params, T, P);

    X_Si = EOS_params.X_le;
    mm_Fe = 55.845;
    mm_Si = 28.0855;
    
    mol_x_Si = wt_Si_to_mol_Si(X_Si);
    mol_frac = [1-mol_x_Si mol_x_Si];
    mm = [mm_Fe mm_Si];
    rho = rho1 * sum(mol_frac .* mm) / mm_Fe;
end

function mol_x_Si = wt_Si_to_mol_Si(X_Si)
    mm_Fe = 55.845;
    mm_Si = 28.0855;
    mol_x_Si = X_Si * mm_Fe / (mm_Si * (1-X_Si) + mm_Fe * X_Si);
end