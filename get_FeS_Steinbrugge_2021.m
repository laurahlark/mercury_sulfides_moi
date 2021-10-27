function [EOS_params_l, EOS_l, EOS_params_s, EOS_s]=get_FeS_Steinbrugge_2021(X_S)
    mm = [55.845 87.9070] * 1e-3;

    Fe.molar_V = 6.88 * 1e-6;
    Fe.K = 148e9;
    Fe.K_d = 5.8;
    Fe.alpha_0 = 9e-5;
    Fe.T_0 = 298;
    Fe.ag_param = 5.1;
    Fe.kappa = .56;
    Fe.rho_0 = mm(1) / Fe.molar_V;

    FeS.molar_V = 22.97 * 1e-6;
    FeS.K = 17.02e9;
    FeS.K_d = 5.92;
    FeS.alpha_0 = 11.9e-5;
    FeS.T_0 = 1650;
    FeS.ag_param = 4.07505;
    FeS.kappa = 1.4;
    FeS.rho_0 = mm(2) / FeS.molar_V;

    EOS_params_l = [];
    for i=1:length(X_S)
        comp.mm = mm;
        comp.mol_frac = [1 - wt_S_to_mol_FeS(X_S(i)) wt_S_to_mol_FeS(X_S(i))];
        comp.X_le = X_S(i);
        comp.Fe_params = Fe;
        comp.FeS_params = FeS;

        % This is from Steinbrugge, but matches the values
        % and equations given by Knibbe.
        W11=-9.91275;
        W12=0.731385;
        W21=-1.32521;
        W22=1.72716;
        a = [W11 W21];
        b = [W12 W22];
        comp.fn_W = @(P) a + b * log(3/2 + P/1e9);
        EOS_params_l = [EOS_params_l, comp];
    end

    EOS_l = @EOS_liq;

    % All solids are pure iron.
    Fe.molar_V = 6.82 * 1e-6;
    Fe.K = 163.4e9;
    Fe.K_d = 5.38;
    Fe.alpha_0 = 7e-5;
    Fe.T_0 = 300;
    Fe.ag_param = 5.5;
    Fe.kappa = 1.4;
    Fe.rho_0 = mm(1) / Fe.molar_V;

    EOS_params_s = [];
    for i=1:length(X_S)
        EOS_params_s = [EOS_params_s, Fe];
    end

    EOS_s = @EOS_sol;
end

function P=EOS_Vinet(rho_norm, K, K_d)
    P = 3 * K .* rho_norm.^(2/3) .* (1-rho_norm.^(-1/3)) .* ...
        exp(3/2 * (K_d - 1) * (1 - rho_norm.^(-1/3)));
end

function rho=EOS_Komabayashi(EOS_params, T, P)
    rho_norm = .8:.001:2;
    P_all = EOS_Vinet(rho_norm, EOS_params.K, EOS_params.K_d);
    
    rho_P_0_over_rho_0 = interp1(P_all, rho_norm, P);
    alpha_P = EOS_params.alpha_0 * exp(...
        - EOS_params.ag_param/EOS_params.kappa * ...
        (1 - rho_P_0_over_rho_0.^-EOS_params.kappa));

    rho_P_T_over_rho_0 = rho_P_0_over_rho_0 ./ exp(alpha_P * (T - EOS_params.T_0));
    rho = rho_P_T_over_rho_0 * EOS_params.rho_0;
end

function V_ex=fn_V_excess(mol_frac, W)
    x2 = mol_frac(2);
    V_ex = x2 * (1 - x2) * (x2 * W(1) + (1-x2) * W(2)) * 1e-6;
end

function rho=EOS_liq(EOS_params, T, P)
    rho1 = EOS_Komabayashi(EOS_params.Fe_params, T, P);
    rho2 = EOS_Komabayashi(EOS_params.FeS_params, T, P);
    V1 = EOS_params.mm(1) ./ rho1;
    V2 = EOS_params.mm(2) ./ rho2;
    
    rho = nan(size(P));
    for i=1:length(P)
        V_ideal = sum([V1(i) V2(i)] .* EOS_params.mol_frac);
        V_excess = fn_V_excess(EOS_params.mol_frac, EOS_params.fn_W(P(i)));
        V_mix = V_ideal + V_excess;
        rho(i) = sum(EOS_params.mm .* EOS_params.mol_frac) / V_mix;
    end
end

function rho=EOS_sol(EOS_params, T, P)
    % Simply pure iron.
    rho = EOS_Komabayashi(EOS_params, T, P);
end

function mol_x_FeS = wt_S_to_mol_FeS(X_S)
    mol_x_S = wt_S_to_mol_S(X_S);
    mol_x_FeS = mol_x_S / (1 - mol_x_S);
end

function mol_x_S = wt_S_to_mol_S(X_S)
    mm_Fe = 55.845;
    mm_S = 32.06;
    mol_x_S = X_S * mm_Fe / (mm_S * (1-X_S) + mm_Fe * X_S);
end