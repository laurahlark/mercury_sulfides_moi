function [EOS_params_l, EOS_l, EOS_params_s, EOS_s]=get_FeSiC_Knibbe_2021(X_C, X_Si)
    % Fe FeSi Fe3C
    mm = [55.845 83.9325 179.55] * 1e-3;

    Fe.molar_V = 6.88 * 1e-6;
    Fe.K = 148e9;
    Fe.K_d = 5.8;
    Fe.alpha_0 = 9e-5;
    Fe.T_0 = 298;
    Fe.ag_param = 5.1;
    Fe.kappa = .56;
    Fe.rho_0 = mm(1) / Fe.molar_V;

    FeSi.molar_V = 16.145 * 1e-6;
    FeSi.K = 72.586e9;
    FeSi.K_d = 7.555;
    FeSi.alpha_0 = 10.767e-5;
    FeSi.ag_param = 4.9857;
    FeSi.kappa = .56;
    FeSi.T_0 = 1723;
    FeSi.rho_0 = mm(2) / FeSi.molar_V;

    Fe3C.molar_V = 26.68 * 1e-6;
    Fe3C.K = 75.66e9;
    Fe3C.K_d = 7.98;
    Fe3C.alpha_0 = 9.59e-5;
    Fe3C.T_0 = 1723;
    Fe3C.ag_param = 9.43;
    Fe3C.kappa = .56;
    Fe3C.rho_0 = mm(3) / Fe3C.molar_V;
    
    EOS_params_l = [];
    for i=1:length(X_Si)
        comp.mm = mm;
        comp.X_le = X_Si(i);
        comp.mol_frac = wt_Si_C_to_mol_frac(X_Si(i), X_C(i));
        comp.params_Fe = Fe;
        comp.params_FeSi = FeSi;
        comp.params_Fe3C = Fe3C;
        EOS_params_l = [EOS_params_l, comp];
    end

    EOS_l = @EOS_liq;
    
    % Now do for the solid. Knibbe says:
    % ?the density profile of the inner core (if present) is characterized 
    % by an Fe-Si-C ideal solid mixing model with solid fcc Fe 
    % (Komabayashi et al., 2014), Fe3C (Litasov et al., 2013), and 
    % ,DO3 Fe-16 wt% Si (Fischer et al., 2012) endmembers
    
    % These are fcc Fe, DO3 Fe-16wt.%Si
    mm = [55.845 55.845*(1-.2747)+28.0855*.2747] * 1e-3; %  179.55 (Fe3C)

    Fe.molar_V = 6.82 * 1e-6;
    Fe.K = 163.4e9;
    Fe.K_d = 5.38;
    Fe.alpha_0 = 7e-5;
    Fe.T_0 = 300;
    Fe.ag_param = 5.5;
    Fe.kappa = 1.4;
    Fe.rho_0 = mm(1) / Fe.molar_V;

    Fe16Si.molar_V = 6.799e-6;
    Fe16Si.K = 193.4e9;
    Fe16Si.K_d = 4.91;
    Fe16Si.g_param_0 = 1.89;
    Fe16Si.debye_0 = 417;
    Fe16Si.kappa = 1;
    Fe16Si.rho_0 = mm(2) / Fe16Si.molar_V;

    EOS_params_s = [];
    for i=1:length(X_Si)
        comp.mm = mm;
        comp.mol_frac = [1-wt_Si_to_mol_Fe16Si(X_Si(i)) wt_Si_to_mol_Fe16Si(X_Si(i))];
        comp.params_Fe = Fe;
        comp.params_Fe16Si = Fe16Si;
        EOS_params_s = [EOS_params_s, comp];
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

function rho=EOS_liq(EOS_params, T, P)
    rho1 = EOS_Komabayashi(EOS_params.params_Fe, T, P);
    rho2 = EOS_Komabayashi(EOS_params.params_FeSi, T, P);
    rho3 = EOS_Komabayashi(EOS_params.params_Fe3C, T, P);
    V1 = EOS_params.mm(1) ./ rho1;
    V2 = EOS_params.mm(2) ./ rho2;
    V3 = EOS_params.mm(3) ./ rho3;
    partial_V1 = EOS_params.mol_frac(1)*V1;
    partial_V2 = EOS_params.mol_frac(2)*V2;
    partial_V3 = EOS_params.mol_frac(3)*V3;
    Vmix = partial_V1 + partial_V2 + partial_V3;
    rho = sum(EOS_params.mm .* EOS_params.mol_frac) ./ Vmix;
end

function P=EOS_bm3(EOS_params, rho_norm) 
    P = 3/2*EOS_params.K * (rho_norm.^(7/3) - rho_norm.^(5/3)) .* ...
        (1 + 3/4*(EOS_params.K_d - 4) * (rho_norm.^(2/3) - 1));
end

function rho=EOS_Fischer(EOS_params, T, P)
    rho_norm_all = .8:.01:1.8;
    P_300 = EOS_bm3(EOS_params, rho_norm_all);
    
    V_norm = 1 ./ rho_norm_all;
    g_param = EOS_params.g_param_0 * V_norm .^ EOS_params.kappa;
    debye = EOS_params.debye_0 * exp(...
        EOS_params.g_param_0 * (1 - V_norm.^EOS_params.kappa) / ...
        EOS_params.kappa);

    V = V_norm * EOS_params.molar_V;
    x = linspace(min(debye/T), max(debye/T), 5);
    debye_fn_300 = debye_fn(debye/300);
    debye_fn_x = debye_fn(x);
    E = (T .* interp1(x, debye_fn_x, debye/T) - 300 * debye_fn_300);
    P_th = g_param ./ V .* E;
    
    rho = interp1(P_300  + P_th, rho_norm_all, P) * EOS_params.rho_0;
end

function rho=EOS_sol(EOS_params, T, P)
    rho1 = EOS_Komabayashi(EOS_params.params_Fe, T, P);
    rho2 = EOS_Fischer(EOS_params.params_Fe16Si, T, P);

    V1 = EOS_params.mm(1) ./ rho1;
    V2 = EOS_params.mm(2) ./ rho2;
    partial_V1 = EOS_params.mol_frac(1)*V1;
    partial_V2 = EOS_params.mol_frac(2)*V2;
    Vmix = partial_V1 + partial_V2;
    rho = sum(EOS_params.mm .* EOS_params.mol_frac) ./ Vmix;
end

% Molar E as fn of T.
function E=debye_fn(x_outer)
    R = 8.31432;
    n = 1; %6.0221409e23;
    persistent FeSiCS_x__;
    persistent FeSiCS_d3__;
    eps = .0000001;
    if isempty(FeSiCS_x__)
        FeSiCS_x__ = linspace(eps, 4, 400);
        FeSiCS_d3__ = FeSiCS_x__.^3 .* log(1 - exp(-FeSiCS_x__)) ...
            - 3 * FeSiCS_x__.^2 .* polylog(2, exp(-FeSiCS_x__)) ...
            - 6 * FeSiCS_x__ .* polylog(3, exp(-FeSiCS_x__)) ...
            - 6 * polylog(4, exp(-FeSiCS_x__));
    end

    d3 = @(x) interp1(FeSiCS_x__, FeSiCS_d3__, x);
    E = 9 * R * n * x_outer.^-3 .* (d3(x_outer) - d3(eps));
end

function mol_x_Fe16Si=wt_Si_to_mol_Fe16Si(X_Si)
    mm_Fe = 55.845;
    mm_Si = 28.0855;
    X_Si_in_Fe16Si = .2747;  % From Fischer 2016
    
    mol_x_Si = X_Si * mm_Fe / (mm_Si * (1-X_Si) + mm_Fe * X_Si);
    mol_x_Fe16Si = mol_x_Si / X_Si_in_Fe16Si;
end

% Gives mol frac of Fe, FeSi, and Fe3C.
function mol_frac=wt_Si_C_to_mol_frac(X_Si, X_C)
    mm_Fe = 55.845;
    mm_Si = 28.0855;
    mm_C = 12.011;
    
    mol_x_Si = X_Si * mm_Fe / (mm_Si * (1-X_Si) + mm_Fe * X_Si);
    mol_x_C = X_C * mm_Fe / (mm_C * (1-X_C) + mm_Fe * X_C);

    mol_x_FeSi = mol_x_Si / (1 - mol_x_Si - 3*mol_x_C);
    mol_x_Fe3C = mol_x_C / (1 - mol_x_Si - 3*mol_x_C);
    mol_x_Fe = 1 - mol_x_FeSi - mol_x_Fe3C;
    mol_frac = [mol_x_Fe mol_x_FeSi mol_x_Fe3C];
end
