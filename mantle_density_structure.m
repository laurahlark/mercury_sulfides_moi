% Generates the density profiles in Figure S1.

% Sulfide EOS.
EOS_params_MgS.molar_mass = molar_mass(["Mg" 1; "S" 1]);
EOS_params_MgS.T_0 = 298;
EOS_params_MgS.rho_0 = 2684;
EOS_params_MgS.K = 78.9e9;
EOS_params_MgS.K_d = 3.71;
EOS_params_MgS.alpha_0 = 6e-5;

EOS_params_CaS.molar_mass = molar_mass(["Ca" 1; "S" 1]);
EOS_params_CaS.T_0 = 298;
EOS_params_CaS.rho_0 = 2594;
EOS_params_CaS.K = 56e9;
EOS_params_CaS.K_d = 5.3;
EOS_params_CaS.alpha_0 = 4.8e-5;

% From Circone and Agee 1996
% https://doi.org/10.1016/0016-7037(96)00117-2
EOS_params_Fo.K = 128.54e9;
EOS_params_Fo.dKdT = -.02176;
EOS_params_Fo.K_d = 5.3;
EOS_params_Fo.alpha_0 = 3.1131e-5;
EOS_params_Fo.alpha_1 = 6.5693e-9;
EOS_params_Fo.alpha_2 = -5.8733e-1;
EOS_params_Fo.rho_0 = 3229;
EOS_params_Fo.T_0 = 298;

EOS_params_En.K = 100.9e9;
EOS_params_En.dKdT = -.00534;
EOS_params_En.K_d = 5;
EOS_params_En.alpha_0 = 4.847e-5;
EOS_params_En.alpha_1 = 0;
EOS_params_En.alpha_2 = 0;
EOS_params_En.rho_0 = 3206;
EOS_params_En.T_0 = 298;

P = 0:1e8:5.5e9;
[rho_cond_MgS, rho_conv_MgS, rcda_MgS, rcva_MgS] = get_rho_profile(...
    @(T, P) EOS_bm3_rho(EOS_params_MgS, T, P), P);
[rho_cond_CaS, rho_conv_CaS, rcda_CaS, rcva_CaS] = get_rho_profile(...
    @(T, P) EOS_bm3_rho(EOS_params_CaS, T, P), P);

[rho_cond_Fo, rho_conv_Fo, rcda_Fo, rcva_Fo] = get_rho_profile(...
    @(T, P) EOS_bm3_CA_rho(EOS_params_Fo, T, P), P);
[rho_cond_En, rho_conv_En, rcda_En, rcva_En] = get_rho_profile(...
    @(T, P) EOS_bm3_CA_rho(EOS_params_En, T, P), P);


%% Functions

function P=EOS_bm3(EOS_params, T, rho_norm_0)
    K = EOS_params.K;
    K_d = EOS_params.K_d;
    alpha_0 = EOS_params.alpha_0;
    dT = T - EOS_params.T_0;
    
    rho_norm = rho_norm_0 * (1 + alpha_0 * dT);
    P = 3/2*K * (rho_norm.^(7/3) - rho_norm.^(5/3)) .* ...
        (1 + 3/4*(K_d - 4) * (rho_norm.^(2/3) - 1));% + alpha_0 * K * dT;
end
function rho=EOS_bm3_rho(EOS_params, T, P) 
    rho_norm_all = .5:.001:1.2;
    P_all = EOS_bm3(EOS_params, T, rho_norm_all);
    rho = interp1(P_all, rho_norm_all, P) * EOS_params.rho_0;
end


function P=EOS_bm3_CA(EOS_params, T, rho_T_P_over_0)
    K = EOS_params.K + EOS_params.dKdT * (T - EOS_params.T_0);
    K_d = EOS_params.K_d;
    
    int_alpha = EOS_params.alpha_0 * (T - EOS_params.T_0) + ...
        EOS_params.alpha_1/2 * (T^2 - EOS_params.T_0^2) - ...
        EOS_params.alpha_2 * (1/T - 1/EOS_params.T_0);
    rho_T_P_over_T_0 = rho_T_P_over_0 * (1 + int_alpha);
    
    P = 3/2*K * (rho_T_P_over_T_0.^(7/3) - rho_T_P_over_T_0.^(5/3)) .* ...
        (1 + 3/4*(K_d - 4) * (rho_T_P_over_T_0.^(2/3) - 1));
end
function rho=EOS_bm3_CA_rho(EOS_params, T, P) 
    rho_norm_all = .5:.001:1.2;
    P_all = EOS_bm3_CA(EOS_params, T, rho_norm_all);
    rho = interp1(P_all, rho_norm_all, P) * EOS_params.rho_0;
end

function [rho_conductive, rho_convective, rho_cond_avg, rho_conv_avg]=get_rho_profile(EOS, P)
    R_c = 1990e3;
    P_c = 5.5e9;

    T_CMB = 1920;
    [r_conductive, ~, T_conductive] = fn_conductive_heat_loss_spherical(...
        440, T_CMB, [4 3 .2], [2.8 6.4 6.4]*1e-8, [R_c 2402e3 2435e3 2440e3]);
    [r_convective, ~, T_convective] = fn_conductive_heat_loss_spherical(...
        440, 1725, [4 3 .2], [2.8 6.4 6.4]*1e-8, [R_c+300e3 2402e3 2435e3 2440e3]);
    
    T_convective = [T_CMB 1726 T_convective];
    r_convective = [R_c R_c+299e3 r_convective];

    P_of_r = @(r) P_c * (2440e3 - r) / (2440e3 - R_c);
    rho_all_conductive = arrayfun(...
        @(i) EOS(T_conductive(i), P_of_r(r_conductive(i))), 1:length(T_conductive));
    rho_all_convective = arrayfun(...
        @(i) EOS(T_convective(i), P_of_r(r_convective(i))), 1:length(T_convective));
    
    rho_conductive = interp1(P_of_r(r_conductive), rho_all_conductive, P);
    rho_convective = interp1(P_of_r(r_convective), rho_all_convective, P);

    rho_cond_avg = V_mean(r_conductive, rho_all_conductive);
    rho_conv_avg = V_mean(r_convective, rho_all_convective);
end

function X_mean = V_mean(r, X)
    X_total = trapz(r, 4*pi*r.^2.*X);
    V = 4/3*pi*(max(r)^3 - min(r)^3);
    X_mean = X_total/V;
end
