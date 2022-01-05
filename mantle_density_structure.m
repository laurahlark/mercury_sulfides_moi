% Generates data (and plots) for mantle density structure and figure S1.
addpath('figures/')

% Sulfide EOS.
EOS_params_MgS.molar_mass = molar_mass(["Mg" 1; "S" 1]);
EOS_params_MgS.T_0 = 298;
EOS_params_MgS.rho_0 = 2661;  % 4 MgS/cell molar_mass / 140.68 A^3/cell
EOS_params_MgS.K = 78.9e9;
EOS_params_MgS.K_d = 3.71;
EOS_params_MgS.alpha_0 = 6e-5;
EOS_MgS = @(T, P) EOS_bm3_rho(EOS_params_MgS, T, P);

EOS_params_CaS.molar_mass = molar_mass(["Ca" 1; "S" 1]);
EOS_params_CaS.T_0 = 298;
EOS_params_CaS.rho_0 = 2607;  % 4 CaS/cell molar_mass / Unit cell from work used by Peiris.
EOS_params_CaS.K = 56e9;
EOS_params_CaS.K_d = 5.3;
EOS_params_CaS.alpha_0 = 4.8e-5;
EOS_CaS = @(T, P) EOS_bm3_rho(EOS_params_CaS, T, P);

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
EOS_Fo = @(T, P) EOS_bm3_CA_rho(EOS_params_Fo, T, P);

EOS_params_En.K = 100.9e9;
EOS_params_En.dKdT = -.00534;
EOS_params_En.K_d = 5;
EOS_params_En.alpha_0 = 4.847e-5;
EOS_params_En.alpha_1 = 0;
EOS_params_En.alpha_2 = 0;
EOS_params_En.rho_0 = 3206;
EOS_params_En.T_0 = 298;
EOS_En = @(T, P) EOS_bm3_CA_rho(EOS_params_En, T, P);

[r_conductive, ~, T_conductive] = fn_conductive_heat_loss_spherical(...
    440, 1920, [4 3 .2], [0 0 0]*1e-8, [1990e3 2402e3 2435e3 2440e3]);
min_T_crmb_unheated = T_conductive(r_conductive==2402e3);
[r_conductive, ~, T_conductive] = fn_conductive_heat_loss_spherical(...
    440, 1920, [4 3 .2], [0 6.8 6.8]*1e-8, [1990e3 2402e3 2435e3 2440e3]);
min_T_crmb_heated = T_conductive(r_conductive==2402e3);
[r_conductive, ~, T_conductive] = fn_conductive_heat_loss_spherical(...
    440, 1920, [4 3 .2], [.3*6.8 6.8 6.8]*1e-8, [1990e3 2402e3 2435e3 2440e3]);
min_T_crmb_heated_all = T_conductive(r_conductive==2402e3);
fprintf("CrMB temperature is %g (unheated) or %g (heated crust) or %g (heated mantle) " + ...
    "for conductive profile.\n", ...
    min_T_crmb_unheated, min_T_crmb_heated, min_T_crmb_heated_all);

P = [.5e9 5.5e9];
T = [800 1400 1900 2100];
rho_MgS = fn_get_rho_profiles_P_T(EOS_MgS, P, T);
rho_CaS = fn_get_rho_profiles_P_T(EOS_CaS, P, T);
rho_En = fn_get_rho_profiles_P_T(EOS_En, P, T);
rho_Fo = fn_get_rho_profiles_P_T(EOS_Fo, P, T);

% What are the effects of temperature?
fprintf("Effects of temperature over 1300 (500) K: \n");
drho_dT_MgS = (rho_MgS(:, 1)-rho_MgS(:, end))/rho_MgS(1, 1);
drho_dT_MgS2 = (rho_MgS(:, T==1400)-rho_MgS(:, T==1900))/rho_MgS(1, T==1400);
fprintf("\tMgS: %g-%g (%g-%g) percent\n", drho_dT_MgS*100, drho_dT_MgS2*100);
drho_dT_CaS = (rho_CaS(:, 1)-rho_CaS(:, end))/rho_CaS(1, 1);
drho_dT_CaS2 = (rho_CaS(:, T==1400)-rho_CaS(:, T==1900))/rho_CaS(1, T==1400);
fprintf("\tCaS: %g-%g (%g-%g) percent\n", drho_dT_CaS*100, drho_dT_CaS2*100);
drho_dT_En = (rho_En(:, 1)-rho_En(:, end))/rho_En(1, 1);
drho_dT_En2 = (rho_En(:, T==1400)-rho_En(:, T==1900))/rho_En(1, T==1400);
fprintf("\tEn: %g-%g (%g-%g) percent\n", drho_dT_En*100, drho_dT_En2*100);
drho_dT_Fo = (rho_Fo(:, 1)-rho_Fo(:, end))/rho_Fo(1, 1);
drho_dT_Fo2 = (rho_Fo(:, T==1400)-rho_Fo(:, T==1900))/rho_Fo(1, T==1400);
fprintf("\tFo: %g-%g (%g-%g) percent\n", drho_dT_Fo*100, drho_dT_Fo2*100);

fprintf("Effects of pressure over 5 GPa at 800-2100 K: \n");
drho_dP_MgS = (rho_MgS(end, [1 end])-rho_MgS(1, [1 end]))/rho_MgS(1, 1);
fprintf("\tMgS: %g-%g percent\n", drho_dP_MgS*100);
drho_dP_CaS = (rho_CaS(end, [1 end])-rho_CaS(1, [1 end]))/rho_CaS(1, 1);
fprintf("\tCaS: %g-%g percent\n", drho_dP_CaS*100);
drho_dP_En = (rho_En(end, [1 end])-rho_En(1, [1 end]))/rho_En(1, 1);
fprintf("\tEn: %g-%g percent\n", drho_dP_En*100);
drho_dP_Fo = (rho_Fo(end, [1 end])-rho_Fo(1, [1 end]))/rho_Fo(1, 1);
fprintf("\tFo: %g-%g percent\n", drho_dP_Fo*100);

max_drho_comp = max(max((rho_Fo-rho_CaS)./rho_CaS))*100;
min_drho_comp = min(min((rho_En-rho_MgS)./rho_MgS))*100;
fprintf("Compositional density change range %g-%g sulfide to silicate\n", ...
    min_drho_comp, max_drho_comp);
max_drho_comp = max(max((rho_Fo-rho_En)./rho_En))*100;
min_drho_comp = min(min((rho_Fo-rho_En)./rho_En))*100;
fprintf("\tand %g-%g silicate to silicate\n", min_drho_comp, max_drho_comp);
max_drho_comp = max(max((rho_MgS-rho_CaS)./rho_CaS))*100;
min_drho_comp = min(min((rho_MgS-rho_CaS)./rho_CaS))*100;
fprintf("\tand %g-%g sulfide to sulfide\n", min_drho_comp, max_drho_comp);

R_c = 1990e3;
P_c = 5.5e9;
P_of_r = @(r) P_c * (2440e3 - r) / (2440e3 - R_c);

T_hi = @(r) 1800;
T_lo = @(r) 1300;

r = flip(1990e3:1e3:2440e3);
[rho_cond_MgS, rho_sol_MgS] = fn_get_rho_profiles(...
    r, EOS_MgS, P_of_r, T_lo, T_hi);
[rho_cond_CaS, rho_sol_CaS] = fn_get_rho_profiles(...
    r, EOS_CaS, P_of_r, T_lo, T_hi);
[rho_cond_En, rho_sol_En] = fn_get_rho_profiles(...
    r, EOS_En, P_of_r, T_lo, T_hi);
[rho_cond_Fo, rho_sol_Fo] = fn_get_rho_profiles(...
    r, EOS_Fo, P_of_r, T_lo, T_hi);

fprintf("Volume-averaged densities for 1300 K (1800 K) profile:\n");
fprintf("\tMgS: %g (%g)\n", -V_mean(r, rho_cond_MgS), -V_mean(r, rho_sol_MgS));
fprintf("\tCaS: %g (%g)\n", -V_mean(r, rho_cond_CaS), -V_mean(r, rho_sol_CaS));
fprintf("\tEn: %g (%g)\n", -V_mean(r, rho_cond_En), -V_mean(r, rho_sol_En));
fprintf("\tFo: %g (%g)\n", -V_mean(r, rho_cond_Fo), -V_mean(r, rho_sol_Fo));

T_solidus = @(P) 1421 + 177 .* P ./ 1e9 - 12.2 * (P ./ 1e9).^2;
T_solidus_of_r = @(r) T_solidus(P_of_r(r));
[r_conductive, ~, T_conductive] = fn_conductive_heat_loss_spherical(...
    440, 1920, [4 3 .2], [0 6.8 6.8]*1e-8, [1990e3 2402e3 2435e3 2440e3]);
T_conductive_of_r = @(r) interp1(r_conductive, T_conductive, r);

r = flip(1990e3:1e3:2440e3);
[rho_cond_MgS, rho_sol_MgS] = fn_get_rho_profiles(...
    r, EOS_MgS, P_of_r, T_conductive_of_r, T_solidus_of_r);
[rho_cond_CaS, rho_sol_CaS] = fn_get_rho_profiles(...
    r, EOS_CaS, P_of_r, T_conductive_of_r, T_solidus_of_r);
[rho_cond_En, rho_sol_En] = fn_get_rho_profiles(...
    r, EOS_En, P_of_r, T_conductive_of_r, T_solidus_of_r);
[rho_cond_Fo, rho_sol_Fo] = fn_get_rho_profiles(...
    r, EOS_Fo, P_of_r, T_conductive_of_r, T_solidus_of_r);

fprintf("Volume-averaged densities for conductive (solidus) profile:\n");
fprintf("\tMgS: %g (%g)\n", -V_mean(r, rho_cond_MgS), -V_mean(r, rho_sol_MgS));
fprintf("\tCaS: %g (%g)\n", -V_mean(r, rho_cond_CaS), -V_mean(r, rho_sol_CaS));
fprintf("\tEn: %g (%g)\n", -V_mean(r, rho_cond_En), -V_mean(r, rho_sol_En));
fprintf("\tFo: %g (%g)\n", -V_mean(r, rho_cond_Fo), -V_mean(r, rho_sol_Fo));


colors = lines(4);

close all
figure(1);
set(gcf, 'position', [0 0 400 600]);
hold on
plot(T_conductive_of_r(r), r/1e3, 'k-', 'linewidth', 2, 'DisplayName', 'Conductive');
plot(T_solidus_of_r(r), r/1e3, 'k--', 'linewidth', 2, 'DisplayName', 'Solidus');

patch([0 1e10 1e10 0 0], [2402 2402 2440 2440 2402], 'k', 'EdgeAlpha', 0, 'FaceAlpha', .1);
hold off
set(gca, 'fontsize', 20);
ylabel("Radius (km)");
xlabel("Temperature (K)");
xlim([400 2100]);
ylim([1990 2440]);

saveas(gcf, 'figures/mantle_T.png');

figure(2);
set(gcf, 'position', [400 0 500 600]);
hold on
plot(rho_cond_MgS, r/1e3, '-', 'color', colors(1, :), 'linewidth', 2);
plot(rho_sol_MgS, r/1e3, '--', 'color', colors(1, :), 'linewidth', 2);
plot(rho_cond_CaS, r/1e3, '-', 'color', colors(2, :), 'linewidth', 2);
plot(rho_sol_CaS, r/1e3, '--', 'color', colors(2, :), 'linewidth', 2);
plot(rho_cond_Fo, r/1e3, '-', 'color', colors(3, :), 'linewidth', 2);
plot(rho_sol_Fo, r/1e3, '--', 'color', colors(3, :), 'linewidth', 2);
plot(rho_cond_En, r/1e3, '-', 'color', colors(4, :), 'linewidth', 2);
plot(rho_sol_En, r/1e3, '--', 'color', colors(4, :), 'linewidth', 2);

patch([0 1e10 1e10 0 0], [2402 2402 2440 2440 2402], 'k', 'EdgeAlpha', 0, 'FaceAlpha', .1);
hold off
set(gca, 'fontsize', 20);
ylabel("Radius (km)");
xlabel("Density (kg/m^3)");
ylim([1990 2440]);
xlim([2400 3300]);

saveas(gcf, 'figures/mantle_rho.png');

%%

X_all = [.58 .35 .07];
V_CaS = EOS_params_CaS.molar_mass / (EOS_params_CaS.rho_0/1e3);
V_Na2S = molar_mass(["Na" 2; "S" 1]) / 1.86;
V_MnS = molar_mass(["Mn" 1; "S" 1]) / 3.99;
V_all = [V_CaS V_Na2S V_MnS];
EOS_params_S3.molar_mass = 75.2;
EOS_params_S3.T_0 = 298;
EOS_params_S3.rho_0 = 1000 * EOS_params_S3.molar_mass / sum(V_all .* X_all);
EOS_params_S3.K = 56e9;
EOS_params_S3.K_d = 5.3;
EOS_params_S3.alpha_0 = 4.8e-5;
EOS_S3 = @(T, P) EOS_bm3_rho(EOS_params_S3, T, P);


R_c = 1990e3;
P_c = 5.5e9;
P_of_r = @(r) P_c * (2440e3 - r) / (2440e3 - R_c);

T_hi = @(r) 1800;
T_lo = @(r) 1300;

r = flip(1990e3:1e3:2440e3);
[rho_cond_S3, rho_sol_S3] = fn_get_rho_profiles(...
    r, EOS_S3, P_of_r, T_lo, T_hi);
fprintf("Volume-averaged densities for 1300 K (1800 K) profile:\n");
fprintf("\tS3: %g (%g)\n", -V_mean(r, rho_cond_S3), -V_mean(r, rho_sol_S3));

T_solidus = @(P) 1421 + 177 .* P ./ 1e9 - 12.2 * (P ./ 1e9).^2;
T_solidus_of_r = @(r) T_solidus(P_of_r(r));
[r_conductive, ~, T_conductive] = fn_conductive_heat_loss_spherical(...
    440, 1920, [4 3 .2], [0 6.8 6.8]*1e-8, [1990e3 2402e3 2435e3 2440e3]);
T_conductive_of_r = @(r) interp1(r_conductive, T_conductive, r);

r = flip(1990e3:1e3:2440e3);
[rho_cond_S3, rho_sol_S3] = fn_get_rho_profiles(...
    r, EOS_S3, P_of_r, T_conductive_of_r, T_solidus_of_r);
fprintf("Volume-averaged densities for conductive (solidus) profile:\n");
fprintf("\tS2: %g (%g)\n", -V_mean(r, rho_cond_S3), -V_mean(r, rho_sol_S3));



%% Functions

function [rho_conductive, rho_solidus]=fn_get_rho_profiles(...
    r, EOS, P_of_r, T_conductive_of_r, T_solidus_of_r)
    rho_solidus = nan(size(r));
    rho_conductive = nan(size(r));
    for i=1:length(r)
        rho_solidus(i) = EOS(T_solidus_of_r(r(i)), P_of_r(r(i)));
        rho_conductive(i) = EOS(T_conductive_of_r(r(i)), P_of_r(r(i)));
    end
end

function rho=fn_get_rho_profiles_P_T(EOS, P, T)
    rho = nan(length(P), length(T));
    for i=1:length(T)
        rho(:, i) = EOS(T(i), P);
    end
end

function P=EOS_bm3(EOS_params, T, rho_norm_0)
    K = EOS_params.K;
    K_d = EOS_params.K_d;
    alpha_0 = EOS_params.alpha_0;
    dT = T - EOS_params.T_0;
    
    rho_norm = rho_norm_0 * (1 + alpha_0 * dT);
    P = 3/2*K * (rho_norm.^(7/3) - rho_norm.^(5/3)) .* ...
        (1 + 3/4*(K_d - 4) * (rho_norm.^(2/3) - 1));
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

function X_mean = V_mean(r, X)
    X_total = trapz(r, 4*pi*r.^2.*X);
    V = 4/3*pi*(max(r)^3 - min(r)^3);
    X_mean = X_total/V;
end
