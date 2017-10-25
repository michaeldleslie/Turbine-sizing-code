clear variables;

%% Validate Y_p1
clear variables;
sc = 0.1:0.01:1.2;
alpha = [10 15 20 25 30 40 50];
for i = 1:length(alpha)
    alpha2 = alpha(i);
    yp1 = Y_p1_nozzle_profile_loss(alpha2, sc);
    plot(sc, yp1, sc, yp1, 'o');
    hold on;
end
hold off;
xlabel('s/c');
ylabel('Yp1');
grid;

%% Validate Y_p2
clear 
sc = 0.1:0.01:1.1;
alpha = [10 15 20 25 30 35 40 50];
for i = 1:length(alf)
    alfa2 = alf(i);
    yp2 = Y_p2_impulse_airfoil_profile_loss(alfa2, sc);
    plot(sc, yp2, sc, yp2, 'o');
    hold on;
end
hold off;
xlabel('s/c');
ylabel('Yp2');
grid;

%% Find final profile Y_P loss value
interpolation_coeff = (90 - alpha1) / (90 - alpha2);
Y_P = Y_p1 + interpolation_coeff^2 * (Y_p2 - Y_p1);

%% Correction coefficients
% Modernity correction
K_mod = 1;
% Inclination correction
K_inc = 1;
% Mach correction
K_M = 1;
% Compressibility correction
K_p = compressibility_correction(Mach_1, Mach_2);
% Reynolds correction
K_Re = reynolds_correction(rho, C_2, chord, dynamic_viscosity, roughness);
% Trailing edge correction
K_Te = 1;
% Secondary flow loss coefficient

%% Secondary flow loss coefficient
% Validation parameters
alpha1 = 90;
alpha2 = 18.43;
L_c_ratio = 0.484;
s_c_opt = 0.706;

Y_S = secondary_flow_loss(alpha1, alpha2, alpha1, L_c_ratio, s_c_opt);

%% Supersonic flow loss coefficient Y_EX
Y_EX = 0;

