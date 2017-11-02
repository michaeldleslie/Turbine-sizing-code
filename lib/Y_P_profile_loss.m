function Y_P = Y_P_profile_loss(alpha_inlet, alpha_outlet, s_c_design, ...
    Mach_inlet, Mach_outlet, rho_outlet, velocity_outlet, chord, ...
    temperature_outlet, t_max)
%Y_P_profile_loss Calculates profile loss Y_P

% Take absolute value of alphas
alpha_inlet = abs(alpha_inlet);
alpha_outlet = abs(alpha_outlet);


%% Adjustment factors
% K_mod is refinement coefficient, K_mod = to one for both stator and rotor
K_mod = 1;

% K_inc is incidence correction, we are assuming operation at design point
% so K_inc = 1 for both stator and rotor
K_inc = 1;

% K_M is correction for Mach effects, we are operating sub-Mach 1 so we can
% set as K_M = 1
K_M = 1;

% K_p is correction for compressibility
K_p = K_p_compressibility_correction(Mach_inlet, Mach_outlet);

% K_Re is correction for Reynolds effects
roughness = 5e-6;
% Sutherland's law to find viscosity
C_1 = 1.458e-6;
C_2 = 110.4;
dynamic_viscosity_outlet = (C_1 * temperature_outlet ^ 1.5)/(temperature_outlet + C_2);
K_Re = K_Re_reynolds_correction(rho_outlet, velocity_outlet, chord, dynamic_viscosity_outlet, roughness);

% K_TE is correction for trailing edge thickness, set to K_TE = 1 for both
% stator and rotor
K_TE = 1;

%% Find Y_p1 and Y_p2
% Execute subfunctions for Y_p1 and Y_p2
Y_p1 = Y_p1_nozzle_profile_loss(alpha_outlet, s_c_design);
Y_p2 = Y_p2_impulse_airfoil_profile_loss(alpha_outlet, s_c_design);

% Calculate interpolation coefficient
interpolation_coefficient = (90 - alpha_inlet) / (90 - alpha_outlet);


%% Find final Y_P value
Y_P = K_mod * K_inc * K_M * K_p * K_Re * K_TE * ...
    (Y_p1 + interpolation_coefficient^2 * (Y_p2 - Y_p1)) * ...
    (t_max/(0.02*chord))^interpolation_coefficient;


end