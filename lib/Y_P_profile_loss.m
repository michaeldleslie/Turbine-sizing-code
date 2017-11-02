function Y_P = Y_P_profile_loss(alpha_inlet, alpha_outlet, s_c_design)
%Y_P_profile_loss Calculates profile loss Y_P by interpolating between
%values for Y_p1 and Y_p2

% Take absolute value of alpha2
alpha_outlet = abs(alpha_outlet);

% Execute subfunctions for Y_p1 and Y_p2
Y_p1 = Y_p1_nozzle_profile_loss(alpha_outlet, s_c_design);
Y_p2 = Y_p2_impulse_airfoil_profile_loss(alpha_outlet, s_c_design);

% Calculate interpolation coefficient
interpolation_coeff = (90 - alpha_inlet) / (90 - alpha_outlet);

% Find final Y_P value
Y_P = Y_p1 + interpolation_coeff^2 * (Y_p2 - Y_p1);

end