function Y_P = Y_P_profile_loss(alpha1, alpha2, s_c_design)
%Y_P_profile_loss Calculates profile loss Y_P by interpolating between
%values for Y_p1 and Y_p2

% Take absolute value of alpha2
alpha2 = abs(alpha2);

% Execute subfunctions for Y_p1 and Y_p2
Y_p1 = Y_p1_nozzle_profile_loss(alpha2, s_c_design);
Y_p2 = Y_p2_impulse_airfoil_profile_loss(alpha2, s_c_design);

% Calculate interpolation coefficient
interpolation_coeff = (90 - alpha1) / (90 - alpha2);

% Find final Y_P value
Y_P = Y_p1 + interpolation_coeff^2 * (Y_p2 - Y_p1);

end