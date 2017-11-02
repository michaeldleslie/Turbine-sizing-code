function Y_S = Y_S_secondary_flow_loss(alpha_inlet, alpha_outlet, beta_inlet, L_c_ratio, s_c_opt)
%Y_S_secondary_flow_loss Calculates secondary flow loss coefficient

% Calculate lift coefficient C_L
alpha_m = abs(atand(2 / (cotd(alpha_outlet) + cotd(alpha_inlet))));
C_L = abs(2 * s_c_opt * (cotd(alpha_inlet) - cotd(alpha_outlet)) * sind(alpha_m));

% Calculate Ainley loading parameter Z
Z = (C_L / s_c_opt)^2 * (sind(alpha_outlet)^2 / sind(alpha_m)^3);

% Use Z in absolute value
Z = abs(Z);

% Calculate cascade aspect correction F_AR
if L_c_ratio > 2
    F_AR = L_c_ratio.^-1;
else
    F_AR = 0.5 * (2 * L_c_ratio.^-1)^0.7;
end

% Calculate Y_S (note, not using Eqn. (35))
Y_S = 0.0334 * F_AR * Z * abs(sind(alpha_outlet)) / abs(sind(beta_inlet));
end