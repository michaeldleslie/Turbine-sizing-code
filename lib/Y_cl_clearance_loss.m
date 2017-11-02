function Y_cl = Y_cl_clearance_loss( alpha_inlet, alpha_outlet, s_c_ratio, L_c_ratio, tip_clrnc_L_ratio, ishr )
%Y_cl_clearance_loss Calculates the tip clearance loss coefficient
%   F_B is 0.47 for unshrouded airfoils and 0.36 for shrouded
if ishr == 0
    F_B = 0.47;
elseif ishr == 1
    F_B = 0.36;
end

%   Calculate lift coefficient C_L
alpha_m = abs(atand(2 / (cotd(alpha_outlet) + cotd(alpha_inlet))));
C_L = abs(2*(s_c_ratio)*(cotd(alpha_inlet) - cotd(alpha_outlet))*sind(alpha_m));

%   Calculate Ainley loading parameter Z
Z = (C_L/s_c_ratio)^2 * (sind(alpha_outlet)^2/sind(alpha_m)^3);

% Use Z in absolute value
Z = abs(Z);

%   Calculate Y_cl
Y_cl = F_B * Z * L_c_ratio^-1 * tip_clrnc_L_ratio^0.78;

end

