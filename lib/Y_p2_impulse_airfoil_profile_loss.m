function Y_p2 = Y_p2_impulse_airfoil_profile_loss(alpha2, s_c_design)
%Y_p2_impulse_airfoil_profile_loss Calculates impulse airfoil profile loss Y_p2

% Take absolute value of alpha2
alpha2 = abs(alpha2);

% Determine location of s_c_min
    s_c_min = 0.224 + 1.575 * (alpha2 / 90) - (alpha2 / 90)^2; 

% Calculate difference between actual and design s_c_min
X = s_c_design - s_c_min;

% Calculate algebraic coefficients
A = 0.242 - alpha2 / 151 + (alpha2/127)^2;
if alpha2 <= 30
    B = 0.3 + (30 - alpha2) / 50;
else
    B = 0.3 + (30 - alpha2) / 275;
end
C = 0.88 - alpha2 / 42.4 + (alpha2 / 72.8)^2;

% Calculate Y_p2 coefficient
Y_p2 = A + B*X.^2 - C*X.^3;
end

