function Y_p1 = Y_p1_nozzle_profile_loss(alpha_outlet, s_c_design)
%Y_p1_nozzle_profile_loss Calculates nozzle profile loss Y_p1

% Take absolute value of alpha2
alpha_outlet = abs(alpha_outlet);

% Determine location of s_c_min 
if alpha_outlet <= 30
    s_c_min = 0.46 + alpha_outlet ./ 77;
else
    s_c_min = 0.614 + alpha_outlet ./ 130;
end

% Calculate difference between actual and design s_c_min
X = s_c_design - s_c_min;

% Calculate algebraic coefficients
if alpha_outlet <= 27
    A = 0.025 + (27 - alpha_outlet) ./ 530;
else
    A = 0.025 + (27 - alpha_outlet) ./ 3085;
end
B = 0.1583 - alpha_outlet ./ 1640;
C = 0.08 .* ((alpha_outlet ./ 30).^2 - 1);
n = 1 + alpha_outlet ./ 30;

% Calculate Y_p1 coefficient
if alpha_outlet <= 30
    Y_p1 = A + B .* X.^2 + C .* X.^3;
else
    Y_p1 = A + B .* abs(X).^n;
end
end

