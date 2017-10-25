function K_p = K_p_compressibility_correction(Mach_1, Mach_2)
%K_p_compressibility_correction Calculates compressiblity correction K_p

% Calculate modified inlet Mach
Mach_bar_1 = (Mach_1 + 0.566 - abs(0.566 - Mach_1)) / 2;

% Calculate modified discharge Mach
Mach_bar_2 = (Mach_2 + 1 - abs(Mach_2 - 1)) / 2;

% Calculate correction factor K_p
X = (2 * Mach_bar_1) / (Mach_bar_1 + Mach_bar_2 + abs(Mach_bar_2 - Mach_bar_1));
K_1 = 1 - 0.625 * (Mach_bar_2 - 0.2 + abs(Mach_bar_2 - 0.2));
K_p = 1 - (1 - K_1) * X^2;

end