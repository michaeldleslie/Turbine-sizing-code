function K_Re = K_Re_reynolds_correction(rho_outlet, velocity_outlet, chord, dynamic_viscosity_outlet, roughness)
%reynolds_correction Calculates Reynolds correction

Re_c = rho_outlet * velocity_outlet * chord / dynamic_viscosity_outlet;

Re_r = 100 * chord / roughness;

if Re_c <= 1E5
    K_Re = sqrt(1E5 / Re_c);
elseif Re_c <= 5E5
    K_Re = 1;
else
    K_Re = 1 + ((log10(5E5)/ log10(Re_r))^2.58 - 1) * (1 - 5E5 / Re_c);
end

end

