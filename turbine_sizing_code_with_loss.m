%% Initialization
clc
clear variables;
% Add lib to path
addpath(genpath(pwd));


%% Set Up Inital Guess Values & Global Variables
% Total-to-total efficiency
eff_tt = 0.85;
flow_coeff = 0.2;
stage_loading = 1.6;
Rc = 0.6;
t01 = 500;
p01 = 13842074;
p3 = 11011992;
cp = 15636.11;
gamma = 1.3994;
z = 1.00;
mass_flow = 13.61;
power_required = 6215*745.7; %Units in Watts
speed = 50000;
r_gas = 4123.311;
nb_stator = 90;
nb_rotor = 100;
Kloss_N = 0.98;
Kloss_R = 0.9604;
blockage = 0.9;
QQ = true;
while (QQ)
    %% Check if required power is achieved
    
    
    % Velocity ratio guess (Equation (9))
    uc0_ratio = sqrt(1/(flow_coeff^2 - (Rc-1)*(4/eff_tt)));
    % Find isentropic discharge temperature
    t3is = t01 * (p3/p01)^((gamma-1)/gamma);
    % Find spouting velocity
    c0 = sqrt(2*cp*(t01-t3is));
    U = c0 * uc0_ratio;
    ca = flow_coeff * U;
    delta_his = (c0^2/2) - (ca^2/2);
    work_int = delta_his * eff_tt;
    power = work_int * mass_flow;
    
    if (power < power_required)
        fprintf('Power requirements not met.\n');
        return;
    end
    
    %% Calculate velocity triangles and alpha2 for station 2 and 3
    % Calculates velocity triangles and alpha2
    
    % Inlet
    alpha2 = atand((U*ca)/work_int);
    c2 = ca / sind(alpha2);
    c2u = ca * cotd(alpha2);
    w2u = c2u - U;
    w2a = ca;
    c2a = ca;
    w2 = sqrt(w2u^2 + w2a^2);
    alpha2_p = atand(w2a/w2u);
    
    % Discharge
    w3a = ca;
    c3a = ca;
    w3u = -U;
    w3 = sqrt(w3a^2 + w3u^2);
    c3 = c3a;
    alpha3_p = atand(w3a/w3u);
    
    
    %% Determination of thermodynamic properties at inlet
    % Guess is that t2 is equal to static temp t01
    t2 = t01;
    % Initial value to check convergence against
    t2_old = 0;
    % Loop until convergence is achieved
    while abs(t2-t2_old) > 0.0001
        t2_old = t2;
        a2 = sqrt(gamma*r_gas*t2);
        m2 = c2/a2;
        t2 = t01/(1+((gamma-1)/2)*m2^2);
    end
    
    %% Find relative conditions
    % Inlet
    mw2 = w2/a2;
    tw2 = t2*(1+((gamma-1)/2)*mw2^2);
    
    % Outlet
    t3 = tw2;
    t3_old = 0;
    
    % Iterate until t3 is converged
    while abs(t3 -t3_old) > 0.0001
        t3_old = t3;
        a3 = sqrt(gamma*r_gas*t3);
        mw3 = w3 / a3;
        t3 = tw2 / (1 + ((gamma-1)/2) * mw3^2);
    end
    
    
    %% Initial guess of inlet and outlet pressures
    
    % Total pressure guess from station 1 to 2
    p02 = Kloss_N*p01;
    % Convert to static absolute
    p2 = p02 / (1 + ((gamma-1)/2) * m2^2)^(gamma/(gamma-1));
    % Convert to static relative
    pw2 = p2 * (1 + ((gamma-1)/2) * mw2^2)^(gamma/(gamma-1));
    % Relative pressure from station 2 to 3
    pw3 = Kloss_R * pw2;
    % Verification pressure
    p3_veri = pw3 / (1 + ((gamma-1)/2) * mw3^2)^(gamma/(gamma-1));
    
    
    
    %% Density Calc at Rotor Inlet
    rho2 = p2/(z*r_gas*t2);
    
    
    %% Density Calc at Rotor Exit
    rho3 = p3/(z*r_gas*t3);
    
    
    %% Annulus Flow Areas
    area2 = mass_flow/(rho2*c2a);
    area3 = mass_flow/(rho3*c3a);
    
    
    %% Sizing Calcs
    r_mean = (60 * U)/(2*pi*speed);
    %BL = 0.9; %Blockage Value
    l_Stator = area2/(2*pi*r_mean*blockage);
    l_Rotor = area3/(2*pi*r_mean*blockage);
    %% S/C Ratio Calc Stator
    s_c_opt = 0.427 + alpha2/58 - (alpha2/93)^2;
    
    
    %% s/bz Calc NOTE: ASK WHERE THE VALUE FOR ALPHA1 COMES FROM
    zweif = 0.8;
    s_bz = zweif/(2 * sind(alpha2)^2 * (cotd(alpha2) - cotd(-87.129)));
    
    
    %% Pitch Calculation and Chord/Axial Chord
    s = (2*pi*r_mean)/(24);
    chord = s/s_c_opt;
    axial_chord = s/s_bz;
    
    
    %% Stagger Angle Stator
    stagger_angle = asind(axial_chord/chord);
    
    
    %% S/C Ratio Calc Rotor %%Reminder: alpha2 =~ alpha3p
    s_c_optRot = 0.427 + alpha2/58 - (alpha2/93)^2;
    s_c_impRot = 0.224 + (1.575 - abs(alpha3_p)/90) * (abs(alpha3_p)/90);
    interpol_param = (90-abs(alpha2_p))/(90-abs(alpha3_p));%seemed a bit high
    s_c_optActRot = s_c_optRot + (s_c_impRot - s_c_optRot)*...
        (interpol_param*interpol_param);%make sure to keep sign
    s_bzRot = zweif/(2 * sind(alpha3_p)^2 * (cotd(alpha2_p) - cotd(alpha3_p)));
    sRot = (2*pi*r_mean)/nb_rotor;
    chordRot = sRot/s_c_optActRot;
    axial_chordRot = sRot/s_bzRot;
    
    
    %% Stagger Angle Rotor
    stagger_angleRot = asind(axial_chordRot/chordRot);
    
    %% Store old efficiency
    % store initial efficiency guess into separate variable
    eff_tt_old = eff_tt;
    
    
    
    %% Perform Loss Calculations
    % roughness = (100 * axial_chordRot) / Re_rot;
    
    % Calculate M1
    t1 = t01 - (c0^2/(2*cp));
    a1 = sqrt(gamma*r_gas*t1);
    m1 = c0/a1;
    % Typical values for roughness 5E-6, c/e 10,000-20,000
    roughness = 5E-6;
    % Absolute Pressure at nozzle inlet
    p1 = p01 / (1 + ((gamma-1)/2) * m1^2)^(gamma/(gamma-1));
    % Density at nozzle inlet
    rho1 = p1/(z*r_gas*t1);
    % Sutherland?s Law
    dynamic_visc = (5.16E-5)/rho2;
    % % Rotor Losses
    % K_p_R = K_p_compressibility_correction(mw2,mw3)
    % K_Re_R = K_Re_reynolds_correction(rho2,w2,axial_chordRot,dynamic_vis, roughness);
    % %Y_cl_R = Y_cl_clearance_loss(alpha2_p,alpha3_p,s_c_optRot,l3/axial
    % %Y_p1 = Y_p1_nozzle_profile_loss(alpha2,s_c_opt); % 0 for rotors?
    %
    % Y_p2_R = Y_p2_impulse_airfoil_profile_loss
    % Y_P_R = Y_P_profile_loss
    % Y_S_R = Y_S_secondary_flow_loss
    
    % Stator
    Y_p1 = Y_p1_nozzle_profile_loss(alpha2,s_c_opt); % 0 for rotors?
    Y_p2 = Y_p2_impulse_airfoil_profile_loss(alpha2,s_c_opt);
    K_mod_Stator = 1;
    K_inc_Stator = 1;
    K_M_Stator = 1;
    K_P_Stator = K_p_compressibility_correction(m1,m2);
    K_Re_Stator = K_Re_reynolds_correction(rho2,c2,axial_chord,dynamic_visc, roughness);
    K_TE_Stator = 1;
    Y_S_Stator = Y_S_secondary_flow_loss(90,alpha2,90,l_Stator/chord, s_c_opt);
    Y_cl_Stator = 0;
QQ = false;    
end