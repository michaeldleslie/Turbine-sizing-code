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
alpha1 = 90;
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
eff_tt_old = 0;
count1 = 0;
while (abs(eff_tt - eff_tt_old) > 0.001)
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
    
    tw3 = t3*(1+((gamma-1)/2)*mw3^2);
    
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
    p3_verification = pw3 / (1 + ((gamma-1)/2) * mw3^2)^(gamma/(gamma-1));
    
    
    
    %% Density calculations
    rho2 = p2/(z*r_gas*t2);
    rho3 = p3/(z*r_gas*t3);
    
    
    %% Annulus flow areas
    area2 = mass_flow/(rho2*c2a);
    area3 = mass_flow/(rho3*c3a);
    
    
    %% Sizing Calcs
    r_mean = (60 * U)/(2*pi*speed);
    l_Stator = area2/(2*pi*r_mean*blockage);
    l_Rotor = area3/(2*pi*r_mean*blockage);
    
    
    %% S/C ratio calc stator
    s_c_opt = 0.427 + alpha2/58 - (alpha2/93)^2;
    
    
    %% s/bz calculation
    zweif = 0.8;
    s_bz = zweif/(2 * sind(alpha2)^2 * (cotd(alpha2) - cotd(-87.129)));
    
    
    %% Pitch Calculation and Chord/Axial Chord
    s = (2*pi*r_mean)/(24);
    chord = s/s_c_opt;
    axial_chord = s/s_bz;
    % Ask about t_max
    t_max = 0.1*chord;
    
    %% Stagger angle stator
    stagger_angle = asind(axial_chord/chord);
    
    
    %% S/C ratio calc rotor
    s_c_optRot = 0.427 + alpha2/58 - (alpha2/93)^2;
    s_c_impRot = 0.224 + (1.575 - abs(alpha3_p)/90) * (abs(alpha3_p)/90);
    interpol_param = (90-abs(alpha2_p))/(90-abs(alpha3_p));%seemed a bit high
    s_c_optActRot = s_c_optRot + (s_c_impRot - s_c_optRot)*...
        (interpol_param*interpol_param);%make sure to keep sign
    s_bzRot = zweif/(2 * sind(alpha3_p)^2 * (cotd(alpha2_p) - cotd(alpha3_p)));
    sRot = (2*pi*r_mean)/nb_rotor;
    chordRot = sRot/s_c_optActRot;
    axial_chordRot = sRot/s_bzRot;
    
    
    %% Stagger angle rotor
    stagger_angleRot = asind(axial_chordRot/chordRot);
    
    %% Store old efficiency
    % store initial efficiency guess into separate variable
    eff_tt_old = eff_tt;
    
    
    %% Perform loss calculations
    
    
    % Find necessary quantities
    % Calculate M1
    t1 = t01 - (c0^2/(2*cp));
    a1 = sqrt(gamma*r_gas*t1);
    m1 = c0/a1;
    
    a3 = sqrt(gamma*r_gas*t3);
    m3 = c3 / a3;
    % Absolute Pressure at nozzle inlet
    p1 = p01 / (1 + ((gamma-1)/2) * m1^2)^(gamma/(gamma-1));
    % Density at nozzle inlet
    rho1 = p1/(z*r_gas*t1);
    
    % Stator
    Y_P_stator = Y_P_profile_loss(alpha1, alpha2, s_c_opt, m1, m2, rho2, c2, chord, t2, t_max);
    Y_S_stator = Y_S_secondary_flow_loss(90, alpha2, 90, l_Stator/chord, s_c_opt);
    Y_cl_stator = 0;
    Y_stator = Y_P_stator + Y_S_stator + Y_cl_stator;
    
    % Rotor Losses
    Y_S_rotor = Y_S_secondary_flow_loss(alpha2_p, alpha3_p, alpha2_p, l_Rotor/chord, s_c_optRot);
    Y_cl_rotor = Y_cl_clearance_loss(alpha2_p, alpha3_p, s_c_optRot, l_Rotor/chord, 0.00025/l_Rotor, 0);
    Y_P_rotor = Y_P_profile_loss(alpha2_p, alpha3_p, s_c_opt, mw2, mw3, rho3, w3, chord, tw3, t_max);
    Y_rotor = Y_P_rotor + Y_S_rotor + Y_cl_rotor;
    
    
    %% Find K_loss values
    funcM2 = (1 + (gamma - 1) / 2 * m2^2)^(gamma / (gamma - 1));
    p2_iteration = p01 / (funcM2 * (Y_stator + 1) - Y_stator);
    p02_iteration = p2_iteration * funcM2;
    funcMw3 = (1 + (gamma - 1) / 2 * mw3^2)^(gamma / (gamma - 1));
    funcM3 = (1 + (gamma - 1) / 2 * m3^2)^(gamma / (gamma - 1));
    p3_iteration = pw2 / (funcMw3 * (Y_rotor + 1) - Y_rotor);
    pw3_iteration = p3_iteration * funcMw3;
    p03_iteration = p3_iteration * funcM3;
    
    Kloss_N = p02_iteration / p01;
    Kloss_R = pw3_iteration / pw2;
    
    
%     c2_iteration = c2 * (p3_iteration/p3)^2;
%     ca_iteration = c2_iteration * sind(alpha2);
    
    count2 = 0;
    %% Start inner loop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while abs(p3 - p3_verification) > 100
        c2_iteration = c2 * (p3_iteration/p3)^2;
        ca_iteration = c2_iteration * sind(alpha2);
        %% Check if required power is achieved     
        % Velocity ratio guess (Equation (9))
        uc0_ratio = sqrt(1/(flow_coeff^2 - (Rc-1)*(4/eff_tt)));
        % Find isentropic discharge temperature
        t3is = t01 * (p3/p01)^((gamma-1)/gamma);
        % Find spouting velocity
        c0 = sqrt(2*cp*(t01-t3is));
        U = c0 * uc0_ratio;
        ca = ca_iteration;
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
        c2 = c2_iteration;
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
        
        tw3 = t3*(1+((gamma-1)/2)*mw3^2);
        
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
        p3_verification = pw3 / (1 + ((gamma-1)/2) * mw3^2)^(gamma/(gamma-1));
        
        
        
        %% Density calculations
        rho2 = p2/(z*r_gas*t2);
        rho3 = p3/(z*r_gas*t3);
        
        
        %% Annulus flow areas
        area2 = mass_flow/(rho2*c2a);
        area3 = mass_flow/(rho3*c3a);
        
        
        %% Sizing Calcs
        r_mean = (60 * U)/(2*pi*speed);
        l_Stator = area2/(2*pi*r_mean*blockage);
        l_Rotor = area3/(2*pi*r_mean*blockage);
        
        
        %% S/C ratio calc stator
        s_c_opt = 0.427 + alpha2/58 - (alpha2/93)^2;
        
        
        %% s/bz calculation
        zweif = 0.8;
        s_bz = zweif/(2 * sind(alpha2)^2 * (cotd(alpha2) - cotd(-87.129)));
        
        
        %% Pitch Calculation and Chord/Axial Chord
        s = (2*pi*r_mean)/(24);
        chord = s/s_c_opt;
        axial_chord = s/s_bz;
        
        
        %% Stagger angle stator
        stagger_angle = asind(axial_chord/chord);
        
        
        %% S/C ratio calc rotor
        s_c_optRot = 0.427 + alpha2/58 - (alpha2/93)^2;
        s_c_impRot = 0.224 + (1.575 - abs(alpha3_p)/90) * (abs(alpha3_p)/90);
        interpol_param = (90-abs(alpha2_p))/(90-abs(alpha3_p));%seemed a bit high
        s_c_optActRot = s_c_optRot + (s_c_impRot - s_c_optRot)*...
            (interpol_param*interpol_param);%make sure to keep sign
        s_bzRot = zweif/(2 * sind(alpha3_p)^2 * (cotd(alpha2_p) - cotd(alpha3_p)));
        sRot = (2*pi*r_mean)/nb_rotor;
        chordRot = sRot/s_c_optActRot;
        axial_chordRot = sRot/s_bzRot;
        
        
        %% Stagger angle rotor
        stagger_angleRot = asind(axial_chordRot/chordRot);
        
        %% Store old efficiency
        % store initial efficiency guess into separate variable
        eff_tt_old = eff_tt;
        
        
        %% Perform loss calculations
      
        % Find necessary quantities
        % Calculate M1
        t1 = t01 - (c0^2/(2*cp));
        a1 = sqrt(gamma*r_gas*t1);
        m1 = c0/a1;
        
        a3 = sqrt(gamma*r_gas*t3);
        m3 = c3 / a3;
        % Absolute Pressure at nozzle inlet
        p1 = p01 / (1 + ((gamma-1)/2) * m1^2)^(gamma/(gamma-1));
        % Density at nozzle inlet
        rho1 = p1/(z*r_gas*t1);
        
        % Stator
        Y_P_stator = Y_P_profile_loss(alpha1, alpha2, s_c_opt, m1, m2, rho2, c2, chord, t2, t_max);
        Y_S_stator = Y_S_secondary_flow_loss(90, alpha2, 90, l_Stator/chord, s_c_opt);
        Y_cl_stator = 0;
        Y_stator = Y_P_stator + Y_S_stator + Y_cl_stator;
        
        % Rotor Losses
        Y_S_rotor = Y_S_secondary_flow_loss(alpha2_p, alpha3_p, alpha2_p, l_Rotor/chord, s_c_optRot);
        Y_cl_rotor = Y_cl_clearance_loss(alpha2_p, alpha3_p, s_c_optRot, l_Rotor/chord, 0.00025/l_Rotor, 0);
        Y_P_rotor = Y_P_profile_loss(alpha2_p, alpha3_p, s_c_opt, mw2, mw3, rho3, w3, chord, tw3, t_max);
        Y_rotor = Y_P_rotor + Y_S_rotor + Y_cl_rotor;
        
        
        %% Find K_loss values
        funcM2 = (1 + (gamma - 1) / 2 * m2^2)^(gamma / (gamma - 1));
        p2_iteration = p01 / (funcM2 * (Y_stator + 1) - Y_stator);
        p02_iteration = p2_iteration * funcM2;
        funcMw3 = (1 + (gamma - 1) / 2 * mw3^2)^(gamma / (gamma - 1));
        funcM3 = (1 + (gamma - 1) / 2 * m3^2)^(gamma / (gamma - 1));
        p3_iteration = pw2 / (funcMw3 * (Y_rotor + 1) - Y_rotor);
        pw3_iteration = p3_iteration * funcMw3;
        p03_iteration = p3_iteration * funcM3;
        
        Kloss_N = p02_iteration / p01;
        Kloss_R = pw3_iteration / pw2;
                
%         c2_iteration = c2 * (p3_iteration/p3)^2;
%         ca_iteration = c2_iteration * sind(alpha2); 
        count2 = count2 + 1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Recalculate efficiency
    t03 = t3*(1+((gamma-1)/2)*m3^2);
    eff_tt = (1 - t03/t01) / 1 - (p03_iteration/p01)^((gamma-1)/gamma);
    
    
    %% Free Vortex
    % All calculations for rotor
    
    w3u = w3u + U;
    
    r_tip = r_mean - l_Rotor/2;
    r_hub = r_mean + l_Rotor/2;
       
    w2u_tip = w2u * r_mean / r_tip;
    w2u_hub = w2u * r_mean / r_tip;
    w3u_tip = w3u * r_mean / r_tip;
    w3u_hub = w3u * r_mean / r_tip;
    
    alpha2_p_hub = 90 - atand(w2u_hub/c3a);
    alpha2_p_tip = 90 - atand(w2u_tip/c3a);
    alpha3_p_hub = 90 - atand(w3u_hub/c3a);
    alpha3_p_tip = 90 - atand(w3u_tip/c3a);
    count1 = count1 + 1;
end
