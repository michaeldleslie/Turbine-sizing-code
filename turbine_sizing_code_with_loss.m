%% Initialization
clc
clear all; %#ok<CLALL>
% Add lib to path
addpath(genpath(pwd));

%% Set Up Inital Guess Values & Global Variables
% Total-to-total efficiency
eff_tt = 0.85;
flow_coefficient = 0.4;
stage_loading = 1.6;
Rc = 0.4;
alpha1 = 90;
t01 = 500;
p01 = 13842074;
p3 = 11011992;
cp = 15636.11;
gamma = 1.3994;
z = 1.00;
mass_flow = 11.67;
power_required = 6215*745.7; %Units in Watts
speed = 50000;
r_gas = 4123.311;
nb_stator = 90;
nb_rotor = 100;
K_loss_N = 0.98;
K_loss_R = 0.9604;
blockage = 0.9;
eff_tt_old = 0;
count1 = 0;
m1 = 0.05;
r_trailingedge = 0.00025;


while (abs(eff_tt - eff_tt_old) > 0.001)
    %% Check if required power is achieved
    
    % Velocity ratio guess (Equation (9))
    uc0_ratio = sqrt(1/(flow_coefficient^2 - (Rc-1)*(4/eff_tt)));
    % Find isentropic discharge temperature
    t3_isentropic = t01 * (p3/p01)^((gamma-1)/gamma);
    % Find spouting velocity
    c0 = sqrt(2*cp*(t01-t3_isentropic));
    U = c0 * uc0_ratio;
    ca = flow_coefficient * U;
    delta_h_isentropic = (c0^2/2) - (ca^2/2);
    work_intensive = delta_h_isentropic * eff_tt;
    power = work_intensive * mass_flow;
    
    if (power < power_required)
        fprintf('Power requirements not met.\n');
        return;
    end
    
    %% Calculate velocity triangles and alpha2 for station 2 and 3
    % Calculates velocity triangles and alpha2
    
    % Inlet
    alpha2 = atand((U*ca)/work_intensive);
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
    % Equation 42
    p02 = K_loss_N*p01;
    % Convert to static absolute
    p2 = p02 / (1 + ((gamma-1)/2) * m2^2)^(gamma/(gamma-1));
    % Convert to static relative
    pw2 = p2 * (1 + ((gamma-1)/2) * mw2^2)^(gamma/(gamma-1));
    % Relative pressure from station 2 to 3
    % Equation 45
    pw3 = K_loss_R * pw2;
    % Verification pressure
    p3_verification = pw3 / (1 + ((gamma-1)/2) * mw3^2)^(gamma/(gamma-1));
    
    
    
    %% Density calculations
    rho2 = p2/(z*r_gas*t2);
    rho3 = p3/(z*r_gas*t3);
    
    
    %% Annulus flow areas
    annulus_area_stator = mass_flow/(rho2*c2a);
    annulus_area_rotor = mass_flow/(rho3*c3a);
    
    
    %% Sizing Calcs
    r_mean = (60 * U)/(2*pi*speed);
    l_stator = annulus_area_stator/(2*pi*r_mean*blockage);
    l_rotor = annulus_area_rotor/(2*pi*r_mean*blockage);
    
    
    %% Stator geometry calculations
    s_c_opt_stator = 0.427 + alpha2/58 - (alpha2/93)^2;
    zweiffel_coefficient = 0.8;
    s_bz_stator = zweiffel_coefficient/(2 * sind(alpha2)^2 * (cotd(alpha2) - cotd(alpha1)));
    s_stator = (2*pi*r_mean)/(nb_stator);
    chord_stator = s_stator/s_c_opt_stator;
    axial_chord_stator = s_stator/s_bz_stator;
    t_max_stator = 0.1*chord_stator;
    stagger_angle_stator = asind(axial_chord_stator/chord_stator);
    
    
    %% Rotor geometry calculations
    s_rotor = (2*pi*r_mean)/(nb_rotor);
    s_c_0_rotor = 0.427 + abs(alpha3_p)/58 - (abs(alpha3_p)/93)^2;
    s_c_1_rotor = 0.224 + (1.575 - abs(alpha3_p)/90) * (abs(alpha3_p)/90);
    interpolation_parameter = (90-abs(alpha2_p))/(90-abs(alpha3_p));
    s_c_opt_rotor = s_c_0_rotor + (s_c_1_rotor - s_c_0_rotor)*...
        (abs(interpolation_parameter)*interpolation_parameter);%make sure to keep sign
    s_bz_rotor = zweiffel_coefficient/(2 * sind(alpha3_p)^2 * (cotd(alpha2_p) - cotd(alpha3_p)));
    chord_rotor = s_rotor/s_c_opt_rotor;
    axial_chord_rotor = s_rotor/s_bz_rotor;
    t_max_rotor = 0.1*chord_rotor;
    stagger_angle_rotor = asind(axial_chord_rotor/chord_rotor);
    
    %% Store old efficiency
    % store initial efficiency guess into separate variable
    eff_tt_old = eff_tt;
    
    
    %% Perform loss calculations
    
    % Find necessary quantities
    % Calculate M1
    t1 = t01 - (ca^2/(2*cp));
    
    a3 = sqrt(gamma*r_gas*t3);
    m3 = c3 / a3;
    % Absolute Pressure at nozzle inlet
    p1 = p01 / (1 + ((gamma-1)/2) * m1^2)^(gamma/(gamma-1));
    % Density at nozzle inlet
    rho1 = p1/(z*r_gas*t1);
    
    % Stator
    Y_P_stator = Y_P_profile_loss(alpha1, alpha2, s_c_opt_stator, m1, m2, rho2, c2, chord_stator, t2, t_max_stator);
    Y_S_stator = Y_S_secondary_flow_loss(90, alpha2, alpha1, l_stator/chord_stator, s_c_opt_stator);
    Y_cl_stator = 0;
    Y_stator = Y_P_stator + Y_S_stator + Y_cl_stator;
    
    % Rotor Losses
    Y_S_rotor = Y_S_secondary_flow_loss(alpha2_p, alpha3_p, alpha2_p, l_rotor/chord_rotor, s_c_opt_rotor);
    Y_cl_rotor = Y_cl_clearance_loss(alpha2_p, alpha3_p, s_c_opt_rotor, l_rotor/chord_rotor, 0.02, 1);
    Y_P_rotor = Y_P_profile_loss(alpha2_p, alpha3_p, s_c_opt_rotor, mw2, mw3, rho3, w3, chord_rotor, tw3, t_max_rotor);
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
    
    K_loss_N = p02_iteration / p01;
    K_loss_R = pw3_iteration / pw2;
    
    count2 = 0;
    %% Start inner loop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    while abs(p3 - p3_verification) > 100
        c2_iteration = c2 * (p3_iteration/p3)^2;
        % Equation 11
        ca_iteration = c2_iteration * sind(alpha2);
        ca = ca_iteration;
        
        
        %% Calculate velocity triangles and alpha2 for station 2 and 3
        % Calculates velocity triangles and alpha2
        
        % Inlet
        % Do not recalculate alpha2 on each inner loop iteration
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
        % Equation 42
        p02 = K_loss_N*p01;
        % Convert to static absolute
        p2 = p02 / (1 + ((gamma-1)/2) * m2^2)^(gamma/(gamma-1));
        % Convert to static relative
        pw2 = p2 * (1 + ((gamma-1)/2) * mw2^2)^(gamma/(gamma-1));
        % Relative pressure from station 2 to 3
        % Equation 45
        pw3 = K_loss_R * pw2;
        % Verification pressure
        p3_verification = pw3 / (1 + ((gamma-1)/2) * mw3^2)^(gamma/(gamma-1));
        
        
        
        %% Density calculations
        rho2 = p2/(z*r_gas*t2);
        rho3 = p3/(z*r_gas*t3);
        
        
        %% Annulus flow areas
        annulus_area_stator = mass_flow/(rho2*c2a);
        annulus_area_rotor = mass_flow/(rho3*c3a);
        
        
        %% Sizing Calcs
        r_mean = (60 * U)/(2*pi*speed);
        l_stator = annulus_area_stator/(2*pi*r_mean*blockage);
        l_rotor = annulus_area_rotor/(2*pi*r_mean*blockage);
        
        
        %% Stator geometry calculations
        s_c_opt_stator = 0.427 + alpha2/58 - (alpha2/93)^2;
        zweiffel_coefficient = 0.8;
        s_bz_stator = zweiffel_coefficient/(2 * sind(alpha2)^2 * (cotd(alpha2) - cotd(alpha1)));
        s_stator = (2*pi*r_mean)/(nb_stator);
        chord_stator = s_stator/s_c_opt_stator;
        axial_chord_stator = s_stator/s_bz_stator;
        t_max_stator = 0.1*chord_stator;
        stagger_angle_stator = asind(axial_chord_stator/chord_stator);
        
        
        %% Rotor geometry calculations
        s_rotor = (2*pi*r_mean)/(nb_rotor);
        s_c_0_rotor = 0.427 + abs(alpha3_p)/58 - (abs(alpha3_p)/93)^2;
        s_c_1_rotor = 0.224 + (1.575 - abs(alpha3_p)/90) * (abs(alpha3_p)/90);
        interpolation_parameter = (90-abs(alpha2_p))/(90-abs(alpha3_p));
        s_c_opt_rotor = s_c_0_rotor + (s_c_1_rotor - s_c_0_rotor)*...
            (abs(interpolation_parameter)*interpolation_parameter);%make sure to keep sign
        s_bz_rotor = zweiffel_coefficient/(2 * sind(alpha3_p)^2 * (cotd(alpha2_p) - cotd(alpha3_p)));
        chord_rotor = s_rotor/s_c_opt_rotor;
        axial_chord_rotor = s_rotor/s_bz_rotor;
        t_max_rotor = 0.1*chord_rotor;
        stagger_angle_rotor = asind(axial_chord_rotor/chord_rotor);
        
        
        %% Perform loss calculations
        
        
        % Find necessary quantities
        % Calculate M1
        t1 = t01 - (ca^2/(2*cp));
        
        a3 = sqrt(gamma*r_gas*t3);
        m3 = c3 / a3;
        % Absolute Pressure at nozzle inlet
        p1 = p01 / (1 + ((gamma-1)/2) * m1^2)^(gamma/(gamma-1));
        % Density at nozzle inlet
        rho1 = p1/(z*r_gas*t1);
        
        % Stator
        Y_P_stator = Y_P_profile_loss(alpha1, alpha2, s_c_opt_stator, m1, m2, rho2, c2, chord_stator, t2, t_max_stator);
        Y_S_stator = Y_S_secondary_flow_loss(90, alpha2, alpha1, l_stator/chord_stator, s_c_opt_stator);
        Y_cl_stator = 0;
        Y_stator = Y_P_stator + Y_S_stator + Y_cl_stator;
        
        % Rotor Losses
        Y_S_rotor = Y_S_secondary_flow_loss(alpha2_p, alpha3_p, alpha2_p, l_rotor/chord_rotor, s_c_opt_rotor);
        Y_cl_rotor = Y_cl_clearance_loss(alpha2_p, alpha3_p, s_c_opt_rotor, l_rotor/chord_rotor, 0.00025/l_rotor, 1);
        Y_P_rotor = Y_P_profile_loss(alpha2_p, alpha3_p, s_c_opt_rotor, mw2, mw3, rho3, w3, chord_rotor, tw3, t_max_rotor);
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
        
        K_loss_N = p02_iteration / p01;
        K_loss_R = pw3_iteration / pw2;
        
        count2 = count2 + 1;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Calculate stagnation temperatures
    t03 = t3*(1+((gamma-1)/2)*m3^2);
    t02 = t01;
    %% Calculate relative stagnation temp/pressure
    % Find relative stagnation temperatures
    tw03 = tw3*(1+((gamma-1)/2)*mw3^2);
    tw02 = tw2*(1+((gamma-1)/2)*mw2^2);
    % Find relative stagnation pressures
    pw02 = pw2 / ((tw2/tw02)^(gamma/(gamma-1)));
    pw03 = pw3 / ((tw3/tw03)^(gamma/(gamma-1)));
    %% Recalculate efficiency
    % Equation 12
    eff_tt = (1 - (t03/t01)) / (1 - (p03_iteration/p01)^((gamma-1)/gamma));
    count1 = count1 + 1;
    
end


%% Free Vortex

w3u = U;
w3 = sqrt(w3u^2 + w3a^2);
c3u = 0;

% Find tip and hub radii for stator and rotor
r_tip_stator = r_mean + l_stator/2;
r_hub_stator = r_mean - l_stator/2;
r_tip_rotor = r_mean + l_rotor/2;
r_hub_rotor = r_mean - l_rotor/2;

% Scale metal velocity to stator lengths
U_tip = U * r_tip_rotor / r_mean;
U_hub= U * r_hub_rotor / r_mean;

% Scale tangential fluid velocity
c2u_tip = c2u * r_mean / r_tip_stator;
c2u_hub = c2u * r_mean / r_hub_stator;
c3u_tip = c3u * r_mean / r_tip_rotor;
c3u_hub = c3u * r_mean / r_hub_rotor;

% Magnitude of hub/tip velocites
c2_tip = sqrt(c2u_tip^2 + c2a^2);
c2_hub = sqrt(c2u_hub^2 + c2a^2);
c3_tip = sqrt(c3u_tip^2 + c3a^2);
c3_hub = sqrt(c3u_hub^2 + c3a^2);

% Convert to relative frame
w2u_tip = c2u_tip - U_tip;
w2u_hub = c2u_hub - U_hub;
w3u_tip = c3u_tip - U_tip;
w3u_hub = c3u_hub - U_hub;

% Magnitude of relative hub/tip velocities
w2_tip = sqrt(w2u_tip^2 + w2a^2);
w2_hub = sqrt(w2u_hub^2 + w2a^2);
w3_tip = sqrt(w3u_tip^2 + w3a^2);
w3_hub = sqrt(w3u_hub^2 + w3a^2);

% Finalize alpha2 calcs
alpha2_hub = 90*sign(c2u_hub) - atand(c2u_hub/c3a);
alpha2_p_hub = 90*sign(w2u_hub) - atand(w2u_hub/c3a);
alpha2_tip = 90*sign(c2u_tip) - atand(c2u_tip/c3a);
alpha2_p_tip = 90*sign(w2u_tip) - atand(w2u_tip/c3a);

% Finalize alpha3 calcs
alpha3_hub = 90*sign(c3u_hub) - atand(c3u_hub/c3a);
alpha3_p_hub = 90*sign(w3u_hub) - atand(w3u_hub/c3a);
alpha3_tip = 90*sign(c3u_tip) - atand(c3u_tip/c3a);
alpha3_p_tip = 90*sign(w3u_tip) - atand(w3u_tip/c3a);


%% Calculation of the Stator/Rotor Throat Areas
o_stator = (s_stator * blockage * sind(alpha2));

o_rotor = (s_rotor * blockage * sind(abs(alpha3_p)));

throat_area_stator = o_stator * l_stator;
throat_area_rotor = o_rotor * l_rotor;

%% Station 2 thermodynamic conditions for Hub and Tip

% Tip Calculations
t2_tip = t02 - (c2_tip^2/(2*cp));
a2_tip = sqrt(gamma * r_gas * t2_tip);
m2_tip = c2_tip / a2_tip;
p2_tip = p02 / (1 + ((gamma-1)/2)*m2_tip^2)^(gamma/(gamma-1));

% Hub Calculations
t2_hub = t02 - (c2_hub^2/(2*cp));
a2_hub = sqrt(gamma * r_gas * t2_hub);
m2_hub = c2_hub / a2_hub;
p2_hub = p02 / (1 + ((gamma-1)/2)*m2_hub^2)^(gamma/(gamma-1));

% Relative Tip Calculations
mw2_tip = w2_tip / a2_tip;
pw2_tip = p2_tip * (1 + ((gamma-1)/2)*mw2_tip^2)^(gamma/(gamma-1));
tw2_tip = t2_tip * (1 + ((gamma-1)/2)*mw2_tip^2);

% Relative Hub Calculations
mw2_hub = w2_hub / a2_hub;
pw2_hub = p2_hub * (1 + ((gamma-1)/2)*mw2_hub^2)^(gamma/(gamma-1));
tw2_hub = t2_hub * (1 + ((gamma-1)/2)*mw2_hub^2);

% Throat dimension for Hub and Tip - Stator
o2_tip = (2*pi*r_tip_stator/nb_stator) * sind(alpha2_tip) - 2*r_trailingedge;
o2_hub = (2*pi*r_hub_stator/nb_stator) * sind(alpha2_hub) - 2*r_trailingedge;
o2_mean = (o2_hub + o2_tip) / 2;
throat_stator_hub = o2_hub * l_stator;
throat_stator_tip = o2_tip * l_stator;

% Throat dimension for Hub and Tip - Rotor
o3_tip = (2*pi*r_tip_rotor/nb_rotor) * sind(abs(alpha3_p_tip)) - 2*r_trailingedge;
o3_hub = (2*pi*r_hub_rotor/nb_rotor) * sind(abs(alpha3_p_hub)) - 2*r_trailingedge;
o3_mean = (o3_hub + o3_tip) / 2;
throat_rotor_hub = o3_hub * l_rotor;
throat_rotor_tip = o3_tip * l_rotor;

fprintf('Code complete.\n');