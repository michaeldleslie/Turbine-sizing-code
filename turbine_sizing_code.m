%% Initialization
clear variables;
% Add lib to path
addpath(genpath(pwd));

%% Set Up Inital Guess Values & Global Variables
eff_tt = 0.9;
flow_coeff = 0.4;
stage_loading = 1.6;
Rc = 0;
t01 = 753;
p01 = 182385;
p3 = 101325;
cp = 1156.9;
gamma = 1.33;
z = 1.00;
mass_flow = 9;
power = 1200*745.7; %Units in Watts
speed = 8000;
r_gas = 286.9;
%% Variables Defined Later in the Notes
nb_stator = 24;
nb_rotor = 40;
%% Power Check
uc0_ratio_guess = sqrt(1/(flow_coeff^2 - (Rc-1)*(4/eff_tt)));
power_check = 0;
while abs((power*1.0095) - power_check)/power > 0.001 %% over shooting power by 0.95%
    t3is = t01 * (p3/p01)^((gamma-1)/gamma);
    c0 = sqrt(2*cp*(t01-t3is));
    U = c0 * uc0_ratio_guess;
    ca = flow_coeff * U;
    delta_his = (c0^2/2) - (ca^2/2);
    work_int = delta_his * eff_tt;
    power_check = work_int * mass_flow;
    Rc = Rc + 0.0001;
    uc0_ratio_guess = sqrt(1/(flow_coeff^2 - (Rc-1)*(4/eff_tt)));
end
%% Alpha 2 Calc
alpha2 = atand((U*ca)/work_int);
c2 = ca / sind(alpha2);
c2u = ca * cotd(alpha2);
w2u = c2u - U;
w2a = ca;
c2a = ca;
w2 = sqrt(w2u^2 + w2a^2);
alpha2_p = atand(w2a/w2u);

w3a = ca;
c3a = ca;
w3u = -U;
w3 = sqrt(w3a^2 + w3u^2);
c3 = c3a;
alpha3_p = atand(w3a/w3u);
%% Determination of Thermodynamic Properties
t2_guess = t01;
t2 = 0;
while abs(t2-t2_guess)/t2_guess > 0.0001
    a2 = sqrt(gamma*r_gas*t2_guess);
    m2 = c2/a2;
    t2_guess = t01/(1+((gamma-1)/2)*m2^2);
    t2 = t2 + 0.01;
end
%% Relative Stator 2 Outlet Calculations
mw2 = w2/a2;
tw2 = t2*(1+((gamma-1)/2)*mw2^2);
%% Rotor Outlet Relative Conditions
t3_guess = tw2;
t3 = 0;
while abs(t3 -t3_guess)/t3_guess > 0.0001
    a3 = sqrt(gamma*r_gas*t3_guess);
    mw3 = w3 / a3;
    t3_guess = tw2 / (1 + ((gamma-1)/2) * mw3^2);
    t3 = t3 + 0.01;
end
%% P3 Outlet Pressure 2% Inital
p02_init = 0.98*p01;
p2_init = p02_init / (1 + ((gamma-1)/2) * m2^2)^(gamma/(gamma-1));
pw2_init = p2_init * (1 + ((gamma-1)/2) * mw2^2)^(gamma/(gamma-1));
pw3_init = 0.98 * pw2_init;
p3_veri_init = pw3_init / (1 + ((gamma-1)/2) * mw3^2)^(gamma/(gamma-1));
%% P3 Outlet Pressure/Loss Verification
loss_percentage = 0.02;
loss = 1 - loss_percentage;
p3_veri_old = 0;
p3_veri = 100000;
count = 0;
while (abs(p3_veri - p3)/p3) > 0.001
    count = count + 1;
    loss = 1 - loss_percentage;
    p02 = loss*p01;
    p2 = p02 / (1 + ((gamma-1)/2) * m2^2)^(gamma/(gamma-1));
    pw2 = p2 * (1 + ((gamma-1)/2) * mw2^2)^(gamma/(gamma-1));
    pw3 = loss * pw2;
    p3_veri = pw3 / (1 + ((gamma-1)/2) * mw3^2)^(gamma/(gamma-1));
    loss_percentage = loss_percentage + 0.0001;
end
%% Density Calc at Rotor Inlet
rho2 = p2/(z*r_gas*t2);
%% Density Calc at Rotor Exit
rho3 = p3/(z*r_gas*t3);
%% Annulus Flow Areas
area2 = mass_flow/(rho2*c2a);
area3 = mass_flow/(rho3*c3a);
%% Sizing Calcs
r_mean = (60 * U)/(2*pi*speed);
BL = 0.9; %Blockage Value
l2 = area2/(2*pi*r_mean*BL);
l2 = area3/(2*pi*r_mean*BL);
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
%% Equation 63 ?
omega = r_mean*c2u;