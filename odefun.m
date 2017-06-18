% This funciton describes the three ODEs governing the the development of
% atmospheric boundary layer (ABL) within a day, influenced by the 
% feedbacks from vegetation.

% Input:
% t: time, second
% y: variables to solve, containing the height (m), potential temperature
%    (K) and specific humidity (g/kg)
% s1: soil moisture in the top layer, unitless
% s2: soil moisture in the bottom layer, unitless
% tmp: average leaf water potential in the previous day, MPa
% IC: initial conditions of potential temperature (K) and specific humidity
%     (g/kg) of the ABL
% Rad: hourly net shortwave radiaiton, W/m2
% II: interception loss, mm/day
% ABL: lapse rates of potential temperature (K) and specific humidity
%     (g/kg) in the free atmosphere
% Soil: soil properties
% Veg: vegetation properties
% Const: constants
% ca: atmospheric CO2 concentration, ppm

% output:
% dy: change in y

function dy = odefun(t,y,s1,s2,tmp,IC,Rad,II,ABL,Soil,Veg,Const,ca)
dy = zeros(3,1);
h=y(1); % heigh of ABL, m
theta=y(2); % potential temperature, K
q=y(3); % specific humidiy, g/kg

es=611.71*exp(2.501/0.461*10^3*(1/273.15-1/theta)); % saturated water pressure, Pa
VPD = max(min(es-q*Const.P0/0.622,es),1)/1e3; % vapor pressure deficit, kPa
RN = interp1([1:24].*3600,Rad,t,'Linear','extrap'); % net shorwave radiation at time t, W/m2

if RN > 0
    gc =  f_Carbon(theta,RN,VPD,tmp,ca,Veg); % stomatal conductance to CO2, mol/m2/s
    gs = 1.6.*gc.*18.*1e-6; % stomatal conductance to H2O, m/s;
    E = f_Water(gs,VPD,s1,s2,Soil,Veg,Const); % transpiration, m/s
    EVmax = f_calEVmax(theta,RN,VPD,Const,Veg,Soil); % potential evaporation from soil, m/s
    
    % actual soil evaporation as a function of EVmax and soil moisture in
    % the top layer, m/s
    if s1>Soil.sfc1
        V = EVmax;
    elseif s1<Soil.sw1
        V = 0;
    else
        V = (s1-Soil.sw1)/(Soil.sfc1-Soil.sw1)*EVmax;
    end
    E = E+V+II;    
else
    E = 0;
end

H = (RN-Const.lambda*Const.rhow*E); % sensible heat, W/m2
if H<10
    H=10;
    E=(RN-H)./(Const.lambda*Const.rhow);
end

%-------------------------------- ODEs ------------------------------------
% dh/dt, development of ABL height
dy(1) = H/(Const.rhoa*Const.cp*ABL.gamma_theta*h);
% dtheta/dt, energy balance
dy(2)=(H+Const.rhoa*Const.cp*(IC.thetaf0+ABL.gamma_theta*h-theta)*dy(1))/(Const.rhoa*Const.cp*h);
% dq/dt, mass balance
dy(3)=(E*Const.rhow+Const.rhoa*(IC.qf0+ABL.gamma_q*h-q)*dy(1))/(Const.rhoa*h);

end



