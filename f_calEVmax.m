% This function calculates potential soil evaporation from the top layer,
% using adapted Penman Equation

% Input:
% T: air temperature, K
% RN: net shortwave radiation, W/m2
% VPD: vapor pressure deficit, kPa
% Const: constants
% Veg: vegetation properties
% Soil: soil properties

% Output:
% EVmax: potential soil evaporation rate, m/s

function EVmax = f_calEVmax(T,RN,VPD,Const,Veg,Soil) 
    VPD = VPD.*1e3; % Pa
    ga = Veg.ga;
    RN = RN*exp(-Veg.LAI*Veg.k);
    gamma=Const.cp.*Const.P0./0.622./Const.lambda; % Pa/K
    T1 = T-273.15; % degree C
    es = 0.6108*exp(17.27*T1/(T1+237.3)).*1e3; % Pa
    epsilon_r = 4098*es/(237.3+T1)^2; % Pa/K
	rsoil = exp(8.206-4.255*Soil.sfc1); % from CLM 3.5
    E1 = (epsilon_r*RN)/(Const.rhow*Const.lambda*(gamma*(1+rsoil*ga)+epsilon_r));
    E2 = (Const.rhoa*Const.cp*ga*VPD)/(Const.rhow*Const.lambda*(gamma*(1+rsoil*ga)+epsilon_r));
    EVmax = E1+E2; 
    
end