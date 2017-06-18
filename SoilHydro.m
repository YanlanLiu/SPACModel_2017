% This function calculates soil properties based on soil composition and
% moisture, using the empirical relations provided in the following reference:
% Saxton, K. E., et al. "Estimating generalized soil-water characteristics
% from texture." Soil Science Society of America Journal 50.4 (1986): 1031-1036.

% Input:
% pct: Soil composition, 1 by 3, percentages of sand, silt and clay, unitless
% s: relative soil moisture, unitless
% SID: soil texture ID, unitless

% Output:
% psi: soil water potential, MPa
% K: soil conductivity, m/s
% thetas: soil porosity, unitless,
% sh: hygroscopic point, unitless
% sw: wilting point, unitless
% sfc: field capacity, unitless
% bb: nonlinear factor of the soil retention curve, unitless

function [psi,K,thetas,sh,sw,sfc,bb] = SoilHydro(pct,s,SID)
S = pct(1);
C = pct(3);
a = -4.396;
b = -0.0715;
c = -4.880e-4;
d = -4.285e-5;
e = -3.140;
f = -2.22e-3;
g = -3.484e-5;
h = 0.332;
j = -7.251e-4;
k = 0.1276; 
m = -0.108;
n = 0.341;

thetas = h+j*S+k*log10(C);
theta = s.*thetas;
A = exp(a+b*C+c*S^2+d*S^2*C)*100;
B = e+f*C^2+g*S^2+g*S^2*C+g*S^2*C;
bb = -B;
psi = A.*theta.^B;
psie = 100*(m+n*thetas);
theta10 = (10/A)^(1/B);
sfc = (10/A)^(1/B)/thetas;
sw = (3000/A)^(1/B)/thetas;
sh = (10000/A)^(1/B)/thetas;

id = find(psi<10);
for i = 1:length(id)
	psi(id) = 10-(theta(id)-thetas)./(thetas-theta10).*(10-psie);
end
psi = -psi./1e3; % kPa->MPa

KSAT = 10.^[-0.13,0.82,0.30,-0.32,-0.40,-0.20,-0.46,-0.54,0.01,-0.72,-0.86,-0.42].*0.0254*24; % inch/hr -> m/d
K = KSAT(SID).*s.^(2*bb+3)/24/3600; % m/s
end