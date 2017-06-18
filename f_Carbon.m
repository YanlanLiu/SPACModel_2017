% This function solves CO2 flux through stomata
% Input:
% T: temperature, K
% Rn: net shortwave radiation, W/m2
% VPD: vapor pressure deficit, kPa
% psil: average leaf water potential in the previous day, MPa
% ca: atmospheric CO2 concentration, ppm
% Veg: vegetation properties

% Output:
% gc: stomatal conductance to CO2, mol/m2/s
% An: net carbon assimilation per leaf area, umol/m2/s
% ci: intercellular CO2 concentration, ppm
% lambda:  marginal water use efficiency, umol/mol/kPa

function [gc,An,ci,lambda] = f_Carbon(T,Rn,VPD,psil,ca,Veg) 

if Veg.PFT==1
    koptj = 155.76; % umol/m2/s
    Haj = 43.79; % kJ/mol
    Hdj = 200; % kJ/mol
    Toptj = 32.19+273.15; % K
    
    koptv = 174.33; % umol/m2/s
    Hav = 61.21; % kJ/mol
    Hdv = 200; % kJ/mol
    Toptv = 37.74+273.15; % K
    
else
    koptj = 169.66; 
    Haj = 43.36;
    Hdj = 200;
    Toptj = 35.24+273.15;
    
    koptv = 146.39; 
    Hav = 67.72;
    Hdv = 200;
    Toptv = 38.76+273.15;
end
R = 8.31*1e-3; % Gas constant, kJ/mol/K
Vcmax = koptv*Hdv*exp(Hav*(T-Toptv)/T/R/Toptv)/(Hdv-Hav*(1-exp(Hav*(T-Toptv)/T/R/Toptv))); % umol/m2/s
Jmax = koptj*Hdj*exp(Haj*(T-Toptj)/T/R/Toptj)/(Hdj-Haj*(1-exp(Haj*(T-Toptj)/T/R/Toptj)));

TC = T-273.15; % C
Kc = 300*exp(0.074*(TC-25)); %umol/mol
Ko = 300*exp(0.015*(TC-25)); %mmol/mol
Coa = 210; % mmol/mol
a = 1.6; % unitless

hc = 2e-25; % Planck constant times light speed, J*s times m/s
wavelen = 500e-9; % wavelength of light, m
EE = hc/wavelen; % energy of photon, J
NA = 6.02e23; % Avogadro's constant, /mol
Q = Rn/(EE*NA)*1e6; % absorbed photon irradiance, umol photons /m2/s, PAR
cp = 36.9+1.18*(TC-25)+0.036*(TC-25)^2;
kai1 = 0.9; 
kai2 = 0.3;
J = (kai2*Q+Jmax-sqrt((kai2*Q+Jmax)^2-4*kai1*kai2*Q*Jmax))/2/kai1; % umol electrons /m2/s
lambda = ca/400*Veg.lambdaww/101.325*exp(Veg.beta0*psil);

% Rubisco limited photosynthesis
a1 = J/4;
a2 = 2*cp;
Rd = 0.015*Vcmax;
gc1 = -a1*(a2-ca+2*cp)/(a2+ca)^2+sqrt(a*VPD*lambda*a1^2*(ca-cp)*(a2+cp)*(a2+ca-2*a*VPD*lambda)^2*(a2+ca-a*VPD*lambda))...
    /(a*VPD*lambda*(a2+ca)^2*(a2+ca-a*VPD*lambda)); % nonlinear form (Katul et al., AOB, 2010)
gc1 = max(gc1,0);
A = -gc1;
B = gc1*ca-a2*gc1-a1+Rd;
C = ca*a2*gc1+a1*cp+a2*Rd;
ci = (-B-sqrt(B^2-4*A*C))/(2*A);
ci1 = max(real(ci),0);
An1 = real(gc1*(ca-ci1)); 

% RuBP limited photosynthesis
a1 = Vcmax;
a2 = Kc*(1+Coa/Ko);
Rd = 0.015*Vcmax;
gc2 = -a1*(a2-ca+2*cp)/(a2+ca)^2+sqrt(a*VPD*lambda*a1^2*(ca-cp)*(a2+cp)*(a2+ca-2*a*VPD*lambda)^2*(a2+ca-a*VPD*lambda))...
    /(a*VPD*lambda*(a2+ca)^2*(a2+ca-a*VPD*lambda)); % nonlinear form (Katul et al., AOB, 2010)
gc2 = max(gc2,0);
A = -gc2;
B = gc2*ca-a2*gc2-a1+Rd;
C = ca*a2*gc2+a1*cp+a2*Rd;
ci = (-B-sqrt(B^2-4*A*C))/(2*A);
ci2 = max(real(ci),0);
An2 = real(gc2*(ca-ci2)); 

if min(An1,An2)>0
if An1<An2
    An = An1;
    gc = gc1;
    ci = ci1;
else
    An = An2;
    gc = gc2;
    ci = ci2;
end
else
    An = 0;
    gc = 0;
    ci = ca;
end

gc = real(gc);

end