% This function calculates soil water balance

% Input
% s1: soil volumetric moisture in the top layer, unitless
% s2: soil volumetric moisture in the bottom layer, unitless
% PP: precipitation, m/s
% F1: soil-root flux in the top layer, m/s
% F2: soil-root flux in the bottom layer, m/s
% EVmax: potential soil evaporation, m/s
% Soil: soil properties
% Const: constants

% Output
% ds1: change in s1, unitless
% ds1: change in s1, unitless
% F1: modified soil-root flux in the top layer considering sh, m/s
% F2: modified soil-root flux in the bottom layer considering sh, m/s
% V1: soil evaporation from the top layer, m/s
% L2: leakage from the bottom layer, m/s
% RT: water flux between the two soil layers, m/s
% Q: surface runoff due to saturation excess, m/s

function [ds1,ds2,F1,F2,V1,L2,RT,Q] = f_calds(s1,s2,PP,F1,F2,EVmax,Soil,Const)

   [psis1,K1] = SoilHydro(Soil.pct1,s1,Soil.SID1);
   [psis2,K2] = SoilHydro(Soil.pct2,s2,Soil.SID2);

    K = (Soil.Zr1+Soil.Zr2)/(Soil.Zr1/K1+Soil.Zr2/K2); % resistance in series
    RT = K*((psis1-psis2).*1e6/Const.rhow/Const.g)/((Soil.Zr1+Soil.Zr2)/2);

    if s1>Soil.sfc1
        V1 = EVmax;
    elseif s1<Soil.sw1
        V1 = 0;
    else
        V1 = (s1-Soil.sw1)/(Soil.sfc1-Soil.sw1)*EVmax;
    end
  
    Loss1 = (F1+V1+RT); % m/s
    dt = 3600;
	ds1 = (PP-Loss1)/(Soil.n1*Soil.Zr1.*1)*dt;    
    

        
    if s2 > Soil.sfc2
        [~,Ksat2] = SoilHydro(Soil.pct2,1,Soil.SID2);
        beta = 2*(Soil.b2)+4;
        L2 = Ksat2/(exp(beta*(1-Soil.sfc2))-1)*(exp(beta*(s2-Soil.sfc2))-1); % deep percolation
    else
        L2 = 0;
    end
    V2 = 0;
    Loss2 = (F2+V2+L2); 
    dt = 3600; % time step in seconds
    ds2 = (RT-Loss2)/(Soil.n2*Soil.Zr2.*1)*dt;
    

    ss1 = s1+ds1;
    ss2 = s2+ds2;
    Q = 0;
    if ss1>Soil.sfc1
        ds1 = Soil.sfc1-s1;
        Q = (ss1-Soil.sfc1)*Soil.n1*Soil.Zr1/dt; % saturated excess
    elseif ss1<Soil.sh1
        ds1 = Soil.sh1-s1;
        F1 = 0;
        RT = 0;
        V1 = 0;
       
    end
    if ss2>Soil.sfc2
        ds2 = Soil.sfc2-s2;
    elseif ss2<Soil.sh2
        ds2 = Soil.sh2-s1;
    end

end