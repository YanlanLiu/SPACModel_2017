% This function calculates water fluxes and states through the SPAC

% Input:
% gs: stomatal conductance to H2O, m/s
% VPD: vapor pressure deficit, kPa
% s1: soil moisture in the top layer, unitless
% s2: soil moisture in the bottom layer, unitless
% Soil: soil properties
% Veg: vegetation properties
% Const: constants

% Output
% E: transpiration, m/s (per ground area)
% F1: soil-root flux in the top layer, m/s
% F2: soil-root flux in the bottom layer, m/s
% psil: leaf water potential, MPa
% psix: xylem water potential, MPa
% psir: root water potential, MPa 
% psis: soil water potential, MPa
% flag: convergence check, 0-> converge, -1-> not converge, i.e., even the
%       maximum supply cannot meet demand

function [E,F1,F2,psil,psix,psir,psis,flag] = f_Water(gs,VPD,s1,s2,Soil,Veg,Const)
E = gs*VPD/Const.P0*1e3.*Veg.LAI;
E = real(E);

x = 0:-0.1:-8;
gp = Veg.gpmax.*(1./(1+(x./Veg.psix50).^Veg.aa)); % vulnerability curve, plant conductance, m/s
RAI1=Veg.RAIW1; % root area index, unitless
RAI2=Veg.RAIW2; 
[psis1,Kroot1] = SoilHydro(Soil.pct1,s1,Soil.SID1); % unsaturated soil water potential (MPa) and conductivity (m/s)
[psis2,Kroot2] = SoilHydro(Soil.pct2,s2,Soil.SID2);
gsr1 = Kroot1.*RAI1^0.5./Const.g./pi./Const.rhow./Soil.Zr1.*10^6; % soil-root conductance, m/s
gsr2 = Kroot2.*RAI2^0.5./Const.g./pi./Const.rhow./Soil.Zr1.*10^6;
gsr = gsr1+gsr2; % resistance in parallel 
psis = (gsr1.*psis1+gsr2.*psis2)./(gsr1+gsr2); 
gsrp = gsr.*gp./(gsr+gp);
et = gsrp.*(psis-(x+Const.rhow.*Const.g.*Veg.H/1e6)); % water flux within plant, m/s

% check convergence
if max(et-E)>0 % converge
    if et(1)> E
        idl = find(et<E,1,'first');
        psil = interp1(real(et(idl:-1:idl-1)),x(idl:-1:idl-1),real(E));
    else
        idl = find(et>E,1,'first');
        psil = interp1(real(et(idl:-1:idl-1)),x(idl:-1:idl-1),real(E));
    end
    flag = 0;
else % not converge, even the maximum supply cannot meet demand
    psil = x(et==max(et));
    flag = -1; 
end

if flag == 0
    psir = psis-E./gsr;
else
    psir = psis-max(et)./gsr;
end

F1 = gsr1*(psis1-psir);
F2 = gsr2*(psis2-psir);

psix = (psil + psir)/2; % water potential at the middle, MPa

end
