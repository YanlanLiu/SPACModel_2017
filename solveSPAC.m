% This function solves all the states and fluxes of SPAC.

% Input:
% dir: output directory;
% CLM: climate time series and conditions;
% Soil: soil properties;
% Veg: vegetation properties;
% IC: initial conditions
% Const: constants used in the model
% NumofDay: model time length in days
% dLAI: daily LAI;
% P: precipitation time series, or an indicator to generate stochastic
% precipitation when sum(P)<0.

% Output: .txt files in the output directory
% AN: net carbon assimilation per ground area, umol/m2/day
% CI: inter cellular CO2 concentration, ppm
% ET: transpiration, m/s
% FF1: soil-root flux in the top soil layer, m/s
% FF2: soil-root flux in the bottom soil layer, m/s
% GS: stomatal conductance to water, m/s
% HAM: Hamiltonian in the stomatal optimization model, i.e., net carbon 
%      gain minus water loss, mol/m2/s
% II.txt: interception loss, mm/day
% LAMBDA: marginal water use efficiency, umol/mol/kPa
% LL2: deep percolation from the bottom soil layer, m/s
% PPT: precipitation, mm/day
% PSIL: leaf water potential, MPa
% PSIR: root water potential, MPa
% PSIS: soil water potential, MPa
% PSIX: xylem water potential, MPa
% QQ: surface runoff due to saturation excess, m/s
% RTT: water flux from the top to the bottom soil layer, m/s
% SM1: soil volumetric moisture in the top soil layer, unitless
% SM2: soil volumetric moisture in the bottom soil layer, unitless
% TEMP: air temperature, K
% VPDD: vapor pressure deficit, kPa
% VV: soil evaporation, m/s

function h = solveSPAC(dir,CLM,Soil,Veg,IC,Const,NumofDay,LAI,P)

s10 = IC.s10; % initial soil moisture in the top layer
s20 = IC.s20; % initial soil moisture in the bottom layer

% Intialization
TEMP = [];VPDD = [];RNN = [];
ET = [];VV1 = [];LL1 = [];LL2 = [];FF1 = [];FF2 = [];RTT = [];QQ = [];
SM1 = []; SM2 = [];PSIX = [];PSIL = [];GS = [];ANN = [];PSIR = [];PSIS = [];Lambda = []; HAM = []; CI = [];Flag = [];
WUEwindow = 24;

% Generate stochastic precipitation
if sum(P)<0 
	sos = 91; eos = 304; gslen = eos-sos+1; % define the growing season (DOY)
    alphan = CLM.alpha; % mean precipitation intensity in the non-growing season, m
    alphag = CLM.alpha; % mean precipitation intensity in the growing season, m
    lambdag = CLM.Rp*CLM.P/gslen/alphag; % mean precipitation frequency in the growing season, /day
    lambdan = (1-CLM.Rp)*CLM.P/(365-gslen)/alphag; % mean precipitation frequency in the non-growing, /day season
    if lambdan > 0.95
        lambdan = 0.95;
        alphan = (1-CLM.Rp)*CLM.P/(365-gslen)/lambdan;
    elseif lambdag > 0.95
        lambdag = 0.95;
        alphag = CLM.Rp*CLM.P/gslen/lambdag;
    end
    if (lambdan>1 || lambdag>1)
        h = -1;
        return;
    else
        P = stochasticP(alphan,lambdan,alphag,lambdag,sos,eos,NumofDay);
    end
end

Porg = P;
II = zeros(length(P),1);
dlmwrite([dir,'PPT.txt'],Porg);

% Calculate soil properties based on soil texture and composition
[~,~,Soil.n1,Soil.sh1,Soil.sw1,Soil.sfc1,Soil.b1] = SoilHydro(Soil.pct1,1,Soil.SID1);
[~,~,Soil.n2,Soil.sh2,Soil.sw2,Soil.sfc2,Soil.b2] = SoilHydro(Soil.pct2,1,Soil.SID2);

for j = 1:NumofDay

if j == 1
    tmppsil = zeros(WUEwindow,1);
end

tday = mod(j,365);
if tday == 0
    tday = 365;
end
Veg.LAI = LAI(tday);
Veg.delta = Veg.LAI*1e-4*1e3; % interception storage, mm
P(j) = max(Porg(j)-Veg.delta,0); % throughfall, mm
II(j) = Porg(j)-P(j); % interception loss, mm

PP = zeros(24,1);
PP(unidrnd(24)) = P(j)*1e-3/3600; % mm within the particular hour -> m/s 

tmp = nanmean(tmppsil); % average leaf water potential in the previous day, MPa

T2ES = @(T) 0.6108.*exp(17.27.*(T-273.15)./(T-273.15+237.3)).*1e3; % saturated water pressure, Pa

IC.thetaf0 = CLM.itheta(tday); % initial potential temperature, K

if P(j)>2 % initial conditions for rainy days
    IC.qf0 = CLM.iqr(tday);
    RN = CLM.Rnr((tday-1)*24+1:tday*24);
else % initial conditions for sunny days
    IC.qf0 = CLM.iqs(tday);
    RN = CLM.Rns((tday-1)*24+1:tday*24);
end

IC.qf0 = min(IC.qf0,T2ES(IC.thetaf0)*0.662/Const.P0);
dtheta = CLM.dtheta(tday);
ABL.gamma_theta = CLM.b0_theta+CLM.b1_theta*dtheta; % lapse rate of potential temperature in the free atmosphere, K/m
ABL.gamma_q = CLM.b0_q+CLM.b1_q*IC.qf0; % lapse rate of specific humidity in the free atmosphere, g/kg/m

% Solve ODEs governing the development of the atmospheric boundary layer,
% with feedbacks from vegetation

T_rain = find(PP>0); % raining time
if (~isempty(T_rain) && RN(T_rain)>0)
    Trange = [0,T_rain].*3600; 
    options = odeset('Refine',10);
    [T1,Y] = ode45(@(t,y) odefun(t,y,s10,s20,tmp,IC,RN,II(j)./12./3600./1e3,ABL,Soil,Veg,Const,CLM.ca),Trange,[IC.h0,IC.thetaf0,IC.qf0],options); 
    ta1=Y(:,2);
    q1=Y(:,3);
    
    % ABL collapses when it rains, starting again from the initial condition
    Trange = [T_rain,24].*3600; 
    options = odeset('Refine',10);
    [T2,Y] = ode45(@(t,y) odefun(t,y,s10,s20,tmp,IC,RN,II(j)./12./3600./1e3,ABL,Soil,Veg,Const,CLM.ca),Trange,[IC.h0,IC.thetaf0,IC.qf0],options); 
    ta2=Y(:,2);
    q2=Y(:,3);
    T = [T1;T2(2:end)];
    ta = [ta1;ta2(2:end)];
    q = [q1;q2(2:end)];
    es = T2ES(ta);
    vpd0=max(es-Const.P0/0.622.*q,0);
    
else
    Trange = [0,1].*24.*3600;
    options = odeset('Refine',10);
    [T,Y] = ode45(@(t,y) odefun(t,y,s10,s20,tmp,IC,RN,II(j)./12./3600./1e3,ABL,Soil,Veg,Const,CLM.ca),Trange,[IC.h0,IC.thetaf0,IC.qf0],options); 
    ta=Y(:,2);
    q=Y(:,3);
    es = T2ES(ta);
    vpd0=max(es-Const.P0/0.622.*q,0);
end

TT = interp1(T,ta,[1:24].*3600,'linear','extrap'); % Hourly air temperature, K
VPD = interp1(T,vpd0,[1:24].*3600,'linear','extrap')./1e3; % Hourly vapor pressure deficit, kPa

% Initialize daily variables
s1 = zeros(24,1); s1(1) = s10; 
s2 = zeros(24,1); s2(1) = s20; 
psis = zeros(24,1); 
psil = zeros(24,1);
psix = zeros(24,1);
psir = zeros(24,1);
gs = zeros(24,1);
V1 = zeros(24,1);
E = zeros(24,1); 
F1 = zeros(24,1);
F2 = zeros(24,1);
L2 = zeros(24,1);
Q = zeros(24,1);
RT = zeros(24,1);
An = zeros(24,1);
lambda = zeros(24,1);
Ham = zeros(24,1);
ci = zeros(24,1);

for i = 1:24
    EVmax = f_calEVmax(TT(i),RN(i),VPD(i),Const,Veg,Soil);
    if RN(i)>0
        tmp = nanmean(tmppsil);
        % Solve carbon flux and states
        [gc,An(i),ci(i),lambda(i)] = f_Carbon(TT(i),RN(i),VPD(i),tmp,CLM.ca,Veg);
        An(i) = An(i)*Veg.LAI;
        gs(i) = 1.6.*gc.*18.*1e-6; % m/s;
        
        % Solve water fluxes and states
        [E(i),F1(i),F2(i),psil(i),psix(i),psir(i),psis(i)] = f_Water(gs(i),VPD(i),s1(i),s2(i),Soil,Veg,Const);
        E(i) = max(E(i),0);
        fe = E(i).*1e6/18; % mol/m2/s; LAI embedded
        fc = An(i)*1e-6; % mol/m2/s;
        ll = lambda(i).*1e-6.*101.325; % umol/mol/kPa -> unitless;
        Ham(i) = fc-ll.*fe;
        if Ham(i)<0
            gc = 0;
            An(i) = 0;
            gs(i) = 1.6.*gc.*18.*1e-6; % m/s;
            [E(i),F1(i),F2(i),psil(i),psix(i),psir(i),psis(i)] = f_Water(gs(i),VPD(i),s1(i),s2(i),Soil,Veg,Const);
            E(i) = max(E(i),0);
            
        end
	else
        E(i) = 0;
        F1(i) = 0;
        F2(i) = 0;
        V1(i) = 0;
        EVmax = 0;
    end
    
	tmppsil = [tmppsil(2:end);psil(i)];
    
    % Calculate soil moisture based on soil water balance
    [ds1,ds2,F1(i),F2(i),V1(i),L2(i),RT(i),Q(i)] = f_calds(s1(i),s2(i),PP(i),F1(i),F2(i),EVmax,Soil,Const);
    s1(i+1) = min(s1(i)+ds1,1);
    s2(i+1) = min(s2(i)+ds2,1);

end

s10 = s1(end);
s20 = s2(end);
s1 = s1(1:end-1);
s2 = s2(1:end-1);

% Arrange data and output
SM1 = [SM1;s1'];
SM2 = [SM2;s2'];
PSIX = [PSIX;psix'];
ET = [ET;E'];
VV1 = [VV1;V1'];
FF1 = [FF1;F1'];
FF2 = [FF2;F2'];
RTT = [RTT;RT'];
LL2 = [LL2;L2'];
QQ = [QQ;Q'];
GS = [GS;gs'];
TEMP = [TEMP;TT'];
RNN = [RNN;RN'];
VPDD = [VPDD;VPD'];
ANN = [ANN;An']; 
PSIL = [PSIL;psil'];
PSIR = [PSIR;psir'];
PSIS = [PSIS;psis'];
Lambda = [Lambda;lambda'];
HAM = [HAM;Ham'];
CI = [CI;ci'];

if (mod(j,365)==0)
   dlmwrite([dir,'TEMP.txt'],TEMP,'-append');
   dlmwrite([dir,'VPDD.txt'],VPDD,'-append');
   dlmwrite([dir,'SM1.txt'],SM1,'-append');
   dlmwrite([dir,'SM2.txt'],SM2,'-append');
   dlmwrite([dir,'GS.txt'],GS,'-append');
   dlmwrite([dir,'ET.txt'],ET,'-append');
   dlmwrite([dir,'VV.txt'],VV1,'-append');
   dlmwrite([dir,'FF1.txt'],FF1,'-append');
   dlmwrite([dir,'FF2.txt'],FF2,'-append');
   dlmwrite([dir,'LL2.txt'],LL2,'-append');
   dlmwrite([dir,'QQ.txt'],QQ,'-append');
   dlmwrite([dir,'RTT.txt'],RTT,'-append');
   dlmwrite([dir,'PSIX.txt'],PSIX,'-append');
   dlmwrite([dir,'PSIL.txt'],PSIL,'-append');
   dlmwrite([dir,'PSIR.txt'],PSIR,'-append');
   dlmwrite([dir,'PSIS.txt'],PSIS,'-append');
   dlmwrite([dir,'AN.txt'],ANN,'-append');
   dlmwrite([dir,'LAMBDA.txt'],Lambda,'-append');
   dlmwrite([dir,'HAM.txt'],HAM,'-append');
   dlmwrite([dir,'CI.txt'],CI,'-append');
   TEMP = [];VPDD = [];RNN = [];
   ET = [];VV1 = [];LL2 = [];FF1 = [];FF2 = [];RTT = [];QQ = [];
   SM1 = []; SM2 = [];PSIX = [];PSIL = [];GS = [];ANN = [];PSIR = [];PSIS = [];Lambda = [];HAM = [];CI = [];
   
end
   
end
dlmwrite([dir,'II.txt'],II);
h = 0;

end