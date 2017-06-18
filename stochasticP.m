% This function generates a time series of stochastic precipitation as a
% marked Poisson process

% Input
% alphan: mean precipitation intensity of the non-growing season, mm
% lambdan: mean frequency of precipitation in the non-growing season, /day
% alphag: mean precipitation intensity of the growing season, mm
% lambdag: mean frequency of precipitation in the growing season, /day 
% sos: start day of the growing season, DOY
% eos: end day of the growing season, DOY
% NumofDay: length of the time series, day

% Output
% P: precipitation time series, mm/day;

function P = stochasticP(alphan,lambdan,alphag,lambdag,sos,eos,NumofDay)

% non-growing season
h = exprnd(alphan,NumofDay,1); % Precipitation intensity follows an exponential distribution
tau = poissrnd(1/lambdan,NumofDay,1); % Interval between two precipitation events follow a Poisson distribution
cumtau = cumsum(tau);
id = find(cumtau>NumofDay,1);
h = h(1:id-1);
cumtau = cumtau(1:id-1);
cumtau(cumtau==0) = 1;
P = zeros(NumofDay,1);
for i = 1:length(h)
    P(cumtau(i)) = h(i)+P(cumtau(i));
end

% growing season
h = exprnd(alphag,NumofDay,1);
tau = poissrnd(1/lambdag,NumofDay,1);
cumtau = cumsum(tau);
id = find(cumtau>NumofDay,1);
h = h(1:id-1);
cumtau = cumtau(1:id-1);
cumtau(cumtau==0) = 1;
P2 = zeros(NumofDay,1);
for i = 1:length(h)
    P2(cumtau(i)) = h(i)+P2(cumtau(i));
end

% combine together
gslen = eos-sos+1;
for i = 1:floor(NumofDay/365)
    P((i-1)*365+sos:(i-1)*365+eos) = P2((i-1)*gslen+1:i*gslen); % mm/d
end

end