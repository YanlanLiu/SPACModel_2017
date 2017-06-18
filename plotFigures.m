function plotFigures(dir,Veg)

GS = importdata([dir,'GS.txt']); 
figure(1);clf; % plot hourly stomatal conductance
plot(reshape(GS',[],1),'-k');
set(gca,'xlim',[180,190].*24,'FontSize',16,'FontName','Times New Roman')
xlabel('Hour');
ylabel('g_{s,H_2O} (m/s)')

PSIS = importdata([dir,'PSIS.txt']);
PSIR = importdata([dir,'PSIR.txt']);
PSIX = importdata([dir,'PSIX.txt']);
PSIL = importdata([dir,'PSIL.txt']);
figure(2);clf; % plot water potential from soil, root, stem to leaf
plot(min(PSIS,[],2));
hold on;
plot(min(PSIR,[],2));
plot(min(PSIX,[],2));
plot(min(PSIL,[],2));
xlim = get(gca,'xlim');
plot(xlim,Veg.psix50.*[1,1],'--k');
set(gca,'xlim',xlim,'FontSize',16,'FontName','Times New Roman')
legend('\psi_{soil}','\psi_{root}','\psi_{xylem}','\psi_{leaf}','\psi_{50}',...
    'Location','SouthWest')
xlabel('Day');
ylabel('\psi_{min} (MPa)')



end