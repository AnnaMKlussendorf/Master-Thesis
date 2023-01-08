function [COREstack_second, COREstack_first, COREstack_minor, COREstack_major] = myclasscomparison(age_trans, COREd18Otimed, CORE, lowlim, uplim, sm, sm_max, sm_min)

% set color codes
cfirst = [0.8510 0.7294 0.6039];
csecond = [0.4902 0.2039 0.1412];
cfirst = [ 0.9294    0.6941    0.1255];
csecond = [0.4902 0.2039 0.1412];
cminor = [0.7686 0.5333 0.9882];
cmajor = [0.3412 0.0157 0.3255];

% Different classifications of DO events
allDO = 1:numel(age_trans);
secondhalf = 1:23; % to compare with Buizert et al. - second half of last ice age
firsthalf = 24:numel(age_trans);
HDOs = [2 3 5 10 14 16 22 25 26 27 29 30]; % DOs following a Heinrich Event (H-1 to H-6)
majorDO = HDOs; %  DOs following a Heinrich Event (H-1 to H-6) + 19.2, 20c, 21.1e showing large temperature variations and long interstadials
minorDO = setdiff(allDO, majorDO);


COREstack_all = mean(COREd18Otimed(allDO,:),'omitnan');
COREstack_second = mean(COREd18Otimed(secondhalf,:),'omitnan');
COREstack_first = mean(COREd18Otimed(firsthalf,:),'omitnan');
COREstack_minor = mean(COREd18Otimed(minorDO,:),'omitnan');
COREstack_major = mean(COREd18Otimed(majorDO,:),'omitnan');

t = -1200:1200;

% FIGURE WITH ALL
figure('units','centimeter','position',[0 2 40 15], 'Visible','off')
plot(-t,COREstack_all,'Color',[0 0.4470 0.7410 0.3], 'HandleVisibility','off'); hold on
plot(-t,smoothdata(COREstack_all,'movmean',sm),'Color',[0 0.4470 0.7410],'LineWidth',1.2)
plot(-t,COREstack_second,'Color',[csecond 0.3],'HandleVisibility','off')
plot(-t,smoothdata(COREstack_second,'movmean',sm),'Color',csecond)
plot(-t, COREstack_first,'Color',[cfirst 0.3],'HandleVisibility','off')
plot(-t,smoothdata(COREstack_first,'movmean',sm_min),'Color',cfirst)
plot(-t,COREstack_major,'Color',[cmajor 0.3],'HandleVisibility','off')
plot(-t,smoothdata(COREstack_major,'movmean', sm_max),'Color',cmajor)
plot(-t, COREstack_minor,'Color',[cminor 0.3],'HandleVisibility','off')
plot(-t, smoothdata(COREstack_minor,'movmean',sm),'Color',cminor)
xline(0,'Color',[0.5 0.5 0.5],'LineWidth',2,'Alpha',.4, 'HandleVisibility','off'); hold off %[0.9290 0.6940 0.1250]
legend('All DOs','10-65 ka b2k','65-115 ka b2k','Major','Minor','Location','northeastoutside')
ylabel(['\delta^{18}O [',char(8240),']']) 
xlabel('Time [Years]')
xlim([t(1) t(end)]) 
ylim([lowlim uplim])
ax = gca;
ax.XTick = t(1):200:t(end);
ax.XMinorTick='on';
ax.XAxis.MinorTickValues = t(1):50:t(end);
grid minor
ax.MinorGridLineStyle = '-';
ax.MinorGridAlpha = 0.1;
title(CORE)
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
text(NW(1), NW(2), CORE, 'VerticalAlignment','top', 'HorizontalAlignment',...
    'left', 'FontSize',10, 'color', 'k')



% 2 PLOTS FOR DIRECT COMPARISON OF TWO GROUPS
figure('units','centimeter','position',[0 2 40 15]) 
subplot(1,2,1) % first vs second half
plot(-t,COREstack_second,'Color',[csecond 0.3],'HandleVisibility','off'); hold on
plot(-t,smoothdata(COREstack_second,'movmean',sm),'Color',csecond)
plot(-t, COREstack_first,'Color',[cfirst 0.3],'HandleVisibility','off')
plot(-t,smoothdata(COREstack_first,'movmean',sm_min),'Color',cfirst)
xline(0,'Color',[0.5 0.5 0.5],'LineWidth',2,'Alpha',.4, 'HandleVisibility','off'); hold off %[0.9290 0.6940 0.1250]
legend('10-65 ka b2k','65-115 ka b2k')
ax = gca;
ax.XTick = t(1):200:t(end);
ax.XAxis.MinorTickValues = t(1):50:t(end);
ax.FontSize = 8;
ylabel(['\delta^{18}O [',char(8240),']'],'FontSize',10) 
xlabel('Time [Years]','FontSize',10)
xlim([-800 800]) 
ylim([lowlim uplim])
set(gca, 'YMinorTick','on','XMinorTick','on')
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
text(NW(1), NW(2), CORE, 'VerticalAlignment','top', 'HorizontalAlignment',...
    'left', 'FontSize',8, 'color', 'k')

subplot(1,2,2) % minor vs major
plot(-t,COREstack_major,'Color',[cmajor 0.3],'HandleVisibility','off'); hold on
plot(-t,smoothdata(COREstack_major,'movmean',sm_max),'Color',cmajor)
plot(-t, COREstack_minor,'Color',[cminor 0.3],'HandleVisibility','off')
plot(-t, smoothdata(COREstack_minor,'movmean',sm),'Color',cminor)
xline(0,'Color',[0.5 0.5 0.5],'LineWidth',2,'Alpha',.4, 'HandleVisibility','off'); hold off %[0.9290 0.6940 0.1250]
legend('Major','Minor')
ax = gca;
ax.XTick = t(1):200:t(end);
ax.XAxis.MinorTickValues = t(1):50:t(end);
ax.FontSize = 8;
ylabel(['\delta^{18}O [',char(8240),']'], 'FontSize',10) 
xlabel('Time [Years]','FontSize',10)
xlim([-800 800])
ylim([lowlim uplim])
set(gca, 'YMinorTick','on','XMinorTick','on')
fig_name = strcat('C:\Users\annak\OneDrive - University of Copenhagen\Desktop\Speciale\Figures\classcomp_',CORE);

print(gcf,fig_name,'-dpng','-r600');

