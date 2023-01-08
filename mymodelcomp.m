function mymodelcomp(coord, T, lon, lat, time, CO2, xlims, ylims) %, stacks_comp) %, d18Otimed)

%CORES = {'NGRIP', 'GRIP', 'NEEM', 'GISP2', 'WDC', 'TAL', 'DF', 'EDC', 'EDML'};

sm = 20;

% COLOR CODES 
cNGRIP = [0 0.4470 0.7410];
cGRIP = [0.7350 0.0780 0.1840];
cNEEM = [0.9290 0.6940 0.1250];
cGISP = [0.4660 0.6740 0.1880];
cWDC = [0.4940 0.1840 0.5560];
cTAL = [.2 .6 .4];
cDF = [0.6235 0.8235 0.9882];
cEDC = [0.8902 0.4627 0.1961];
cEDML = [0.8902 0.6039 0.9608];

% find indices for the drilling sites in temperature 
val_lon = zeros(length(table2array(coord)));
val_lat = zeros(length(table2array(coord)));
idx_lon = zeros(length(table2array(coord)));
idx_lat = zeros(length(table2array(coord)));

for i = 1:length(table2array(coord))
    [val_lon(i),idx_lon(i)] = min(abs(table2array(coord(i,2))-lon));
    [val_lat(i),idx_lat(i)] = min(abs(table2array(coord(i,3))-lat));
end 

% TEMPERATURE FOR THE DIFFERENT DRILLING SITES FOR THE WHOLE RUN  [°C]
T_NGRIP = double(squeeze(T(idx_lon(1),idx_lat(1),:)-273.15));  % SQUEEZE extract only the temperarute as function of time (i.e. 788x1 single)
T_GRIP = double(squeeze(T(idx_lon(2),idx_lat(2),:)-273.15));   % DOUBLE turns single into double
T_NEEM = double(squeeze(T(idx_lon(3),idx_lat(3),:)-273.15)); 
T_GISP = double(squeeze(T(idx_lon(4),idx_lat(4),:)-273.15));
T_WDC = double(squeeze(T(idx_lon(5),idx_lat(5),:)-273.15)); 
T_TAL = double(squeeze(T(idx_lon(6),idx_lat(6),:)-273.15)); 
T_DF = double(squeeze(T(idx_lon(7),idx_lat(7),:)-273.15)); 
T_EDC = double(squeeze(T(idx_lon(8),idx_lat(8),:)-273.15)); 
T_EDML = double(squeeze(T(idx_lon(9),idx_lat(9),:)-273.15)); 

fac = 1.8;

% TEMPERATURE PROFILES NORMALISED: change xlims when zooming into one
% single event
figure('units','centimeter','position',[0 2 40 15], 'Name', 'Model Only Normalised');
%figure
plot(time,(T_NGRIP(:)-min(T_NGRIP))/(max(T_NGRIP)-min(T_NGRIP))*fac,'Color',[cNGRIP 0.2],'HandleVisibility','off'); hold on
plot(time,smoothdata((T_NGRIP(:)-min(T_NGRIP))/(max(T_NGRIP)-min(T_NGRIP)),'movmean',10)*fac,'Color',cNGRIP)
plot(time,(T_GRIP(:)-min(T_GRIP))/(max(T_GRIP)-min(T_GRIP))*fac,'Color',[cGRIP 0.2],'HandleVisibility','off')
plot(time,smoothdata((T_GRIP(:)-min(T_GRIP))/(max(T_GRIP)-min(T_GRIP)),'movmean',10)*fac,'Color',cGRIP)
plot(time,(T_NEEM(:)-min(T_NEEM))/(max(T_NEEM)-min(T_NEEM))*fac,'Color',[cNEEM 0.2],'HandleVisibility','off')
plot(time,smoothdata((T_NEEM(:)-min(T_NEEM))/(max(T_NEEM)-min(T_NEEM))*fac,'movmean',10),'Color',cNEEM)
plot(time,(T_GISP(:)-min(T_GISP))/(max(T_GISP)-min(T_GISP))*fac,'Color',[cGISP 0.2],'HandleVisibility','off')
plot(time,smoothdata((T_GISP(:)-min(T_GISP))/(max(T_GISP)-min(T_GISP))*fac,'movmean',10),'Color',cGISP)
plot(time,(T_WDC(:)-min(T_WDC))/(max(T_WDC)-min(T_WDC))+0.3,'Color',[cWDC 0.2],'HandleVisibility','off')
plot(time,smoothdata((T_WDC(:)-min(T_WDC))/(max(T_WDC)-min(T_WDC))+0.3,'movmean',10),'Color',cWDC)
plot(time,(T_TAL(:)-min(T_TAL))/(max(T_TAL)-min(T_TAL))+0.3,'Color',[cTAL 0.2],'HandleVisibility','off')
plot(time,smoothdata((T_TAL(:)-min(T_TAL))/(max(T_TAL)-min(T_TAL))+0.3,'movmean',10),'Color',cTAL)
plot(time,(T_DF(:)-min(T_DF))/(max(T_DF)-min(T_DF))+0.3,'Color',[cDF 0.2],'HandleVisibility','off')
plot(time,smoothdata((T_DF(:)-min(T_DF))/(max(T_DF)-min(T_DF))+0.3,'movmean',10),'Color',cDF)
plot(time,(T_EDC(:)-min(T_EDC))/(max(T_EDC)-min(T_EDC))+0.3,'Color',[cEDC 0.2],'HandleVisibility','off')
plot(time,smoothdata((T_EDC(:)-min(T_EDC))/(max(T_EDC)-min(T_EDC))+0.3,'movmean',10),'Color',cEDC)
plot(time,(T_EDML(:)-min(T_EDML))/(max(T_EDML)-min(T_EDML))+0.3,'Color',[cEDML 0.2],'HandleVisibility','off')
plot(time,smoothdata((T_EDML(:)-min(T_EDML))/(max(T_EDML)-min(T_EDML))+0.3,'movmean',10),'Color',cEDML)
set(gca,'ytick',[])
ax = gca; 
ax.FontSize = 8;
xlim([time(1) time(end)])
%xlim([xlims(1) xlims(2)]) % zoom into one event.
%xlim([4000 6000]) % zoom into one event.
xlabel('Time [Years]','FontSize',10)
%title(CO2)
legend('NGRIP', 'GRIP', 'NEEM', 'GISP2', 'WDC', 'TAL', 'DF', 'EDC', 'EDML','Location','northeastoutside')
set(gca,'XMinorTick','on')
set(gca,'FontName','Times New Roman')
print(gcf,'C:\Users\annak\OneDrive - University of Copenhagen\Desktop\Speciale\Figures\model210ppm.png','-dpng','-r600')


%%
% 1. FIND TRANSITION
sm = 5;

% finds steepest point for one single event
k1 = find(time >= xlims(1),1);
k2 = find(time >= xlims(2),1);
%timeI = time(k1:k2);
T_sm = smoothdata(T_NGRIP(k1:k2),'movmean',10);

% find differences over larger interval: returns indices and depth of transitions
L = 5;
diff_T = zeros(length(T_sm),1);
for j = L+1:length(T_sm)-L
    diff_T(j) = T_sm(j+L)-T_sm(j-L);
end 
[~,idx_transI] = maxk(diff_T,1); % finds index of transition in interval k1:k2
%time_transI = timeI(idx_transI);
idx_trans = idx_transI + (k1-1); % index of transition in whole set
time_trans = time(idx_trans); 

% 2. STACK THE CORES: normalise Antarctica, stack, smooth
n_AA = 5; 
n_GL = 3; % coarse grid resolution: GRIP = GISP2

T_GL_stack = (T_NGRIP + T_GRIP + T_NEEM)/n_GL;
T_AA_stack = (normalize(T_WDC,'range') + normalize(T_TAL,'range') + normalize(T_DF,'range') + normalize(T_EDC,'range') + normalize(T_EDML,'range'))/n_AA;
T_AA_stack_sm = smoothdata(T_AA_stack,'movmean',sm);


[max_AA, idxmaxAA] = max(T_AA_stack_sm(idx_trans:idx_trans+30));
idxmaxAA = idxmaxAA+(idx_trans-1);
time_bp = time(idxmaxAA);
bp = time_bp-time_trans;

% CROP DATA
[~,k1] = min(abs(time-xlims(1)));
[~,k2] = min(abs(time-xlims(2)));

b = 0.2;
figure('Name','MODEL COMPARISON NEW_02-10-22; normalise - stack - smooth')
a = fill([time_trans, time_trans, time(idxmaxAA), time(idxmaxAA)],[-80, -30, -30, -80], [0.9608 0.9020 0.7020], 'HandleVisibility','off'); hold on 
a.EdgeColor = [0.9608 0.9020 0.7020];
colororder([cNGRIP ;0 0 0]); yyaxis right
plot(time,normalize(T_WDC,'range')*b,'-','Color',[cWDC 0.2],'HandleVisibility','off'); hold on
plot(time,normalize(T_TAL,'range')*b,'-','Color',[cTAL 0.2],'HandleVisibility','off')
plot(time,normalize(T_DF,'range')*b,'-','Color',[cDF 0.2], 'HandleVisibility','off')
plot(time,normalize(T_EDC,'range')*b,'-','Color',[cEDC 0.2],'HandleVisibility','off')
plot(time,normalize(T_EDML,'range')*b,'-','Color',[cEDML 0.2], 'HandleVisibility','off')
plot(time,smoothdata(normalize(T_WDC,'range')*b,'movmean',sm),'-','Color',cWDC)
plot(time,smoothdata(normalize(T_TAL,'range')*b,'movmean',sm),'-','Color',cTAL)
plot(time,smoothdata(normalize(T_DF,'range')*b,'movmean',sm),'-','Color',cDF)
plot(time,smoothdata(normalize(T_EDC,'range')*b,'movmean',sm),'-','Color',cEDC)
plot(time,smoothdata(normalize(T_EDML,'range')*b,'movmean',sm),'-','Color',cEDML)
plot(time, T_AA_stack_sm*b, 'k-','LineWidth',.8)
plot(time(idxmaxAA), max_AA*0.2, 'o','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',4,'HandleVisibility','off'); hold off
%text(time(idxmaxAA)+2,0.23,sprintf('t = %.0f years',bp),'FontName','Times New Roman','FontSize',8) 
axis([xlims(1) xlims(2) 0 0.25])
ax = gca;
ax.YAxis(1).Color = cNGRIP;
ax.YAxis(2).Color = 'k';
ax.FontSize = 8;
set(gca,'XMinorTick','on','YMinorTick','on')
ylabel('Antarctica Normalised Temperature','FontSize',10)
legend('WDC', 'TAL', 'DF', 'EDC', 'EDML','mean','Location','northeastoutside')
 
yyaxis left 
%plot(time, T_NGRIP,'-','Color',[cNGRIP 0.2]); hold on
%plot(time, smoothdata(T_NGRIP,'movmean',sm),'-','Color',cNGRIP); hold off
plot(time, T_GL_stack,'-','Color', [cNGRIP 0.2],'HandleVisibility','off'); hold on 
plot(time, smoothdata(T_GL_stack,'movmean',sm),'-','Color',cNGRIP,'HandleVisibility','off'); hold off
ylabel(['Greenland Temperature [',char(176),'C]'],'FontSize',10)
xlabel('Time [years]','FontSize',10)
axis([xlims(1) xlims(2) ylims(1) ylims(2)])
%text(4790,-42.5,sprintf('t = %.0f years',bp),'FontName','Times New Roman','FontSize',8) 
text(5800,-42,sprintf('t = %.0f years',bp),'FontName','Times New Roman','FontSize',8) 
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'FontName','Times New Roman')
set(gcf, 'Color', 'w');

%% 

b = 0.2;
figure('Name','MODEL COMPARISON NEW_02-10-22; normalise - stack - smooth')
a = fill([time_trans, time_trans, time(idxmaxAA), time(idxmaxAA)],[-80, -30, -30, -80], [0.9608 0.9020 0.7020], 'HandleVisibility','off'); hold on 
a.EdgeColor = [0.9608 0.9020 0.7020];
colororder([cNGRIP ;0 0 0]); yyaxis right
plot(time(k1:k2),normalize(T_WDC(k1:k2),'range')*b,'-','Color',[cWDC 0.2],'HandleVisibility','off'); hold on
plot(time(k1:k2),normalize(T_TAL(k1:k2),'range')*b,'-','Color',[cTAL 0.2],'HandleVisibility','off')
plot(time(k1:k2),normalize(T_DF(k1:k2),'range')*b,'-','Color',[cDF 0.2], 'HandleVisibility','off')
plot(time(k1:k2),normalize(T_EDC(k1:k2),'range')*b,'-','Color',[cEDC 0.2],'HandleVisibility','off')
plot(time(k1:k2),normalize(T_EDML(k1:k2),'range')*b,'-','Color',[cEDML 0.2], 'HandleVisibility','off')
plot(time(k1:k2),smoothdata(normalize(T_WDC(k1:k2),'range')*b,'movmean',sm),'-','Color',cWDC)
plot(time(k1:k2),smoothdata(normalize(T_TAL(k1:k2),'range')*b,'movmean',sm),'-','Color',cTAL)
plot(time(k1:k2),smoothdata(normalize(T_DF(k1:k2),'range')*b,'movmean',sm),'-','Color',cDF)
plot(time(k1:k2),smoothdata(normalize(T_EDC(k1:k2),'range')*b,'movmean',sm),'-','Color',cEDC)
plot(time(k1:k2),smoothdata(normalize(T_EDML(k1:k2),'range')*b,'movmean',sm),'-','Color',cEDML)
plot(time(k1:k2), T_AA_stack_sm(k1:k2)*b, 'k-','LineWidth',.8)
plot(time(idxmaxAA), max_AA*0.2, 'o','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',4,'HandleVisibility','off'); hold off
%text(time(idxmaxAA)+2,0.23,sprintf('t = %.0f years',bp),'FontName','Times New Roman','FontSize',8) 
axis([xlims(1) xlims(2) 0 0.25])
ax = gca;
ax.YAxis(1).Color = cNGRIP;
ax.YAxis(2).Color = 'k';
ax.FontSize = 8;
set(gca,'XMinorTick','on','YMinorTick','on')
ylabel('Antarctica Normalised Temperature','FontSize',10)
legend('WDC', 'TAL', 'DF', 'EDC', 'EDML','mean','Location','northeastoutside')
 
yyaxis left 
%plot(time, T_NGRIP,'-','Color',[cNGRIP 0.2]); hold on
%plot(time, smoothdata(T_NGRIP,'movmean',sm),'-','Color',cNGRIP); hold off
plot(time(k1:k2), T_GL_stack(k1:k2),'-','Color', [cNGRIP 0.2],'HandleVisibility','off'); hold on 
plot(time(k1:k2), smoothdata(T_GL_stack(k1:k2),'movmean',sm),'-','Color',cNGRIP,'HandleVisibility','off'); hold off
ylabel(['Greenland Temperature [',char(176),'C]'],'FontSize',10)
xlabel('Time [years]','FontSize',10)
%axis([xlims(1) xlims(2) ylims(1) ylims(2)])
axis([9 2600 ylims(1) ylims(2)])
%xticks(xlims(1):200:xlims(2))
xticks(0:300:2600)
text(time(idxmaxAA)+50,-40,sprintf('t = %.0f years',bp),'FontName','Times New Roman','FontSize',8) 
text(xlims(2)-450,ylims(2)-0.8,CO2,'FontName','Times New Roman','FontSize',10)
set(gca,'XMinorTick','on','YMinorTick','on')
set(gca,'FontName','Times New Roman')
set(gcf, 'Color', 'w');
