% DETERMINE ONSETS OF DO EVENTS
% Current script to find transitions, compare results to Buizert et al. and
% INTIMATE, Capron et al. 2020 and Myrvoll-Nielsen et al. 2022 as well as stacking of transitions. Uses a NEW DATASET including the entire core.
% uses myfindtransitions_single.m

clear variables; close all; clc
set(0,'DefaultTextInterpreter','tex')
set(0,'DefaultFigureVisible','on')
warning('off','all')

cNGRIP = [0 0.4470 0.7410];

%% DATA: 

% NGRIP5cm.xls including depth and d18O measurements
% excluding the Holocene, starting at 10ka b2k (corresponding to a depth of 1384.34m according to the GICC05 timescale)
NGRIP = readtable('NGRIP5cm.xls');
depth = NGRIP.Var1(6651:end); 
d18O = NGRIP.Var2(6651:end);

% Warming Events: depth and age after INTIMATE and Buizert et al. 
DOs = readtable('Warming_Events.xls');
DO_depth_Intimate = DOs.Depth_m_(1:24); % only Part 1: Holocene to GI-18 
DO_depth_B = DOs.Depth_m__1;
DO_depth_B = DO_depth_B(~isnan(DO_depth_B));
n_B = length(DO_depth_B);

% All warming events without GI-2.1 as it is not found by Buizert et al. up to 60ka b2k 
eventsI = DOs.Var1([1:2,4:24]); % Part 1: Holocene to GI-18
eventsII = DOs.Var1(25:end); % Part 2: GI-19.1 to GI-25c

% BIPOLAR VOLCANOES
volc = readtable('bipolar_volc.xlsx');


%% PART I: FIND TRANSITION FOR EACH EVENT (Second half of last ice age)   

% Limits for Part 1: Holocene to GI-18
depth_lim = [1470 1550 1790 1840 1880 1910 1935 1965 1990 2050 2090 2110 2145 2210 2245 2325 2348 2358 2385 2400 2406 2418 2460; ...
             1495 1650 1795 1870 1900 1950 1960 1978 2020 2080 2110 2135 2165 2240 2265 2351 2361 2375 2400 2408 2418 2425 2475]; % define upper and lower depth limit to zoom into
smI = [20 10 15 10 10 20 15*ones(1,3) 10 10 20 10*ones(1,4) 15 25 5 10*ones(1,3) 20].'; % different smoothing window for each event

% Using my function to find trantions, returns the depth
[depth_transI, diff_d18OI] = myfindtransitions_single(depth_lim, depth, d18O, smI, n_B, DO_depth_B, eventsI); % compare to Buizert et al. (n_B, DO_depth_B)

%% PART II: first half of last ice age (60ka b2k to 120ka b2k)
clc;

DO_depth_Intimate = DOs.Depth_m_(25:end); % Interstadials according to INTIMATE for comparison
n_Intimate = numel(DO_depth_Intimate);

depth_lim = [2500 2530 2575 2680 2690 2745 2875 2895 2910 2930 3000;  % interval incl. the different events
             2510 2545 2585 2690 2695 2755 2895 2900 2930 2950 3005]; % for Part 2: GI-19.1 to 25c           
smII = [15 20 15 15 10*ones(1,3) 20 15*ones(1,2) 10]; % array with smoothing windows

%FIND TRANSITION FOR EACH EVENT 
[depth_transII, diff_d18OII] = myfindtransitions_single(depth_lim, depth, d18O, smII, n_Intimate, DO_depth_Intimate, eventsII);

%% Plot all the transitions onto d18O signal
% creates Figure 6 - Defining onsets of DO events

idx1 = find(depth >= 1945,1);
idx2 = find(depth >= 2020,1);

depthI = depth(idx1:idx2);
%d18Osm = smoothdata(d18O,'movmean',sm(i)); % smoothed data
d18OI_sm = d18Osm(idx1:idx2);
L = 5;

% find differences over larger interval: returns indices and depth of transitions
diff_d18O = zeros(length(d18OI_sm),1);
for j = L+1:length(d18OI_sm)-L
    diff_d18O(j) = d18OI_sm(j+L)-d18OI_sm(j-L);
end 


figure()
T = tiledlayout(5,1,"TileSpacing","tight");
ax1 = nexttile([4 1]);
plot(depth(idx1:idx2), d18O(idx1:idx2), 'Color',[cNGRIP 0.2]); hold on
plot(depth(idx1:idx2), d18Osm(idx1:idx2),'Color',cNGRIP)
xline(depth_trans(:),'Color',[0.5 0.5 0.5])
scatter(depth_trans,d18Osm(idx),20,'o','MarkerFaceColor','w','MarkerEdgeColor','k')
set(gca,'XColor','none')
set(gca,'box','off')
ylabel(['\delta^{18}O [',char(8240),']'])
axis([depth(idx1) depth(idx2) -48 -36])
set(gca,'FontName','Times New Roman')

ax2 = nexttile([1 1]);
for i = [7,8,9]
    a=fill([depth_limI(1,i), depth_limI(1,i), depth_limI(2,i), depth_limI(2,i)],[-3, 3, 3, -3], [ 0.9020 0.9020 0.9020]); hold on
    a.EdgeColor = [0.9020 0.9020 0.9020];
end 
plot(depthI, -diff_d18O,'k')
set(gca,'box','off')
ylabel(['\Delta \delta^{18}O [',char(8240),']'])
xlabel('NGRIP Depth [m]')
xline(depth_trans(:),'Color',[0.5 0.5 0.5])
axis([depth(idx1) depth(idx2) -3 3])
set(gca,'FontName','Times New Roman')

%% INTERPOLATION:

% interpolate on GICC05 time scale with 20-year resolution
GICC05 = readtable('2010-11-19 GICC05modelext for NGRIP.xls');
GICC05_age = GICC05.x_Age;
i10k = find(GICC05_age == 10000); % start right in the beginning of Holocene
GICC05_age = GICC05_age(i10k:end);
GICC05_depth = GICC05.NGRIP1_2(i10k:end);
GICC05_d18O = GICC05.NGRIP(i10k:end);

age = interp1(GICC05_depth,GICC05_age,depth-0.025); % finds age corresponding to the depth and d18O from NGRIP5cm.xls, substract 2.5cm to interpolate onto midpoint

% NGRIP 5cm resolution with GICC05 age, depth (m) and d18O (permil) 
NGRIP = [depth age d18O]; 

% FIND AGE OF TRANSITIONS
% since the indices are known for the midpoints of the transitions, their
% corresponding age can be found. 
depth_trans = [depth_transI; depth_transII];
idx_trans = find(ismember(depth,depth_trans));
age_trans = age(idx_trans); % age of transition midpoint (GICC05)
d18O_trans = d18O(idx_trans);

NGRIP_trans = [depth_trans age_trans d18O_trans];
%save('idx_trans.csv','idx_trans','-ascii')

%% Save Events: include depth, d18O and age of events in excel file

title = {'NGRIP depth (m)','NGRIP d18O (permil)','GICC05 age (years b2k)'};
eventsall = {'Greenland Warming Event','YD-PB', 'GI-1e', 'GI-2.2', 'GI-3', 'GI-4','GI-5.1','GI-5.2','GI-6','GI-7c','GI-8c',...
'GI-9','GI-10','GI-11','GI-12','GI-13c','GI-14e','GI-15.1','GI-15.2','GI-16.1','GI-16.2', 'GI-17.1c','GI-17.2' ...
'GI-18', 'GI-19.1', 'GI-19.2', 'GI-20c', ' GI-21.1e', 'GI-21.2','GI-22g','GI-23.1', ...
           'GI-23.2','GI-24.1c','GI-24.2','GI-25c'};
depth_trans_cell = [title(1)' num2cell(depth_trans)']';
d18O_trans_cell = [title(2)' num2cell(d18O_trans)']';
age_trans_cell = [title(3)' num2cell(age_trans)']';

NGRIP_trans_cell = [eventsall' depth_trans_cell d18O_trans_cell age_trans_cell];
%writecell(NGRIP_trans_cell,'NGRIP_transitions.xls');

%% COMPARE MY RESULTS TO OTHER STUDIES
% GI-2.1 removed (only found in Rasmussen et al., 2014)
% Produces FIG 4.2 and FIG 4.3

x = 1:numel(NGRIP_trans_cell(([1:2,4:end]),1));
DOs.mft_depth = [NGRIP_trans(1:2,1); NaN; NGRIP_trans(3:end,1)];
DOs.mft_age = [NGRIP_trans(1:2,2); NaN; NGRIP_trans(3:end,2)];

mft_depth = DOs.mft_depth([1:2,4:end]);
mft_age = DOs.mft_age([1:2,4:end]);

R_depth = DOs.Depth_m_([1:2,4:end]); %Rasmussen et al. (2014)
R_vs_mft_depth = R_depth-mft_depth;
R_age = DOs.Age_aB2k_([1:2,4:end]);
R_vs_mft_age = R_age-mft_age; % Rasmussen et al. 

B_depth = DOs.Depth_m__1([1:2,4:end]); %Buizert et al. (2015)
B_vs_mft_depth = B_depth-mft_depth;
B_age = DOs.Age_aB2k__1([1:2,4:end]);
B_vs_mft_age = B_age-mft_age;

C_depth = DOs.Depth_C_m_([1:2,4:end]); %Capron et al. (2021)
C_vs_mft_depth = C_depth-mft_depth;
C_age = DOs.Age_C_aB2k_([1:2,4:end]);
C_vs_mft_age = C_age-mft_age;

MN_depth = DOs.Depth_MN_m_([1:2,4:end]); %Myrvoll-Nielsen et al. (2022)
MN_vs_mft_depth = MN_depth-mft_depth;
MN_age = DOs.Age_MV_aB2k_([1:2,4:end]);
MN_vs_mft_age = MN_age-mft_age;

msz = 5;
figure("Units","centimeters",'Position',[1 2 40 15]);
tiledlayout(1,2) %,"TileSpacing","compact")
nexttile
plot(x,R_vs_mft_depth,'o','MarkerEdgeColor','[0.4660 0.6740 0.1880]','MarkerFaceColor','[0.4660 0.6740 0.1880]','MarkerSize',msz); hold on
plot(x,B_vs_mft_depth,'square','MarkerEdgeColor','[0 0.4470 0.7410]','MarkerFaceColor','[0 0.4470 0.7410]','MarkerSize',msz)
plot(x, C_vs_mft_depth,'pentagram','MarkerEdgeColor','[0.9290 0.6940 0.1250]','MarkerFaceColor','[0.9290 0.6940 0.1250]','MarkerSize',msz)
plot(x, MN_vs_mft_depth,'diamond','MarkerEdgeColor','[0.4940 0.1840 0.5560]','MarkerFaceColor','[0.4940 0.1840 0.5560]','MarkerSize',msz)
hold off
xlim([0 numel(x)+1])
xticks(1:numel(x))
labels = {'YD-PB', 'OD-BA', '2.2', '3', '4','5.1','5.2','6','7c','8c',...
    '9','10','11','12','13c','14e','15.1','15.2','16.1','16.2', '17.1c','17.2',...
    '18', '19.1', '19.2', '20c', ' 21.1e', '21.2','22g','23.1','23.2','24.1c','24.2','25c'};
set(gca, 'XTickLabel',labels, 'FontName', 'Times New Roman', 'FontSize',8)
ylabel('\DeltaDepth [m]','FontName','Times New Roman', 'FontSize',10)
set(gca,'FontName','Times New Roman')
set(gca, 'YMinorTick','on')
grid on
set(gca,'YMinorGrid','on','MinorGridAlpha',0.1,'MinorGridLineStyle','-')
%legend('Rasmussen et al. (2014)','Buizert et al. (2015)', 'Capron et al. (2021)','Myrvoll-Nielsen et al. (2022)')

nexttile
plot(x,R_vs_mft_age,'o','MarkerEdgeColor','[0.4660 0.6740 0.1880]','MarkerFaceColor','[0.4660 0.6740 0.1880]','MarkerSize',msz); hold on
plot(x,B_vs_mft_age,'square','MarkerEdgeColor','[0 0.4470 0.7410]','MarkerFaceColor','[0 0.4470 0.7410]','MarkerSize',msz)
plot(x, C_vs_mft_age,'pentagram','MarkerEdgeColor','[0.9290 0.6940 0.1250]','MarkerFaceColor','[0.9290 0.6940 0.1250]','MarkerSize',msz)
plot(x, MN_vs_mft_age,'diamond','MarkerEdgeColor','[0.4940 0.1840 0.5560]','MarkerFaceColor','[0.4940 0.1840 0.5560]','MarkerSize',msz); hold off
legend('Rasmussen et al. (2014)','Buizert et al. (2015)', 'Capron et al. (2021)','Myrvoll-Nielsen et al. (2022)', 'Location','northeast','FontName','Times New Roman')
xlim([0 numel(x)+1])
xticks(1:numel(x))
ylim([-100 150])
yticks(-150:30:150)
labels = {'YD-PB', 'OD-BA', '2.2', '3', '4','5.1','5.2','6','7c','8c',...
    '9','10','11','12','13c','14e','15.1','15.2','16.1','16.2', '17.1c','17.2',...
    '18', '19.1', '19.2', '20c', ' 21.1e', '21.2','22g','23.1','23.2','24.1c','24.2','25c'};
set(gca, 'XTickLabel',labels, 'FontName', 'Times New Roman', 'FontSize',8)
ylabel('\DeltaAge [years]','FontName', 'Times New Roman','FontSize',10)
set(gca,'FontName','Times New Roman')
grid on
set(gca,'YMinorTick','on','YMinorGrid','on','MinorGridAlpha',0.1,'MinorGridLineStyle','-')
%print(gcf,'C:\Users\annak\OneDrive - University of Copenhagen\Desktop\Speciale\Figures\trans_diff.png','-dpng','-r600')


events = eventsall(2:end);
cNGRIP = [0 0.4470 0.7410];
n = 40;

B_age(isnan(B_age)) = 0;
C_age(isnan(C_age)) = 0;
MN_age(isnan(MN_age)) = 0;

figure('PaperUnits', 'centimeters', 'PaperPosition', [0 0 16 17]);
fig = tiledlayout(7,5,'TileSpacing','tight');

for i = 1:numel(NGRIP_trans_cell(([1:2,4:end]),1))
    nexttile
    plot((NGRIP(idx_trans(i)-n:idx_trans(i)+n,2))/1000, NGRIP(idx_trans(i)-n:idx_trans(i)+n,3),'Color',[cNGRIP 0.2],'HandleVisibility','off'); hold on
    plot((NGRIP(idx_trans(i)-n:idx_trans(i)+n,2))/1000, smoothdata(NGRIP(idx_trans(i)-n:idx_trans(i)+n,3),'movmean',5),'Color',cNGRIP,'HandleVisibility','off')
    xline(mft_age(i)/1000,'Color', cAMK)
    xline(R_age(i)/1000,'Color',cRasm)
    xline(B_age(i)/1000,'Color',cBuiz)
    xline(C_age(i)/1000,'Color',cCapr)
    xline(MN_age(i)/1000,'Color',cMVN)
    xlim([(NGRIP(idx_trans(i)-n,2))/1000 (NGRIP(idx_trans(i)+n,2))/1000])
    title(events(i))
    set(gca,'FontSize',5,'FontName','Times New Roman')
    
end 
ylabel(fig, ['\delta^{18}O [',char(8240),']'],'FontName','Times New Roman','FontSize',10);
xlabel(fig, 'GICC05modelext Age [ka b2k]','FontName','Times New Roman','FontSize',10);
lg = legend('This Work','Rasmussen et al. (2014)','Buizert et al. (2015)', 'Capron et al. (2021)','Myrvoll-Nielsen et al. (2022)','FontSize',4); %, 'Position',[10 10 10 10])
lg.Layout.Tile = 35;
lg.ItemTokenSize = [5,18];
%print(gcf,'C:\Users\annak\OneDrive - University of Copenhagen\Desktop\Speciale\Figures\comp_alltrans.png','-dpng','-r600')

