%% BIPOLAR PHASING 
% Analysis of the climate signal across the transitions associated with Dansgaard-Oeschger events by stacking the events aligned at the onset and determining the Antarctic response time for all events, those occurring in the early and the late ice age, respectively, and minor as well as major events. 
% This script uses the following functions:
% - mystacking.m to stack the events for the individual ice core records
% - myclasscom.m to compare stacks of the different classifications
% This script uses the following data:
% -NGRIP: 'NGRIP5cm.xls'
% - GRIP: 'GRIP_d18O.xls'
% - GISP2: 'gisp2_measured.xls'
% - NEEM: 'NEEM_d18O_5cm.xls'
% - WAIS Divide: 'WAIS_project_members_Source_Data.xls'
% - Talos Dome: 'TD_d18O.xls'
% - DOME Fuji:  'DF_d18O.xls'
% - EDC: 'dc5pbag_1500m.xls' and 'EDC96-EDC99_conversion.xls' for conversion
% - EDML: 'EDML_d18O_NEW.xls'
% - TRANSITIONS: 'Onset_depths_all_cores_ver_2022-06-12.xls'
% - BIPOLAR MATCHPOINTS: 'Bipolar match points ver 2022-06-13.xls'
% - GICC05: '2010-11-19 GICC05modelext for NGRIP.xls'
% - BIPOLAR VOLCANOES: 'bipolar_volc.xlsx'
%%
clear variables; close all; clc
warning('off','all')
set(0,'DefaultFigureVisible','on');
set(0,'DefaultAxesFontName','Times New Roman','DefaultTextFontName','Times New Roman');

% SAVE FIGURES BY SETTING SAVE TO 1 + CHANGE FOLDER NAME (last section)!!
save = 0; % set to 2 to save stacks, set to 3 to save d18Otimed + stacks for comparison


CORES = {'NGRIP', 'GRIP', 'NEEM', 'GISP2', 'WDC', 'TAL', 'DOME F', 'EDC', 'EDML'};

% DEFINE COLOUR CODES
cNGRIP = [0 0.4470 0.7410];
cGRIP = [0.7350 0.0780 0.1840];
cNEEM = [0.9290 0.6940 0.1250];
cGISP = [0.4660 0.6740 0.1880];
cWDC = [0.4940 0.1840 0.5560];
cTAL = [0.2 0.6 0.4];
cDF = [0.6235 0.8235 0.9882];
cEDC = [0.8902 0.4627 0.1961];
cEDML = [0.8902 0.6039 0.9608];
% DEFINE ALPHA: a = 0.2;

cfirst = [ 0.9294 0.6941 0.1255];
csecond = [0.4902 0.2039 0.1412];
cminor = [0.7686 0.5333 0.9882];
cmajor = [0.3412 0.0157 0.3255];

% load warming events and volcanic matchpoints
transitions = readtable('Onset_depths_all_cores_ver_2022-06-12.xls');
events = transitions.DO_onset_;
age_trans = transitions.GICC05_age_b2k_;
matchpoints = readtable('Bipolar match points ver 2022-06-13.xls','Format','auto'); 
d_bp_match = transitions.Match_uncertainty_yr_;  % matching uncertainty for each transition


% chose DO Events to focus on
allDO = 1:numel(age_trans); % all events (YD-PB to GI-24.2)
secondhalf = 1:23;  % more recent half of the last ice age (second half), YD-PB to GI-18
firsthalf = 24:numel(age_trans); % first half of last ice age (GI-19 to GI-24.2)
majorDO = [2 3 5 10 14 16 22 25 26 27 29 30]; % DOs following a Heinrich Event (H-1 to H-10)
minorDO = setdiff(allDO, majorDO); 
ofinterest = allDO; % CHANGE THIS TO INVESTIGATE DIFFERENT 'CLASSES'

%% STACKING OF EVENTS FOR INDIVIDUAL ICE CORES
% uses mystacking.m: plots the stacked records across the events along a time vector; gives the change in d18O across the transition
% Overwrites the input for all ice cores besides NGRIP.

% NGRIP

% !!! do not overwrite matchpoints, needed to set Antarctic ice cores on GICC05 timescale
mp(:,1) = vertcat(matchpoints.NGRIP_depth_m_, transitions.NGRIP); % NGRIP matchpoints including transitions 

% load NGRIP data
NGRIP = readtable('NGRIP5cm.xls');
NGRIPdepth = NGRIP.Var1-0.025; % bottom data, substract 2.5 cm to get midpoint of sample
d18O = NGRIP.Var2;

% load age vs depth data on GICC05 timescale (later used for interpolation)
GICC05 = readtable('2010-11-19 GICC05modelext for NGRIP.xls');
GICC05_age = GICC05.x_Age;
GICC05_depth = GICC05.NGRIP1_2;

[NGRIPage, NGRIPdepth, NGRIPstack, NGRIPd18Otimed, NGRIPdeld18O] = mystacking(GICC05_age, GICC05_depth, NGRIPdepth, d18O, events, age_trans, CORES(1), ofinterest, cNGRIP);


% GRIP

Core = readtable('GRIP_d18O.xls');
depth = Core.Depth_m_(2:end)-0.0125; % bottom data, substract 1.25cm to get midpoint of sample
d18O = Core.d18O_permil_(2:end); % remove first value as there is no isotope measurement

mp(:,2) = vertcat(matchpoints.GRIP_depth_m_, transitions.GRIP); % all matchpoints including transitions
idx = find(~isnan(mp(:,1)) & ~isnan(mp(:,2)));

% set on NGRIP depth
NGRIPdepth_CORE = interp1(unique(mp(idx,2)), unique(mp(idx,1)), depth); % NGRIP depth scale for GRIP

[~, ~, GRIPstack, GRIPd18Otimed, GRIPdeld18O] = mystacking(NGRIPage, NGRIPdepth, NGRIPdepth_CORE, d18O, events, age_trans, CORES(2), ofinterest, cGRIP);


% GISP2

% 6.5cm-resolution data; DUPLICATED DATA & MISSING DATA in dataset
Core = readtable('gisp2_measured.xls'); 
depth = Core.Depth-0.0325; % bottom data, substract 3.25cm to get midpoint of sample
d18O = Core.Del18; 
mp(:,2) = vertcat(matchpoints.GISP2_depth_m_, transitions.GISP2); % all matchpoints including transitions
idx = find(~isnan(mp(:,1)) & ~isnan(mp(:,2)));

% remove missing data from dataset (marked as 999999)
idx_99 = find(d18O == 999999);
d18O(d18O == 999999) = [];
depth(idx_99) = [];

% find duplicates in depth and remove them, take average of their corresponding d18O values
[v, w] = unique(depth, 'stable' );
dup_idx = setdiff(1:numel(depth), w); % finds index of the duplicate (second value) 
depth(dup_idx) = [];
avg18O = zeros(length(dup_idx),1);
for i = 1:length(dup_idx)
    avg18O(i) = (d18O(dup_idx(i))+d18O(dup_idx(i)-1))/2;
    d18O(dup_idx(i)-1) = avg18O(i);
end
d18O(dup_idx) = [];

% set on NGRIP depth
NGRIPdepth_CORE = interp1(unique(mp(idx,2)), unique(mp(idx,1)), depth); % NGRIP depth scale for CORE

[~, ~, GISPstack, GISPd18Otimed, GISPdeld18O] = mystacking(NGRIPage, NGRIPdepth, NGRIPdepth_CORE, d18O, events, age_trans, CORES(4), ofinterest, cGISP);


% NEEM 

Core = readtable('NEEM_d18O_5cm.xls');
depth = Core.Depth-0.025; % bottom data, substract 2.5cm to get midpoint of sample
d18O = Core.d18; % remove first value as there is no isotope measurement

mp(:,2) = vertcat(matchpoints.NEEM_depth_m_, transitions.NEEM); % all matchpoints including transitions
idx = find(~isnan(mp(:,1)) & ~isnan(mp(:,2)));

% set on NGRIP depth
NGRIPdepth_CORE = interp1(unique(mp(idx,2)), unique(mp(idx,1)), depth); % NGRIP depth scale for CORE

[~, ~, NEEMstack, NEEMd18Otimed, NEEMdeld18O] = mystacking(NGRIPage, NGRIPdepth, NGRIPdepth_CORE, d18O, events, age_trans, CORES(3), ofinterest, cNEEM);



% WAIS DIVIDE (WDC)

Core = readtable('WAIS_project_members_Source_Data.xls','Sheet','d18O');
d18O = Core.IsotopeData;
depth = (Core.Depths + Core.Var3)/2;

mp(:,2) = vertcat(matchpoints.WD_depth_m_, transitions.WDC); % all matchpoints including transitions
idx = find(~isnan(mp(:,1)) & ~isnan(mp(:,2)));

% Antarctic ice cores on GICC05 timescale
NGRIPdepth_CORE = interp1(unique(mp(idx,2)), unique(mp(idx,1)), depth); % NGRIP depth scale for AA ice core

[WDCage, ~, WDCstack, WDCd18Otimed, WDCdeld18O] = mystacking(NGRIPage, NGRIPdepth, NGRIPdepth_CORE, d18O, events, age_trans, CORES(5), ofinterest, cWDC);

%diffWDC = diff(WDCage);


% TALOS DOME (TAL)

Core = readtable('TD_d18O.xls');
d18O = Core.d18O;
depth = (Core.TopDepth + Core.BottomDepth)/2; 

%mp(:,1) = vertcat(matchpoints.NGRIP_depth_m_, transitions.NGRIP);
mp(:,2) = vertcat(matchpoints.TALDICE_depth_m_, transitions.TAL); % all matchpoints including transitions
idx = find(~isnan(mp(:,1)) & ~isnan(mp(:,2)));

% Antarctic ice cores on GICC05 timescale
NGRIPdepth_CORE = interp1(unique(mp(idx,2)), unique(mp(idx,1)), depth); % NGRIP depth scale for AA ice core

[~, ~, TALstack, TALd18Otimed, TALdeld18O, TALavg_stadial, TALavg_interstadial] = mystacking(NGRIPage, NGRIPdepth, NGRIPdepth_CORE, d18O, events, age_trans, CORES(6), ofinterest, cTAL);



% DOME FUJI (Dome F)

Core = readtable('DF_d18O.xls');
d18O = Core.d18O;
depthDF = (Core.TopDepth + Core.BottomDepth)/2; % mean of top and bottom depth 

mp(:,2) = vertcat(matchpoints.Dome_Fuji_depth_m_, transitions.DF); % all matchpoints including transitions
idx = find(~isnan(mp(:,1)) & ~isnan(mp(:,2)));

% Antarctic ice cores on GICC05 timescale
NGRIPdepth_CORE = interp1(unique(mp(idx,2)), unique(mp(idx,1)), depthDF); % NGRIP depth scale for AA ice core

[~, ~, DFstack, DFd18Otimed, DFdeld18O] = mystacking(NGRIPage, NGRIPdepth, NGRIPdepth_CORE, d18O, events, age_trans, CORES(7), ofinterest, cDF);


% EDC

Core = readtable('dc5pbag_1500m.xls');
depth = Core.Depth-0.025;
d18O = Core.delO18_nopeaks;

% EDC d18O is on the EDC96 depth scale until 762 m depth - needs to be converted.
EDC96_99 = readtable('EDC96-EDC99_conversion.xls');
D96_99 = find(depth>762,1);
EDC99_depth = interp1(EDC96_99.EDC96_depth_m_,EDC96_99.EDC99_depth_m_, depth(1:D96_99));
depth(1:D96_99) = EDC99_depth;

mp(:,2) = vertcat(matchpoints.EDC_depth_m_, transitions.EDC); % all matchpoints including transitions
idx_EDC = find(~isnan(mp(:,1)) & ~isnan(mp(:,2)));

% Antarctic ice cores on GICC05 timescale
NGRIPdepth_CORE = interp1(unique(mp(idx_EDC,2)), unique(mp(idx_EDC,1)), depth); % NGRIP depth scale for AA ice core

[EDCage, ~, EDCstack, EDCd18Otimed, EDCdeld18O] = mystacking(NGRIPage, NGRIPdepth, NGRIPdepth_CORE, d18O, events, age_trans, CORES(8), ofinterest, cEDC);



% EDML 
% 50-cm resolution, measured in top and bottom depth, but already centered

Core = readtable('EDML_d18O_NEW.xls');
depth = Core.DepthIce_snow_m_; % centered depth
d18O = Core.x__18OH2O____SMOW__measured_MassSpectrometerF____;

mp(:,2) = vertcat(matchpoints.EDML_depth_m_, transitions.EDML); % all matchpoints including transitions
idx = find(~isnan(mp(:,1)) & ~isnan(mp(:,2)));

% Antarctic ice cores on GICC05 timescale
NGRIPdepth_CORE = interp1(unique(mp(idx,2)), unique(mp(idx,1)), depth); % NGRIP depth scale for AA ice core

[EDMLage, ~, EDMLstack, EDMLd18Otimed, EDMLdeld18O, ~,~, EDMLstackdeltad18O] = mystacking(NGRIPage, NGRIPdepth, NGRIPdepth_CORE, d18O, events, age_trans, CORES(9), ofinterest, cEDML);

%% DETERMINE BIPOLAR PHASING
% 1. normalise data;  2. stack data; 3. smooth data; 4. determine breakpoint (maximum in Antarctic stack)


t = -1200:1200; % time vector
n_GL = 4; % number Greenlandic ice cores
n_A = 5; % number Antarctic ice cores 
sm = 12; % smoothing window; smoothing procedure: moving average

NGRIPstack_norm = normalize(NGRIPstack,'range');
GRIPstack_norm = normalize(GRIPstack,'range');
GISPstack_norm = normalize(GISPstack,'range');
NEEMstack_norm = normalize(NEEMstack,'range');

GLstack_norm = (NGRIPstack_norm + GRIPstack_norm + GISPstack_norm + NEEMstack_norm)/n_GL;

WDCstack_norm = normalize(WDCstack,'range');
TALstack_norm = normalize(TALstack,'range');
DFstack_norm = normalize(DFstack,'range');
EDCstack_norm = normalize(EDCstack,'range');
EDMLstack_norm = normalize(EDMLstack,'range');

 if isnan(WDCstack_norm(1))
    AAstack_norm = (EDCstack_norm +EDMLstack_norm + DFstack_norm)/3; % NO data for WDC, TAL in first half --> exclude from mean when ofinterest = firsthalf
 else 
    AAstack_norm = (EDCstack_norm + EDMLstack_norm + TALstack_norm + DFstack_norm + WDCstack_norm)/n_A;
 end


 % EXCLUDE DF and EDML
 if isnan(WDCstack_norm(1))
    AAstack_norm = (EDCstack_norm); % NO data for WDC, TAL in first half --> exclude from mean when ofinterest = firsthalf
 else 
    AAstack_norm = (EDCstack_norm + TALstack_norm + WDCstack_norm)/3;
 end

% Determine bipolar phasing (break point)
[smAAmax_norm, idxmaxAA_norm] = max(smoothdata(AAstack_norm,'movmean',sm));


disp('WARNING: uncertainties in plot need to be changed manually!!')

% Produces FIG 4.6 & F.1-F.4
figure('Name','Bipolar Phasing: first stack then smooth')
tiledlayout(2,1,"TileSpacing",'none')
GL = nexttile;
a = fill([0, 0, t(end-idxmaxAA_norm+1), t(end-idxmaxAA_norm+1)],[-3, 3, 3, -3], [0.9608 0.9020 0.7020], 'HandleVisibility','off'); hold on 
a.EdgeColor = [0.9608 0.9020 0.7020];
plot(-t,smoothdata(NGRIPstack_norm,'movmean',sm),'Color',cNGRIP); hold on 
plot(-t,smoothdata(GRIPstack_norm,'movmean',sm),'Color',cGRIP)
plot(-t,smoothdata(NEEMstack_norm,'movmean',sm),'Color',cNEEM)
plot(-t,smoothdata(GISPstack_norm,'movmean',sm),'Color',cGISP)
plot(-t,smoothdata(GLstack_norm,'movmean',sm),'k','LineWidth',.8, 'HandleVisibility','off')
ax = gca; 
ax.FontSize = 8;
ylabel(['Greenlandic \delta^{18}O [',char(8240),']'], 'FontSize',10) 
set(gca,'YMinorTick','on', 'XMinorTick','on')
axis([-400 600 -0.1 1.1])
set(GL, 'XColor','none')
box off
legend('NGRIP','GRIP','NEEM', 'GISP2','mean','Location','northeastoutside')

AA = nexttile; 
a = fill([0, 0, t(end-idxmaxAA_norm+1), t(end-idxmaxAA_norm+1)],[-0.7, 2, 2, -0.7], [0.9608 0.9020 0.7020], 'HandleVisibility','off'); hold on 
a.EdgeColor = [0.9608 0.9020 0.7020];
if isnan(WDCstack(1)) && length(ofinterest) == length(events) % NO DATA FOR WDC AND TAL BEFORE 60ka b2k
    plot(-t,0.6*(smoothdata(DFstack_norm,'movmean',sm))-0.3,'Color',cDF)
    plot(-t,0.6*(smoothdata(EDCstack_norm,'movmean',sm))-0.3,'Color',cEDC)
    plot(-t,(0.6*(smoothdata(EDMLstack_norm,'movmean',sm))-0.3),'Color',cEDML) % -0.09
    plot(-t, 0.6*(smoothdata(AAstack_norm,'movmean',sm))-0.3,'k','LineWidth',0.8)
    legend('DF','EDC','EDML','mean','Location','northeastoutside')
   
else 
    plot(-t,0.6*(smoothdata(WDCstack_norm,'movmean',sm))-0.3,'Color',cWDC); hold on 
    plot(-t,0.6*(smoothdata(TALstack_norm,'movmean',sm))-0.3,'Color', cTAL)
    plot(-t,0.6*(smoothdata(DFstack_norm,'movmean',sm))-0.3,'Color',cDF)
    plot(-t,0.6*(smoothdata(EDCstack_norm,'movmean',sm))-0.3,'Color',cEDC)
    plot(-t,(0.6*(smoothdata(EDMLstack_norm,'movmean',sm))-0.3),'Color',cEDML) % -0.09
    plot(-t, 0.6*(smoothdata(AAstack_norm,'movmean',sm))-0.3,'k','LineWidth',0.8)
    legend('WDC','TAL','DF','EDC','EDML','mean', 'Location','northeastoutside')
end
plot(-t(idxmaxAA_norm), 0.6*smAAmax_norm-0.3, 'o','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',4,'HandleVisibility','off')
text(-t(idxmaxAA_norm),0.6*smAAmax_norm-0.22,sprintf(['t = %.0f ',char(177),' 62 years'],-t(idxmaxAA_norm-1)),'FontSize',8); hold off
ax = gca; 
ax.FontSize = 8;
ylabel(['Antarctic \delta^{18}O [',char(8240),']'],'FontSize',10)  %%     CHANGE UNCERTAINTIES MANUALLY!!!
xlabel('Time [years]', 'FontSize',10)
set(gca,'YMinorTick','on', 'XMinorTick','on')
axis([-400 600 -0.425 0.5])
box off



%% TEST INFLUENCE OF SMOOTHING THE DATA
% Produces Fig. D.1
% To test the influence of smoothing data the stack is smoothed using different windows before determining the breakpoint 
% in the five-core stacks for each individual event. The break points for the individual events are afterwards averaged.

% t = -1200:1200; 
% min_sm = 5;
% max_sm = 120;
% sm = (min_sm:5:max_sm); % smoothing windows
% AAstack = zeros(numel(events),length(EDCd18Otimed)); % preallocate matrix for Antarctic stack of ice cores for different smoothing intervals for each event
% bp = zeros(numel(events), numel(sm)); % preallocate vector for bipolar phasing/delay for each event for each smoothing interval
% mean_bp = zeros(numel(sm),1); % preallocate vector for response time mean
% std_bp = zeros(numel(sm),1); % preallocate vetor for standard deviation
% 
% for i = 1:numel(events)
%     for j = 1:numel(sm)
%         for k = 1:length(EDCd18Otimed)
%              AAstack(1:22,k) = ((WDCd18Otimed(1:22,k)+TALd18Otimed(1:22,k) + DFd18Otimed(1:22,k) + EDCd18Otimed(1:22,k) + EDMLd18Otimed(1:22,k))/5);
%              AAstack(23,k) =  ((WDCd18Otimed(23,k)+ DFd18Otimed(23,k) + EDCd18Otimed(23,k) + EDMLd18Otimed(23,k))/4);
%              AAstack(24:25,k) = ((DFd18Otimed(24:25,k) + EDCd18Otimed(24:25,k) + EDMLd18Otimed(24:25,k))/3);
%              AAstack(26:33,k) = (EDCd18Otimed(26:33,k) + EDMLd18Otimed(26:33,k))/2;
%             [~, idx(i,j)] = max(smoothdata(AAstack(i,1201-300:1201),'movmean',sm(j))); % finds max in interval between t = 0:300 years
%             idx(i,j) = 1201-301+idx(i,j); % finds index of maximum in t vector
%             bp(i,j) = t(end-idx(i,j)); % finds breakpoint (after transition, t=0)
%             std_bp(j) = std(bp(:,j));
%             mean_bp(j) = mean(bp(:,j),'omitnan');
% %           figure(i)
% %           plot(-t, smoothdata(AAstack(i,:),'movmean',10)); hold on
% %           xline(0,'Color',[0.5 0.5 0.5],'LineWidth',2,'Alpha',.4, 'HandleVisibility','off');
% %           xline(bp(i), 'Color', [0.8314 0.4471 0.043], 'LineWidth',2, 'Alpha', .4, 'HandleVisibility', 'off')
% %           xlim([-1200 1200])
% %           xlabel('Time [years]')
% %           ylabel(['\Delta \delta^{18}O [',char(8240),']'])
% %           set(gca, 'FontName', 'Times New Roman')
%         end 
%     end 
% end 
% 
% x2 = [sm, fliplr(sm)];
% spread = [(mean_bp+(std_bp/2))', fliplr(mean_bp-(std_bp/2))']; % change back to std_bp
% figure('units','centimeter','position',[0 1.5 45 15],'Name', 'Smoothing Influence')
% tiledlayout(1,2,"TileSpacing",'compact')
% T_allDOs = nexttile;
% a = fill(x2, spread, [.5 .5 .5],'HandleVisibility','off'); hold on
% set(a,'EdgeColor',[.5 .5 .5],'FaceAlpha',.4,'EdgeAlpha',.4)
% plot(sm, mean_bp,'-', 'Color',[.5 .5 .5 .8],'LineWidth',3,'HandleVisibility','off'); hold on
% for i = 1:numel(events)
%     plot(sm, bp(i,:),'.-','Color', rand(1,3))
% end
% hold off
% ax = gca; 
% ax.FontSize = 8;
% axis([sm(1) sm(end) 0 300])
% xlabel('Smoothing Window [Number of Data Points]', 'FontSize',10)
% ylabel('Response Time [years]','FontSize',10)
% legend('YD-PB', 'GI-1', 'GI-2.2','GI-3','GI-4','GI-5.1','GI-5.2','GI-6','GI-7','GI-8','GI-9','GI-10','GI-11','GI-12','GI-13','GI-14',...
%         'GI-15.1','GI-15.2','GI-16.1','GI-16.2','GI-17.1','GI-17.2','GI-18','GI-19.1','GI-19.2','GI-20','GI-21.1','GI-22.2','GI-23.1','GI-23.1',...
%         'GI-24.1','GI-24.2','Location','Northwestoutside')
% 
% 
% % EXCLUDE SOME EVENTS WITHOUT A TYPICAL SIGNAL
% bp_exc = bp([2:10 13:14 16 18:21 23:26 33],:);
% mean_bp_exc = mean(bp_exc,'omitnan');
% x2 = [sm, fliplr(sm)];
% spread = [(mean_bp+(std_bp/2))', fliplr(mean_bp-(std_bp/2))'];
% 
% T_exc = nexttile;
% a = fill(x2, spread, [.5 .5 .5],'HandleVisibility','off'); hold on
% set(a,'EdgeColor',[.5 .5 .5],'FaceAlpha',.4,'EdgeAlpha',.4)
% plot(sm, mean_bp,'-', 'Color',[.5 .5 .5 .8],'LineWidth',3,'HandleVisibility','off')
% for i = [2:10 13:14 16 18:21 23:26 33]
%     plot(sm, bp(i,:),'.-','Color', rand(1,3))
% end
% hold off
% axis([sm(1) sm(end) 0 300])
% ax = gca; 
% ax.FontSize = 8;
% xlabel('Smoothing Window [Number of Data Points]', 'FontSize',10)
% ylabel('Response Time [years]', 'FontSize',10)
% legend('GI-1', 'GI-2.2','GI-3','GI-4','GI-5.1','GI-5.2','GI-6','GI-7','GI-8','GI-11','GI-12','GI-14',...
%         'GI-15.1','GI-16.1','GI-16.2','GI-17.1','GI-18','GI-19.1','GI-19.2','GI-20','GI-24.2','Location','Northeastoutside')


%% Test smoothing and normalising effect on response time determined in stack of all cores and events
% produces Fig. D.2
% Stacks of all events and ice cores for Greenland and Antarctica
% (unsmoothed, unnormalised)
n_GL = 4; % number Greenlandic ice cores
n_A = 5; % number Antarctic ice cores 
GLstack = (NGRIPstack + GRIPstack + GISPstack + NEEMstack)/n_GL;
if isnan(WDCstack(1)) % NO data for WDC and TAL in first half --> exclude from mean (when ofinterest = firsthalf)
    AAstack = (EDCstack + EDMLstack + DFstack)/3;
    else 
    AAstack = (EDCstack + EDMLstack + TALstack + DFstack + WDCstack)/n_A;
end

idxbp = zeros(length(sm),1); % !!! CHANGE NAME preallocate vector for indeces of break points in stacks for each smoothing window
bp_all = zeros(length(sm),1); %preallocate vector for break points (NON-normalised stack)
bp_all_norm = zeros(length(sm),1); % preallocate vector for break points (normalised stack)

% determine break points for unsmoothed data
[~, idxAA_unsmoothed] = max(AAstack);
bp_unsmoothed = t(end-idxAA_unsmoothed+1);
[~, idxAA_unsmoothed] = max(AAstack_norm);
bp_unsmoothed_norm = t(end-idxAA_unsmoothed+1);

% determine breakpoints for each smoothing window
for j = 1:length(sm)
    [~,idxbp(j)] = max(smoothdata(AAstack,'movmean',sm(j)));
    bp_all(j) = t(end-idxbp(j)+1);
    [~,idxbp(j)] = max(smoothdata(AAstack_norm,'movmean',sm(j)));
    bp_all_norm(j) = t(end-idxbp(j)+1);
end 

bp_all = vertcat(bp_unsmoothed,bp_all);
bp_all_norm = vertcat(bp_unsmoothed_norm,bp_all_norm);

figure('Name','TEST SMOOTHING INFLUENCE')
a = fill([0,0,max_sm,max_sm],[min(bp_all), max(bp_all), max(bp_all), min(bp_all)], [0.7686 0.5333 0.9882],'HandleVisibility','off'); hold on
set(a,'EdgeColor',[0.7686 0.5333 0.9882], 'FaceAlpha', .2, 'EdgeAlpha',.2)
%a = fill([0,0,max_sm, max_sm],[min(bp_all_norm),max(bp_all_norm),max(bp_all_norm),min(bp_all_norm)],[0.3412 0.0157 0.3255], 'HandleVisibility','off'); hold on
%set(a,'EdgeColor',[0.3412 0.0157 0.3255],'FaceAlpha', .3, 'EdgeAlpha',.3)
plot([0 sm],bp_all_norm,'.-', 'Color', [0.3412 0.0157 0.3255], 'LineWidth',1); 
plot([0 sm],bp_all,'.-', 'Color', [0.7686 0.5333 0.9882], 'LineWidth',1); hold off
axis([0 max_sm 0 300])
xticks([0 sm])
ax = gca;
ax.FontSize = 8;
set(gca,'YMinorTick','on')
ylabel('Average Antarctic Response Time [years]', 'FontSize',10)
xlabel('Smoothing Window [Number of Data Points]', 'FontSize',10)
legend('Normalised','Not normalised')


max(bp_all)
min(bp_all)

%% CORES INDIVIDUALLY: Calculate Response Time for each Ice Core to calculate smoothing influence. 

t = -1200:1200;
min_sm = 5;
max_sm = 80;
sm = (min_sm:5:max_sm); % smoothing windows

idxbp = zeros(length(sm),1); % !!! CHANGE NAME preallocate vector for indeces of break points in stacks for each smoothing window
idxbp_norm = zeros(length(sm),1);
maxbp = zeros(length(sm),1);
bp_all = zeros(length(sm),1); %preallocate vector for break points (NON-normalised stack)
bp_all_norm = zeros(length(sm),1); % preallocate vector for break points (normalised stack)

COREstack = DFstack; % CHANGE FOR EACH CORE
COREstack_norm = normalize(COREstack,'range');

% determine break points for unsmoothed data
[~, idxCORE_unsmoothed] = max(COREstack(1201-300:1201));
idxCORE_unsmoothed = 1201-302+idxCORE_unsmoothed;
bp_unsmoothed = t(end-idxCORE_unsmoothed+1);
[~, idxCORE_unsmoothed] = max(COREstack_norm(1201-300:1201));
idxCORE_unsmoothed = 1201-302+idxCORE_unsmoothed;
bp_unsmoothed_norm = t(end-idxCORE_unsmoothed+1);

% determine breakpoints for each smoothing window
for j = 1:length(sm)
    COREstack_sm = smoothdata(COREstack,'movmean',sm(j));
    [~,idxbp(j)] = max(COREstack_sm(1201-300:1201));
    idxbp(j) = 1201-301+idxbp(j);
    bp_all(j) = t(end-idxbp(j)+1);
    COREstack_norm_sm = smoothdata(COREstack_norm,'movmean',sm(j));
    [maxbp(j),idxbp_norm(j)] = max(COREstack_norm_sm(1201-300:1201));
    idxbp_norm(j) = 1201-301+idxbp_norm(j);
    bp_all_norm(j) = t(end-idxbp_norm(j)+1);
end 

CORE_bp_all = vertcat(bp_unsmoothed,bp_all);
CORE_bp_all_norm = vertcat(bp_unsmoothed_norm,bp_all_norm);

% %% CORES INDIVIDUALLY: Plot Bipolar Phasing
% 
% idx_sm = 3; % idx of smoothing choosen (3 for window of 15)
% 
% figure('Name','Individual Response Time')
% tiledlayout(2,1,"TileSpacing",'none')
% GL = nexttile;
% a = fill([0, 0, t(end-idxbp(idx_sm-1)), t(end-idxbp(idx_sm-1))],[-3, 3, 3, -3], [0.9608 0.9020 0.7020], 'HandleVisibility','off'); hold on 
% a.EdgeColor = [0.9608 0.9020 0.7020];
% plot(-t,smoothdata(NGRIPstack_norm,'movmean',sm(idx_sm)),'Color',cNGRIP); hold on 
% ylabel(['Greenlandic \delta^{18}O [',char(8240),']']) 
% set(gca,'YMinorTick','on', 'XMinorTick','on')
% axis([-400 600 -0.1 1.1])
% set(GL, 'XColor','none')
% box off
% legend('NGRIP','Location','northeastoutside')
% 
% AA = nexttile; 
% a = fill([0, 0, t(end-idxbp(idx_sm-1)), t(end-idxbp(idx_sm-1))],[-0.7, 2, 2, -0.7], [0.9608 0.9020 0.7020], 'HandleVisibility','off'); hold on 
% a.EdgeColor = [0.9608 0.9020 0.7020];
% plot(-t,0.6*(smoothdata(COREstack_norm,'movmean',sm(idx_sm)))-0.3,'k'); hold on 
% plot(-t(idxbp(idx_sm-1)), 0.6*maxbp(idx_sm-1)-0.3, 'o','MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',4,'HandleVisibility','off')
% text(-t(idxbp(idx_sm-1)),0.6*maxbp(idx_sm-1)-0.22,sprintf('t = %.0f years',-t(idxbp(idx_sm-1))),'FontSize',8); hold off
% ax = gca;
% ax.FontSize = 8;
% ylabel(['Antarctic \delta^{18}O [',char(8240),']'], 'FontSize',10) 
% xlabel('Time [years]', 'FontSize',10)
% set(gca,'YMinorTick','on', 'XMinorTick','on')
% axis([-400 600 -0.425 0.5])
% box off
% 
% figure('Name','TEST SMOOTHING INFLUENCE')
% a = fill([0,0,max_sm,max_sm],[min(bp_all), max(bp_all), max(bp_all), min(bp_all)], [0.7686 0.5333 0.9882],'HandleVisibility','off'); hold on
% set(a,'EdgeColor',[0.7686 0.5333 0.9882], 'FaceAlpha', .2, 'EdgeAlpha',.2)
% %a = fill([0,0,max_sm, max_sm],[min(bp_all_norm),max(bp_all_norm),max(bp_all_norm),min(bp_all_norm)],[0.3412 0.0157 0.3255], 'HandleVisibility','off'); hold on
% %set(a,'EdgeColor',[0.3412 0.0157 0.3255],'FaceAlpha', .3, 'EdgeAlpha',.3)
% plot([0 sm],CORE_bp_all_norm,'.-', 'Color', [0.3412 0.0157 0.3255], 'LineWidth',1); 
% plot([0 sm],CORE_bp_all,'.-', 'Color', [0.7686 0.5333 0.9882], 'LineWidth',1); hold off
% axis([0 max_sm 0 300])
% xticks([0 sm])
% ax = gca;
% ax.FontSize = 8;
% %set(gca,'YMinorTick','on')
% ylabel('Average Antarctic Response Time [years]', 'FontSize',10)
% xlabel('Smoothing Window [Number of Data Points]', 'FontSize',10)
% legend('Normalised','Not normalised')
% disp(CORE_bp_all(idx_sm))
% %max(CORE_bp_all)-CORE_bp_all(idx_sm);
% %CORE_bp_all(idx_sm)-min(CORE_bp_all);

%% FIGURES
% The following sections are used to present the data, only to create figures.

%% FIG 4.1: ONSET DELAY NGRIP
% Determines delay of midpoint relative to actual onset of warming trend.

[minNGRIP,idxminNGRIP] = min(smoothdata(NGRIPstack(1200:1230),'movmean',1)); % finds minimum within 50 years before midpoint
%[minNGRIP, idxminNGRIP] = min(abs(NGRIPstack-minNGRIP));

figure('Name','Onset Delay NGRIP')
a = fill([0,0, -idxminNGRIP, -idxminNGRIP],[-100, 100, 100, -100], [0.93 0.93 0.93], 'HandleVisibility','off'); hold on 
a.EdgeColor = [0.93 0.93 0.93];
xline(t(end-idxminNGRIP+1), '--','Color',[0.7 0.7 0.7], 'HandleVisibility','off')
xline(0,'Color',[0.85 0.85 0.85], 'HandleVisibility','off')
text(idxminNGRIP-145,-41.5,sprintf('t = %.0f years',idxminNGRIP+1),'FontName','Times New Roman','FontSize',8) 
plot(-t, smoothdata(NGRIPstack,'movmean',23), 'Color',cNGRIP)
axis([-400 400 -43 -38])
yticks(-43:-38)
xticks(-1200:100:1200)
set(gca,'YMinorTick','on', 'XMinorTick','on')
ax = gca; 
ax.FontSize = 8; 
xlabel('Time [years]','FontSize',10)
ylabel(['\delta^{18}O [',char(8240),']'],'FontSize',10)


%% FIG D.2: AVERAGE DURATION OF DO EVENTS

[minGL,idxminGL_norm] = min(smoothdata(GLstack(1200:1250),'movmean',30)); % finds minumim 100 years prior midpoint in Greenland stack
[~, idxminGL_norm] = min(abs(GLstack-minGL));
[maxGL_norm, idxmaxGL_norm] = max(GLstack_norm); % finds maximum after midpoint

figure('Name','Onset Delay Greenland')
%a = fill([t(end-idxminGL_norm), t(end-idxminGL_norm), t(end-idxmaxGL_norm), t(end-idxmaxGL_norm)],[0, 1, 1, 0], [0.93 0.93 0.93], 'HandleVisibility','off'); hold on 
a = fill([-67, -67, 57, 57],[0, 1, 1, 0], [0.93 0.93 0.93], 'HandleVisibility','off'); hold on 
a.EdgeColor = [0.9490 0.9490 0.9490];
%xline(t(end-idxmaxGL_norm), '--','Color', [0.5 0.5 0.5], 'HandleVisibility','off')
%xline(t(end-idxminGL_norm),'--', 'Color', [0.5 0.5 0.5], 'HandleVisibility','off')
xline(0, 'Color', [0.85 0.85 0.85], 'HandleVisibility','off')
plot(-t,smoothdata(NGRIPstack_norm,'movmean',30),'Color',cNGRIP); hold on 
plot(-t,smoothdata(GRIPstack_norm,'movmean',30),'Color',cGRIP)
plot(-t,smoothdata(NEEMstack_norm,'movmean',30),'Color',cNEEM)
plot(-t,smoothdata(GISPstack_norm,'movmean',30),'Color',cGISP)
plot(-t,smoothdata(GLstack_norm,'movmean',30),'k','LineWidth',0.8)
%text(-t(idxminGL_norm)-120,0.2,sprintf('t = %.0f years',-t(idxminGL_norm)),'FontSize',6) 
%text(-t(idxmaxGL_norm)+1,0.7,sprintf('t = %.0f years',-t(idxmaxGL_norm)),'FontSize',6) 
text(-200,0.2,sprintf('t = %.0f years',-67),'FontSize',8) 
text(100,0.92,sprintf('t = %.0f years',57),'FontSize',8) 
xlim([-500 500])
xticks(-1200:100:1200)
ax = gca; 
ax.FontSize = 8;
set(gca,'YMinorTick','on', 'XMinorTick','on')
legend('NGRIP','GRIP','NEEM', 'GISP2','mean','Location','northeastoutside')
ylabel('Greenland Normalised \delta^{18}O', 'FontSize',10)
xlabel('Time [years]','FontSize',10)


%% FIG 4.4: STACKS OF INDIVIDUAL ICE CORES

xlims = [-400 400];
t = -1200:1200;
k1 = find(t == xlims(1));
k2 = find(t == xlims(2));


figure('units','centimeter','position',[5 1 15 20]); 
tiledlayout(6,1,'TileSpacing','none');

T_GL = nexttile; % GREENLAND STACKS
colororder([1 1 1; 0 0 0]); yyaxis right
plot(-t(k1:k2),NGRIPstack(k1:k2),'-','Color',cNGRIP); hold on
plot(-t(k1:k2),GRIPstack(k1:k2),'-','Color',cGRIP)
plot(-t(k1:k2),NEEMstack(k1:k2), '-', 'Color', cNEEM)
plot(-t(k1:k2),GISPstack(k1:k2), '-', 'Color', cGISP)
xline(0,'Color',[0.5 0.5 0.5],'LineWidth',2,'Alpha',.4, 'HandleVisibility','off'); hold off
set(gca,'XColor','none') % x-axis has no color
legend('NGRIP','GRIP','NEEM', 'GISP2','Location','northeastoutside')
ax = gca;
ax.FontSize = 8;
xlim([xlims(1) xlims(2)])

T_WDC = nexttile;
colororder([cWDC; 1 1 1]); yyaxis left;
plot(-t(k1:k2),WDCstack(k1:k2),'Color',cWDC); hold on
xline(0,'Color',[0.5 0.5 0.5],'LineWidth',2,'Alpha',.4, 'HandleVisibility','off'); hold off
set(gca, 'XColor','none')
xlim([xlims(1) xlims(2)])
ax = gca; 
ax.FontSize = 8; 
ylabel('WDC', 'FontSize',10)

T_TAL = nexttile;
colororder([1 1 1; cTAL]); yyaxis right;
plot(-t(k1:k2),TALstack(k1:k2),'Color',cTAL); hold on 
xline(0,'Color',[0.5 0.5 0.5],'LineWidth',2,'Alpha',.4, 'HandleVisibility','off'); hold off
set(gca, 'XColor', 'none', 'YColor',cTAL)
xlim([xlims(1) xlims(2)])
ax = gca;
ax.FontSize = 8;
ylabel('TAL', 'FontSize',10)

T_DF = nexttile; 
colororder([cDF; 1 1 1]); yyaxis left; 
plot(-t(k1:k2),DFstack(k1:k2),'Color',cDF); hold on 
xline(0,'Color',[0.5 0.5 0.5],'LineWidth',2,'Alpha',.4, 'HandleVisibility','off'); hold off
set(gca,'XColor','none')
xlim([xlims(1) xlims(2)])
ax = gca;
ax.FontSize = 8;
ylabel('Dome F','FontSize',10) 

T_EDC = nexttile; 
colororder([1 1 1;cEDC]); yyaxis right;
plot(-t(k1:k2),EDCstack(k1:k2), 'Color',cEDC); hold on 
xline(0,'Color',[0.5 0.5 0.5],'LineWidth',2,'Alpha',.4, 'HandleVisibility','off'); hold off
set(gca, 'XColor','none', 'YColor',cEDC)
xlim([xlims(1) xlims(2)])
ax = gca;
ax.FontSize = 8;
ylabel('EDC', 'FontSize',10)

T_EDML = nexttile; 
colororder([cEDML; 1 1 1]); yyaxis left;
set(gca, 'ytick',[], 'yticklabel',[]);
plot(-t(k1:k2), EDMLstack(k1:k2),'Color',cEDML); hold on
xline(0,'Color',[0.5 0.5 0.5],'LineWidth',2,'Alpha',.4, 'HandleVisibility','off'); hold off
box off
xlim([xlims(1) xlims(2)])
xticks(t(1):100:t(end))
ax = gca;
ax.FontSize = 8;
ylabel('EDML', 'FontSize',10) 
xlabel('Time [Years]', 'FontSize',10)
set(gca,'XMinorTick','on')


%% FIG C.2 and FIG 4.8: CHANGE IN d18O OVER TRANSITION

a = 0.85; %sets alpha

figure('Name','Delta d18O,all')
plot(1:numel(events), NGRIPdeld18O, 'o-','MarkerEdgeColor',cNGRIP, 'MarkerFaceColor', cNGRIP,'Color', cNGRIP, 'MarkerSize',3); hold on
plot(1:numel(events), GRIPdeld18O, 'o-','MarkerEdgeColor',cGRIP, 'MarkerFaceColor', cGRIP, 'Color', cGRIP, 'MarkerSize',3);
plot(1:numel(events), NEEMdeld18O, 'o-','MarkerEdgeColor',cNEEM, 'MarkerFaceColor', cNEEM,'Color', cNEEM, 'MarkerSize',3);
plot(1:numel(events), GISPdeld18O, 'o-','MarkerEdgeColor',cGISP, 'MarkerFaceColor', cGISP,'Color', cGISP, 'MarkerSize',3);

plot(1:numel(events), WDCdeld18O, 'o-','MarkerEdgeColor',cWDC, 'MarkerFaceColor', cWDC, 'Color', cWDC, 'MarkerSize',3); hold on
plot(1:numel(events), TALdeld18O, 'o-','MarkerEdgeColor',cTAL, 'MarkerFaceColor', cTAL,'Color', cTAL, 'MarkerSize',3);
plot(1:numel(events), DFdeld18O, 'o-','MarkerEdgeColor',cDF, 'MarkerFaceColor', cDF, 'Color', cDF, 'MarkerSize',3);
plot(1:numel(events), EDCdeld18O, 'o-','MarkerEdgeColor',cEDC, 'MarkerFaceColor', cEDC,'Color', cEDC, 'MarkerSize',3);
plot(1:numel(events), EDMLdeld18O, 'o-','MarkerEdgeColor',cEDML, 'MarkerFaceColor', cEDML, 'Color', cEDML, 'MarkerSize',3);
yline(0,'Color',[0.5 0.5 0.5],'LineWidth',1,'Alpha',.4, 'HandleVisibility','off')
xlim([0.5 33.5])
xticks(1:numel(events))
labels = {'YD-PB', 'OD-BA', '2.2', '3', '4','5.1','5.2','6','7c','8c',...
    '9','10','11','12','13c','14e','15.1','15.2','16.1','16.2', '17.1c','17.2',...
    '18', '19.1', '19.2', '20c', '21.1e', '21.2','22g','23.1','23.2','24.1c','24.2'};
set(gca, 'XTickLabel',labels, 'FontSize',8) %,'FontSize',9)
set(gca,'YMinorTick','on')
ylabel(['\Delta \delta^{18}O [',char(8240),']'],'FontSize',10)
legend('NGRIP', 'GRIP', 'GISP2','NEEM', 'WDC', 'TAL', 'DOME F', 'EDC', 'EDML','Location','northeastoutside')



%two plots for Greenland and Antarctica seperatly 
Y_GL = [mean(NGRIPdeld18O(allDO),'omitnan') mean(GRIPdeld18O(allDO),'omitnan') mean(GISPdeld18O(allDO),'omitnan') mean(NEEMdeld18O(allDO),'omitnan') ...
    ; mean(NGRIPdeld18O(secondhalf),'omitnan') mean(GRIPdeld18O(secondhalf),'omitnan') mean(GISPdeld18O(secondhalf),'omitnan') mean(NEEMdeld18O(secondhalf),'omitnan') ...
    ; mean(NGRIPdeld18O(firsthalf), 'omitnan') mean(GRIPdeld18O(firsthalf), 'omitnan') mean(GISPdeld18O(firsthalf), 'omitnan') mean(NEEMdeld18O(firsthalf), 'omitnan')...
    ; mean(NGRIPdeld18O(minorDO),'omitnan') mean(GRIPdeld18O(minorDO),'omitnan') mean(GISPdeld18O(minorDO),'omitnan') mean(NEEMdeld18O(minorDO),'omitnan') ...
    ; mean(NGRIPdeld18O(majorDO),'omitnan') mean(GRIPdeld18O(majorDO),'omitnan') mean(GISPdeld18O(majorDO),'omitnan') mean(NEEMdeld18O(majorDO),'omitnan')];

Y_AA = [mean(WDCdeld18O(allDO),'omitnan') mean(TALdeld18O(allDO),'omitnan') mean(DFdeld18O(allDO),'omitnan') mean(EDCdeld18O(allDO),'omitnan') mean(EDMLdeld18O(allDO),'omitnan') ...
    ; mean(WDCdeld18O(secondhalf),'omitnan') mean(TALdeld18O(secondhalf),'omitnan') mean(DFdeld18O(secondhalf),'omitnan') mean(EDCdeld18O(secondhalf),'omitnan') mean(EDMLdeld18O(secondhalf),'omitnan') ...
    ; mean(WDCdeld18O(firsthalf), 'omitnan') mean(TALdeld18O(firsthalf), 'omitnan') mean(DFdeld18O(firsthalf), 'omitnan'), mean(EDCdeld18O(firsthalf), 'omitnan') mean(EDMLdeld18O(firsthalf), 'omitnan') ...
    ; mean(WDCdeld18O(minorDO),'omitnan') mean(TALdeld18O(minorDO),'omitnan') mean(DFdeld18O(minorDO),'omitnan') mean(EDCdeld18O(minorDO),'omitnan') mean(EDMLdeld18O(minorDO),'omitnan') ...
    ; mean(WDCdeld18O(majorDO),'omitnan') mean(TALdeld18O(majorDO),'omitnan') mean(DFdeld18O(majorDO),'omitnan') mean(EDCdeld18O(majorDO),'omitnan') mean(EDMLdeld18O(majorDO),'omitnan')];

figure('units','centimeter','position',[0 2 40 15],"Name",'Change in d18O over transitions')
subplot(1,2,1) % GREENLAND 
b=bar(Y_GL);
b(1).FaceColor = cNGRIP; b(2).FaceColor = cGRIP; b(3).FaceColor = cGISP; b(4).FaceColor = cNEEM;
alpha(a)
ylabel(['\Delta \delta^{18}O [',char(8240),']'],'FontSize',10)
ylim([0 5.5])
xticklabels(gca,{'10-115ka b2k','10-65ka b2k','65-115ka b2k','Minor DOs', 'Major DOs'})
ax = gca; 
ax.FontSize = 8;
legend('NGRIP', 'GRIP', 'GISP2', 'NEEM', 'Location','northwestoutside')
set(gca,'YMinorTick','on')

subplot(1,2,2) % ANTARCTICA
b = bar(Y_AA);
b(1).FaceColor = cWDC; b(2).FaceColor = cTAL; b(3).FaceColor = cDF; b(4).FaceColor = cEDC; b(5).FaceColor = cEDML;
alpha(a)
ylim([-0.2 0.8])
yticks(-0.2:0.1:0.8)
xticklabels(gca,{'10-115ka b2k','10-65ka b2k','65-115ka b2k','Minor DOs', 'Major DOs'})
ylabel(['\Delta \delta^{18}O [',char(8240),']'],'FontSize',10)
legend('WDC', 'TAL', 'DOME F', 'EDC', 'EDML', 'Location', 'northeastoutside')
ax = gca; 
ax.FontSize = 8;
xtickangle(30)
set(gca,'YMinorTick','on')
print(gcf,'C:\Users\annak\OneDrive - University of Copenhagen\Desktop\Speciale\Figures\delta_d18O_bar.png','-dpng','-r600')


% AVERAGE FOR ALL ANTARCTIC ICE CORES
AAdeld18O = [WDCdeld18O TALdeld18O DFdeld18O EDCdeld18O EDMLdeld18O];
AAdeld18O_mean = mean(abs(AAdeld18O),2,'omitnan'); % absolute average change across the transition for each event
% !!!!!!! NEEDED ??????

%% FIG. E.1: CLIMATE SIGNAL ACROSS TRANSITIONS


volc = readmatrix('bipolar_volc.xlsx'); 
volc_age = volc(:,2);

% Plots all events including onsets and bipolar matchpoints
for i = 1:length(events)
    t = (-1200:1200)+age_trans(i);
    figure('units','centimeter','position',[10 1.5 15 15])
    plot(t,(NGRIPd18Otimed(i,:)-min(NGRIPd18Otimed(i,:)))/(max(NGRIPd18Otimed(i,:))-min(NGRIPd18Otimed(i,:)))+4.5,'Color',[cNGRIP 0.3],'HandleVisibility','off'); hold on 
    plot(t,smoothdata((NGRIPd18Otimed(i,:)-min(NGRIPd18Otimed(i,:)))/(max(NGRIPd18Otimed(i,:))-min(NGRIPd18Otimed(i,:))),'movmean',15)+4.5,'Color', cNGRIP);
    plot(t,(GRIPd18Otimed(i,:)-min(GRIPd18Otimed(i,:)))/(max(GRIPd18Otimed(i,:))-min(GRIPd18Otimed(i,:)))+4,'Color',[cGRIP 0.3],'HandleVisibility','off')
    plot(t,smoothdata((GRIPd18Otimed(i,:)-min(GRIPd18Otimed(i,:)))/(max(GRIPd18Otimed(i,:))-min(GRIPd18Otimed(i,:))),'movmean',9)+4,'Color',cGRIP)
    plot(t,(NEEMd18Otimed(i,:)-min(NEEMd18Otimed(i,:)))/(max(NEEMd18Otimed(i,:))-min(NEEMd18Otimed(i,:)))+3.5,'Color',[cNEEM 0.3],'HandleVisibility','off')
    plot(t,smoothdata((NEEMd18Otimed(i,:)-min(NEEMd18Otimed(i,:)))/(max(NEEMd18Otimed(i,:))-min(NEEMd18Otimed(i,:))),'movmean',9)+3.5,'Color',cNEEM)
    plot(t,(GISPd18Otimed(i,:)-min(GISPd18Otimed(i,:)))/(max(GISPd18Otimed(i,:))-min(GISPd18Otimed(i,:)))+3,'Color',[cGISP 0.3],'HandleVisibility','off')
    plot(t,smoothdata((GISPd18Otimed(i,:)-min(GISPd18Otimed(i,:)))/(max(GISPd18Otimed(i,:))-min(GISPd18Otimed(i,:))),'movmean',7)+3,'Color',cGISP)
    plot(t,(WDCd18Otimed(i,:)-min(WDCd18Otimed(i,:)))/(max(WDCd18Otimed(i,:))-min(WDCd18Otimed(i,:)))+2,'Color',[cWDC 0.3], 'HandleVisibility','off')
    plot(t,smoothdata((WDCd18Otimed(i,:)-min(WDCd18Otimed(i,:)))/(max(WDCd18Otimed(i,:))-min(WDCd18Otimed(i,:))), 'movmean', 10)+2,'Color',cWDC)
    plot(t,(TALd18Otimed(i,:)-min(TALd18Otimed(i,:)))/(max(TALd18Otimed(i,:))-min(TALd18Otimed(i,:)))+1.5, 'Color', [cTAL 0.3], 'HandleVisibility','off')
    plot(t,smoothdata((TALd18Otimed(i,:)-min(TALd18Otimed(i,:)))/(max(TALd18Otimed(i,:))-min(TALd18Otimed(i,:))), 'movmean', 10)+1.5, 'Color',  cTAL)
    plot(t,(DFd18Otimed(i,:)-min(DFd18Otimed(i,:)))/(max(DFd18Otimed(i,:))-min(DFd18Otimed(i,:)))+1,'Color', cDF)
    plot(t,(EDCd18Otimed(i,:)-min(EDCd18Otimed(i,:)))/(max(EDCd18Otimed(i,:))-min(EDCd18Otimed(i,:)))+0.5,'Color', [cEDC 0.3], 'HandleVisibility','off')
    plot(t,smoothdata((EDCd18Otimed(i,:)-min(EDCd18Otimed(i,:)))/(max(EDCd18Otimed(i,:))-min(EDCd18Otimed(i,:))),'movmean', 10)+0.5,'Color', cEDC)
    plot(t,(EDMLd18Otimed(i,:)-min(EDMLd18Otimed(i,:)))/(max(EDMLd18Otimed(i,:))-min(EDMLd18Otimed(i,:))), 'Color',cEDML)
    xline(age_trans,'Color',[0.5 0.5 0.5],'LineWidth',1.2,'Alpha',.4); 
    xline(volc_age(:),'Color',[0.8 0.7 0.1],'LineWidth',1.2,'Alpha',.3); hold off
    set(gca,'ytick',[])
    ax=gca; 
    ax.FontSize = 8;
    legend('NGRIP','GRIP','NEEM','GISP2','WDC','TAL','DF','EDC','EDML','Location','northeastoutside') %'WDC','TAL','Dome F',
    xlabel('GICC05modelxt Age [Years b2k]', 'FontSize',10)
    xlim([round(age_trans(i)-1000,-2) round(age_trans(i)+1000,-2)]) 
    ax.XAxis.Exponent = 0;
    ax.XTick = round(age_trans(i),-2)-1000:400:round(age_trans(i)+1000,-2);
    set(ax,'XMinorTick','on')
    NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
    %text(NW(1)+100, NW(2), events(i), 'VerticalAlignment','top','HorizontalAlignment',...
    % 'left', 'FontSize',8,'FontName','Times New Roman')
    text(age_trans(i)-20,NW(2), events(i), 'HorizontalAlignment','right','FontSize',8)
    %text(age_trans(31)-20,NW(2), events(31), 'HorizontalAlignment','right','FontSize',8)
    %print(gcf,'C:\Users\annak\OneDrive - University of Copenhagen\Desktop\Speciale\Figures\climate_GI242.png','-dpng','-r600')
end

%% FIG 4.7 & FIG G.4 : COMPARE CLIMATE BEHAVIOUR ACROSS TRANSITION FOR DIFFERENT GROUPS OF EVENTS FOR SAME ICE CORE
% Plots stacked d18O record of events grouped in first vs second half of last ice age, and minor vs major events for each individual ice core.
% Uses function myclasscomparison.m

sm = 30;
sm_min = 1;
sm_max = 40;
t = -1200:1200;

if length(ofinterest) == numel(age_trans) % needs all transitions, thus, it runs only when all events are of interest
   [NGRIPstack_second, NGRIPstack_first, NGRIPstack_minor, NGRIPstack_major] = myclasscomparison(age_trans, NGRIPd18Otimed,'NGRIP', -45, -36, sm, sm_max, sm);
   [GRIPstack_second, GRIPstack_first, GRIPstack_minor, GRIPstack_major] = myclasscomparison(age_trans, GRIPd18Otimed,'GRIP', -42.5, -36, sm, sm_max, sm);
   [NEEMstack_second, NEEMstack_first, NEEMstack_minor, NEEMstack_major] = myclasscomparison(age_trans, NEEMd18Otimed,'NEEM',-45, -36, sm, sm_max, sm);
   [GISPstack_second, GISPstack_first, GISPstack_minor, GISPstack_major] = myclasscomparison(age_trans, GISPd18Otimed,'GISP2', -42, -36, sm, sm_max, sm_min);
   [WDCstack_second, WDCstack_first, WDCstack_minor, WDCstack_major] = myclasscomparison(age_trans, WDCd18Otimed,'WDC',-39.8, -38, sm, sm_max, sm_min);
   [TALstack_second, TALstack_first, TALstack_minor, TALstack_major] = myclasscomparison(age_trans, TALd18Otimed,'Talos Dome',-41, -38.5, 50, 50, sm_min);
   [DFstack_second, DFstack_first, DFstack_minor, DFstack_major] = myclasscomparison(age_trans, DFd18Otimed,'Dome F', -60.5, -56.5, sm_min, sm_min, sm_min);
   [EDCstack_second, EDCstack_first, EDCstack_minor, EDCstack_major] = myclasscomparison(age_trans, EDCd18Otimed,'EDC', -55.5, -52.5, sm, sm_max, sm);
   [EDMLstack_second, EDMLstack_first, EDMLstack_minor, EDMLstack_major] = myclasscomparison(age_trans, EDMLd18Otimed,'EDML', -48.8, -47, sm_min, sm_min, sm_min);
end 

% Average for Greenland and Antarctica 
GL_seoncd = [NGRIPstack_second; GRIPstack_second; NEEMstack_second; GISPstack_second];
GLstack_second = mean(GL_seoncd,'omitnan');
GL_first = [NGRIPstack_first; GRIPstack_first; NEEMstack_first; GISPstack_first];
GLstack_first = mean(GL_first,'omitnan');
GL_min = [NGRIPstack_minor; GRIPstack_minor; NEEMstack_minor; GISPstack_minor];
GLstack_min = mean(GL_min, 'omitnan');
GL_maj = [NGRIPstack_major; GRIPstack_major; NEEMstack_major; GISPstack_major];
GLstack_maj = mean(GL_maj, 'omitnan');

AA_second = [EDCstack_second; EDMLstack_second]; %[WDCstack_second; TALstack_second; DFstack_second;
AAstack_second = mean(AA_second,'omitnan');
AA_first = [EDCstack_first; EDMLstack_first]; %WDCstack_first; TALstack_first; DFstack_first
AAstack_first = mean(AA_first,'omitnan');
AA_min = [WDCstack_minor; TALstack_minor; DFstack_minor; EDCstack_minor; EDMLstack_minor];
AAstack_min = mean(AA_min, 'omitnan');
AA_maj = [WDCstack_major; TALstack_major; DFstack_major; EDCstack_major; EDMLstack_major];
AAstack_maj = mean(AA_maj,'omitnan');

% First vs secodn half of Last Glacial
figure('units','centimeter','position',[0 2 40 15], 'Name', 'STACKED First vs Second')
subplot(1,2,1) % Greenland
plot(-t,GLstack_second,'Color',[csecond 0.3],'HandleVisibility','off'); hold on
plot(-t,smoothdata(GLstack_second,'movmean',sm),'Color',csecond)
plot(-t, GLstack_first,'Color',[cfirst 0.3],'HandleVisibility','off')
plot(-t,smoothdata(GLstack_first,'movmean',sm_min),'Color',cfirst)
xline(0,'Color',[0.5 0.5 0.5],'LineWidth',2,'Alpha',.4, 'HandleVisibility','off'); hold off 
ax = gca;
ax.XTick = t(1):200:t(end);
set(gca, 'XMinorTick','on','YMinorTick','on')
ax.XAxis.MinorTickValues = t(1):50:t(end);
ax.FontSize = 8;
ylabel(['\delta^{18}O [',char(8240),']'],'FontSize',10) 
xlabel('Time [Years]','FontSize',10)
axis([-800 800 -42 -37]) 
yticks(-43:-37)
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
text(NW(1), NW(2), 'Greenland', 'VerticalAlignment','top', 'HorizontalAlignment',...
    'left', 'FontSize',8, 'color', 'k')

subplot(1,2,2) % Antarctica
plot(-t,AAstack_second,'Color',[csecond 0.3],'HandleVisibility','off'); hold on
plot(-t,smoothdata(AAstack_second,'movmean',sm),'Color',csecond)
plot(-t, AAstack_first,'-','Color',[cfirst 0.3],'HandleVisibility','off')
plot(-t,smoothdata(AAstack_first,'movmean',20),'-','Color',cfirst)
xline(0,'Color',[0.5 0.5 0.5],'LineWidth',2,'Alpha',.4, 'HandleVisibility','off'); hold off %[0.9290 0.6940 0.1250]
legend('10-65 ka b2k','65-115 ka b2k')
ax = gca;
ax.XTick = t(1):200:t(end);
set(gca, 'XMinorTick','on','YMinorTick','on')
ax.XAxis.MinorTickValues = t(1):50:t(end);
ax.FontSize = 8;
ylabel(['\delta^{18}O [',char(8240),']'],'FontSize',10) 
xlabel('Time [Years]','FontSize',10)
axis([-800 800 -52 -50]) 
yticks(-52:0.4:-50)
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
text(NW(1), NW(2), 'Antarctica', 'VerticalAlignment','top', 'HorizontalAlignment',...
    'left', 'FontSize',8, 'color', 'k')



% Minor vs Major DOs
figure('units','centimeter','position',[0 2 40 15], 'Name', 'STACKED Minor vs Major')
subplot(1,2,1) % Greenland
plot(-t,GLstack_maj,'Color',[cmajor 0.3],'HandleVisibility','off'); hold on
plot(-t,smoothdata(GLstack_maj,'movmean',sm_max),'Color',cmajor)
plot(-t, GLstack_min,'Color',[cminor 0.3],'HandleVisibility','off')
plot(-t, smoothdata(GLstack_min,'movmean',sm),'Color',cminor)
xline(0,'Color',[0.5 0.5 0.5],'LineWidth',2,'Alpha',.4, 'HandleVisibility','off'); hold off %[0.9290 0.6940 0.1250]
ax = gca;
ax.XTick = t(1):200:t(end);
ax.XAxis.MinorTickValues = t(1):50:t(end);
ax.FontSize = 8;
set(gca, 'YMinorTick','on','XMinorTick','off')
ylabel(['\delta^{18}O [',char(8240),']'],'FontSize',10) 
xlabel('Time [Years]','FontSize',10)
axis([-800 800 -42 -37])
yticks(-42:-37)
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
text(NW(1), NW(2), 'Greenland', 'VerticalAlignment','top', 'HorizontalAlignment',...
    'left', 'FontSize',8, 'color', 'k')

subplot(1,2,2) % Antarctica 
plot(-t,AAstack_maj,'Color',[cmajor 0.3],'HandleVisibility','off'); hold on
plot(-t,smoothdata(AAstack_maj,'movmean',sm_max),'Color',cmajor)
plot(-t, AAstack_min,'Color',[cminor 0.3],'HandleVisibility','off')
plot(-t, smoothdata(AAstack_min,'movmean',sm),'Color',cminor)
xline(0,'Color',[0.5 0.5 0.5],'LineWidth',2,'Alpha',.4, 'HandleVisibility','off'); hold off %[0.9290 0.6940 0.1250]
legend('Major','Minor')
ax = gca;
set(gca, 'YMinorTick','on','XMinorTick','off')
ax.XTick = t(1):200:t(end);
ax.XAxis.MinorTickValues = t(1):50:t(end);
ax.FontSize = 8;
ylabel(['\delta^{18}O [',char(8240),']'],'FontSize',10) 
xlabel('Time [Years]','FontSize',10)
axis([-800 800 -48.5 -47])
yticks(-48.5:0.3:-47)
NW = [min(xlim) max(ylim)]+[diff(xlim) -diff(ylim)]*0.05;
text(NW(1), NW(2), 'Antarctica', 'VerticalAlignment','top', 'HorizontalAlignment',...
    'left', 'FontSize',8, 'color', 'k')


%% FIG D.3 UNCERTAINTY ESTIMATION
% 1. MATCH UNCERTAINTY

figure('Name','Match Uncertainty')
h = histogram(d_bp_match);
h.BinWidth = 2.5;
ax= gca; 
ax.FontSize = 8; 
xticks(0:5:60)
yticks(1:12)
xlabel('Match Uncertainty [years]', 'FontSize',10)

d_bpmatch_lg = median(d_bp_match);
d_bpmatch_first = median(d_bp_match(firsthalf));
d_bpmatch_sec = median(d_bp_match(secondhalf));
d_bpmatch_minor = median(d_bp_match(minorDO));
d_bpmatch_major = median(d_bp_match(majorDO));
% mean_bpmatch = mean(d_bp_match);
% std_match = std(sigma_bp_match);

%% 2 RESOLUTION
% mean temporal resolution based on equidistant samples, achieved from timescale from interpolation. Estimated in resolution.m

% mean temporal resolution for WDC-TAL-DF-EDC-EDML
s_res_first = [NaN NaN 56 11 46]; % first half
s_res_sec = [17 12 36 9 24]; % second half
s_res_lg = [17 12 56 11 46]; % whole Last Glacial Period


%% 3 SMOOTHING INFLUENCE
% maximum deviation of response time determined for different smoothing windows.
% CHANGE FOR EACH CORE. DEFINE LIMITS AND AREA FOR EACH CORE.
% Produces FIG D.4 for WDC 

AAstack2 = (WDCstack + TALstack + DFstack + EDCstack + EDMLstack)/5;

t = -1200:1200;
sm = 80;
COREstack = WDCstack;
cCORE = cWDC;
res_CORE = (COREstack - smoothdata(COREstack, 'movmean',sm));
std(res_CORE);
sigma = std(res_CORE);
COREstack_sm_up = (smoothdata(COREstack,'movmean',sm)+sigma);
[max_d, idx_max_d] = max(smoothdata(COREstack, 'movmean',sm));
[max_d18O,idx_max_up] = mink(abs(max_d - (smoothdata(COREstack,'movmean',sm)+sigma)),2);

disp(-t(idx_max_up))
d_max_CORE = (abs(t(idx_max_up(1))-t(idx_max_up(2))))/2;

figure('Name','CORE')
subplot(3,8,[1:6 9:14 17:22])
a = fill([-t fliplr(-t)], [(smoothdata(COREstack,'movmean',sm))-std(res_CORE) fliplr(smoothdata(COREstack,'movmean',sm)+std(res_CORE))], cCORE); hold on
a.EdgeColor = 'none';
alpha(0.15)
plot(-t, COREstack,'Color',[cCORE 0.5]); hold on
plot(-t, smoothdata(COREstack,'movmean',sm), 'Color',cCORE,'LineWidth',2)
plot(-t,(smoothdata(COREstack,'movmean',sm)+sigma), 'Color', [cCORE 0.2])
plot(-t,(smoothdata(COREstack,'movmean',sm)-sigma), 'Color', [cCORE 0.2])
%plot(-t(idx_max_d), max_d, 'o','MarkerFaceColor','r','MarkerEdgeColor','w','MarkerSize',4,'HandleVisibility','off')
plot(-t(idx_max_up(1)), COREstack_sm_up(idx_max_up(1)),'o', 'MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',4,'HandleVisibility','off')
plot(-t(idx_max_up(2)), COREstack_sm_up(idx_max_up(2)),'o', 'MarkerFaceColor','k','MarkerEdgeColor','w','MarkerSize',4,'HandleVisibility','off')
yline(max_d)
xlim([-1200 1200])
ylim([-39.65 -38.55])
xticks(-1200:400:1200)
yticks(-39.6:0.1:-38.5)
xline(0,'Color',[0.5 0.5 0.5],'LineWidth',2,'Alpha',.4, 'HandleVisibility','off'); hold off
ax = gca; 
ax.FontSize = 8;
set(gca,'YMinorTick','on', 'XMinorTick','on')
xlabel('Time [years]', 'FontSize',10)
ylabel(['\delta^{18}O [',char(8240),']'],'FontSize',10) 

area = [-100 -38.74 180 -38.64];
inlarge = subplot(3,8,[7 8]);
panpos = inlarge.Position;
delete(inlarge);
inlarge = zoomin(ax,area,panpos);



%% Uncertainty for each Core 
% combining uncertainties with error propagation

s_Cores_first = sqrt((ones(1,numel(s_res_first))*std_err_first).^2 + s_res_first.^2 + s_sm_first.^2);
s_Cores_sec = sqrt((ones(1,numel(s_res_sec))*std_err_sec).^2 + s_res_sec.^2 + s_sm_sec.^2);
s_Cores_lg = sqrt((ones(1,numel(s_res_lg))*std_err_lg).^2 + s_res_lg.^2 + s_sm_lg.^2);

s_Cores_first(isnan(s_Cores_first))=0; % replace NaN by zero
s_allC_first = sqrt(s_Cores_first(1)^2 + s_Cores_first(2)^2 + s_Cores_first(3)^2 + s_Cores_first(4)^2 + s_Cores_first(5)^2); 
s_Cores_sec(isnan(s_Cores_sec))=0; % replace NaN by zero
s_allC_sec = sqrt(s_Cores_sec(1)^2 + s_Cores_sec(2)^2 + s_Cores_sec(3)^2 + s_Cores_sec(4)^2 + s_Cores_sec(5)^2);

s_Cores_lg(isnan(s_Cores_lg))=0; % replace NaN by zero
s_allC_lg = sqrt(s_Cores_lg(1)^2 + s_Cores_lg(2)^2 + s_Cores_lg(3)^2 + s_Cores_lg(4)^2 + s_Cores_lg(5)^2);


