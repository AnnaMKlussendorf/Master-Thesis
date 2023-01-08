function [age, depth, COREstack, COREd18Otimed, COREdeltad18O, avg_stadial, avg_interstadial, COREstackdeltad18O] = mystacking(Chrono_age, Chrono_depth, depth, d18O, events, age_trans, CORE, ofinterest,cCORE)

% not interested in Holocene, start at around 10k a b2k
i10k = find(Chrono_age >= 10000,1);
Chrono_age = Chrono_age(i10k:end);
Chrono_depth = Chrono_depth(i10k:end);
i10k = find(depth >= Chrono_depth(1),1);
depth = depth(i10k:end);
d18O = d18O(i10k:end);

% omit possible NaN-values in depth and age 
idx = find(~isnan(Chrono_depth) & ~isnan(Chrono_age));
age = interp1(Chrono_depth(idx), Chrono_age(idx), depth); % finds age corresponding to the depth and d18O high-resolution data

% stacking 
t = -1200:1200;

% !!!!! OLD ORIGINAL METHOD: INTERPOLATE BEFORE SMOOTHING
COREd18Otimed = NaN(numel(events),length(t)); % preallocate matrix with number of rows = number of events (34x101)
for i = ofinterest
    COREd18Otimed(i,:) = interp1(age(~isnan(age))-age_trans(i),d18O(~isnan(age)),t); %+age_trans(i); % change to MINUS again!!! % age has NaN values, omit when interpolate the climate signal onto t
end 

% % !!!!!!! ACHTUNG TEST: NEEDS TO BE CHANGED AGAIN. TEST WHETHER IT MAKES SENSE TO
% % SMOOTH THE DATA BEFORE INTERPOLATING 
% 
% COREd18Otimed = NaN(numel(events),length(t)); % preallocate matrix with number of rows = number of events (34x101)
% for i = ofinterest
%     COREd18Otimed(i,:) = interp1(age(~isnan(age))-age_trans(i),smoothdata(d18O(~isnan(age),'movmean',10)),t); %+age_trans(i); % change to MINUS again!!! % age has NaN values, omit when interpolate the climate signal onto t
% end 

COREd18Otimed = smoothdata(COREd18Otimed,2,'movmean',10);

COREstack = mean(COREd18Otimed,'omitnan');

figure('visible','off'); % NOT normalised
plot(-t,COREstack,'Color',[cCORE 0.2]); hold on
plot(-t, smoothdata(COREstack,'movmean',20),'Color', cCORE)
xline(0,'Color',[0.5 0.5 0.5],'LineWidth',1.2,'Alpha',.4,'HandleVisibility','off'); hold off
grid on 
xlabel('Years')
ylabel(['\delta^{18}O [',char(8240),']']) 
xlim([t(1) t(end)])
title(CORE)
ax = gca;
ax.XTick = t(1):300:t(end);
ax.XMinorTick='on';
grid minor
ax.MinorGridLineStyle = '-';
ax.MinorGridAlpha = 0.1;
set(gca,'FontName','Times New Roman')

% DETERMINE ABSOLUTE CHANGE IN d18O OVER EACH TRANSITION

COREdeltad18O = zeros(numel(events),1);
for i = 1:numel(events)
    if i == 19
    avg_stadial = mean(COREd18Otimed(i,1201+30:1201+130)); % only 100 years for GI-16.1
    else
    avg_stadial = mean(COREd18Otimed(i,1201+30:1201+280));
    end
    avg_interstadial = mean(COREd18Otimed(i, 1201-75:1201-25)); % 50 to 100
    COREdeltad18O(i) = (avg_interstadial - avg_stadial);    
end 

% DETERMINE ABSOLUTE CHANGE IN d18O OVER STACKED TRANSITIONS 

avg_stadial = mean(COREstack(1,1201+30:1201+130));
avg_interstadial = mean(COREstack(1, 1201-75:1201-25)); % 50 to 100
COREstackdeltad18O = (avg_interstadial - avg_stadial);  



k = 14; % event 12
for i = k
    avg_stadial = mean(COREd18Otimed(i, 1201+30:1201+129));
    avg_interstadial = mean(COREd18Otimed(i, 1201-74:1201-25));
    %COREdeltad18O_19 = (avg_interstadial - avg_stadial);
end 

avg_stadial = ones(100,1)*avg_stadial;
avg_interstadial = ones(50,1)*avg_interstadial;

figure('Visible','off')
%plot(-t,COREd18Otimed(k,:),'Color',[cCORE 0.2]); hold on
plot(-t, smoothdata(COREd18Otimed(k,:),'movmean',40),'Color', cCORE); hold on
plot(-t(1201+30:1201+129), avg_stadial, 'Color', [0.9882 0.8235 0.2196], 'LineWidth',1.2)
plot(-t(1201-74:1201-25),avg_interstadial, 'Color', [0.9882 0.8235 0.2196], 'LineWidth',1.2)
xline(0,'Color',[0.5 0.5 0.5],'LineWidth',1.2,'Alpha',.4,'HandleVisibility','off'); hold off
axis([-600 600 -43 -35.5])
ax = gca; 
ax.FontSize = 8;
xlabel('Time [Years]','FontSize',10)
ylabel(['\delta^{18}O [',char(8240),']'],'FontSize',10) 
set(gca, 'FontName', 'Times New Roman')
set(gca, 'XMinorTick','on', 'YMinorTick', 'on')
%print(gcf,'C:\Users\annak\OneDrive - University of Copenhagen\Desktop\Speciale\Figures\det_Dd18O_GRIP.png','-dpng','-r600')



