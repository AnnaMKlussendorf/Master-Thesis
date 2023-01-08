function [depth_trans,diff_d18O] = myfindtransitions_single(depth_lim, depth, d18O, sm, n, DO_depth, events) % input variables

L = 5;
depth_trans = zeros(length(events),1); 
cNGRIP = [0 0.4470 0.7410];

for i = 1:numel(events)
    k1 = find(depth > depth_lim(1,i),1);
    k2 = find(depth > depth_lim(2,i),1);
    depthI = depth(k1:k2);
    d18Osm = smoothdata(d18O,'movmean',sm(i)); % smoothed data
    d18OI_sm = d18Osm(k1:k2);

    % find differences over larger interval: returns indices and depth of transitions
    diff_d18O = zeros(length(d18OI_sm),1);
    for j = L+1:length(d18OI_sm)-L
    diff_d18O(j) = d18OI_sm(j+L)-d18OI_sm(j-L);
    end 

    %n = 1; 
    [~,idx] = mink(diff_d18O,1);
    depth_trans(i) = depthI(idx);
    %depth_trans_all(i) = depth_trans;
    disp(['Calculated: [' num2str(depth_trans(i).') ']'])
    disp(['for comparison: [' num2str(DO_depth(i).') ']'])
    disp( events(i).' )  

    %Plot d18O record and transitions found 
    figure('Visible','off')
    plot(depthI,d18O(k1:k2),'Color',[cNGRIP 0.2])
    hold on 
    plot(depthI,d18OI_sm,'Color',cNGRIP)
    xline(depth_trans,'r')
    scatter(depth_trans,d18OI_sm(idx),10,'o','MarkerFaceColor','w','MarkerEdgeColor','k')
    for m2 = 1:n
    xline(DO_depth,'k') % transition depth found by Buizert et al. 
    end
    grid on 
    axis([depth_lim(1,i) depth_lim(2,i) -50 -30])
    xlabel('NGRIP Depth [m]')
    ylabel(['\delta^{18}O [',char(8240),']'])
    set(gca,'FontName','Times New Roman')
    title(events(i),'FontName','Times New Roman')
end   
