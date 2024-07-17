function emt_contour(lats, lons, emt,percent)

global GRAPH_EMT_CONTOURS GRAPH_VPL_COLORS

n_lats=length(lats);
n_lons=length(lons);
n_levels=size(GRAPH_EMT_CONTOURS,2);

temp=num2str(GRAPH_EMT_CONTOURS');
ticklabels(n_levels,:)=['> ' temp(n_levels,:)];
for idx=1:n_levels-1
  ticklabels(idx,:)=['< ' temp(idx+1,:)];
end


clf
if nargin > 3
  bartext = ['EMT (m) - ' num2str(percent*100,4) '%'];
else
  bartext = 'EMT (m)';
end
svm_contour(lons,lats,reshape(emt,n_lons,n_lats)', GRAPH_EMT_CONTOURS, ...
            ticklabels, GRAPH_VPL_COLORS, bartext);


title('EMT as a function of user location', 'FontSize', 14);



end

