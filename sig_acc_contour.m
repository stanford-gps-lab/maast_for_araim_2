function sig_acc_contour(lats, lons, sig_acc,percent)

global GRAPH_SIG_ACC_CONTOURS GRAPH_VPL_COLORS

n_lats=length(lats);
n_lons=length(lons);
n_levels=size(GRAPH_SIG_ACC_CONTOURS,2);

temp=num2str(GRAPH_SIG_ACC_CONTOURS',3);
ticklabels(n_levels,:)=['> ' temp(n_levels,:)];
for idx=1:n_levels-1
  ticklabels(idx,:)=['< ' temp(idx+1,:)];
end


clf
if nargin > 3
  bartext = ['\sigma_a_c_c (m) - ' num2str(percent*100,4) '%'];
else
  bartext = '\sigma_a_c_c (m)';
end
svm_contour(lons,lats,reshape(sig_acc,n_lons,n_lats)', GRAPH_SIG_ACC_CONTOURS, ...
            ticklabels, GRAPH_VPL_COLORS, bartext);


title('\sigma_a_c_c as a function of user location', 'FontSize', 14);



end

