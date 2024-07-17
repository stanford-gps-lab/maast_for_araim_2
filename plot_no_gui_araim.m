function plot_no_gui_araim(usrdata,usrvpl,usrhpl,latgrid,longrid,outputs,percent,vhal,pa_mode,emt,sig_acc, latmax)
%*************************************************************************
%*     Copyright c 2024 The board of trustees of the Leland Stanford     *
%*                      Junior University. All rights reserved.          *
%*     This script file may be distributed and used freely, provided     *
%*     this copyright notice is always kept with it.                     *
%*                                                                       *
%*     Questions and comments should be directed to Juan Blanch at:      *
%*     blanch@stanford.edu                                               *
%*************************************************************************
global COL_USR_LL COL_USR_INBND COL_IGP_LL

global GRAPH_AVAIL_FIGNO GRAPH_VPL_FIGNO GRAPH_HPL_FIGNO

                   
                   
                   
init_graph();
init_col_labels();

usrlat = usrdata(:,COL_USR_LL(1));
usrlon = usrdata(:,COL_USR_LL(2));
idx = find(usrlon>=180);
usrlon(idx) = usrlon(idx)-360;

inbnd = usrdata(:,COL_USR_INBND);
inbndidx=find(inbnd);
%igp_mask = igpdata(:,COL_IGP_LL);

nt = size(usrvpl,2);

%Availability
if outputs(1)
    h=figure(GRAPH_AVAIL_FIGNO);
    avail_contour4(latgrid,longrid,usrvpl,usrhpl,emt, sig_acc, inbnd,percent,vhal,pa_mode, latmax);
    set(h,'name','AVAILABILITY CONTOUR');
end

%VPL, HPL
if outputs(2)
    % sort v/hpl for each user and determine vpl at given percentage
    nusr = size(usrdata,1);
    sortvpl = zeros(size(usrvpl));
    sorthpl = zeros(size(usrhpl));
    percentidx = ceil(percent*nt);
    for i = 1:nusr
        sortvpl(i,:) = sort(usrvpl(i,:));
        sorthpl(i,:) = sort(usrhpl(i,:));
    end
    vpl = sortvpl(:,percentidx);
    hpl = sorthpl(:,percentidx);
    %VAL specific plot
    if(pa_mode)
        h=figure(GRAPH_VPL_FIGNO);
        vpl_contour(latgrid,longrid,vpl,percent);
        set(h,'name','VPL CONTOUR');
%        plot(wrsdata(:,COL_USR_LL(2)),wrsdata(:,COL_USR_LL(1)),'m*');
%        plot(usrdata(inbndidx,COL_USR_LL(2)),usrdata(inbndidx,COL_USR_LL(1)),'ko');
%        text(longrid(1)+1,latgrid(1)+2,'* - WRS');
%        text(longrid(1)+1,latgrid(1)+1,'o -  USER');
    end
    h=figure(GRAPH_HPL_FIGNO);
    hpl_contour(latgrid,longrid,hpl,percent);
    set(h,'name','HPL CONTOUR');
end

%EMT
if outputs(3)
    % sort emt for each user and determine emt at given percentage
    nusr = size(usrdata,1);
    sortemt = zeros(size(emt));
    percentidx = ceil(percent*nt);
    for i = 1:nusr
        sortemt(i,:) = sort(emt(i,:));
    end
    emt = sortemt(:,percentidx);
    GRAPH_EMT_FIGNO = 10;
    h=figure(GRAPH_EMT_FIGNO);
    emt_contour(latgrid,longrid,emt,percent);
    set(h,'name','EMT CONTOUR');
end

%Accuracy
if outputs(4)
    % sort sig_acc for each user and determine sig_acc at given percentage
    nusr = size(usrdata,1);
    sortsig_acc = zeros(size(sig_acc));
    percentidx = ceil(percent*nt);
    for i = 1:nusr
        sortsig_acc(i,:) = sort(sig_acc(i,:));
    end
    sig_acc = sortsig_acc(:,percentidx);
    GRAPH_SIGACC_FIGNO = 15;
    h=figure(GRAPH_SIGACC_FIGNO);
    sig_acc_contour(latgrid,longrid,sig_acc,percent);
    set(h,'name','SIG_ACC CONTOUR');
    
    
end

end