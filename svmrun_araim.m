function svmrun_araim(usrcnmpfun, usrfile, svfile, tstart, tend, ...
				tstep, usrlatstep, usrlonstep, outputs, percent, vhal, ...
                pa_mode, dual_freq, latmax, prm, isd, outputs_com_filename)
%*************************************************************************
%*     Copyright c 2007 The board of trustees of the Leland Stanford     *
%*                      Junior University. All rights reserved.          *
%*     This script file may be distributed and used freely, provided     *
%*     this copyright notice is always kept with it.                     *
%*                                                                       *
%*     Questions and comments should be directed to Todd Walter at:      *
%*     twalter@stanford.edu                                              *
%*************************************************************************
%
% SVMRUN    Run the SVM simulation.
%   svmrun(gpsudrefun, geoudrefun, givefun, usrtrpfun, usrcnmpfun,...
%                wrstrpfun, wrsgpscnmpfun, wrsgeocnmpfun,...
%                wrsfile, usrfile, igpfile, svfile, tstart, tend, tstep,...
%                usrlatstep, usrlonstep, outputs, percent);
% Inputs:
%   gpsudrefun    -   function to run for udre calculation for GPS satellites
%   geoudrefun    -   function to run for udre calculation for GEO satellites
%   givefun       -   function to run for give calculation
%   usrcnmpfun    -   function to run for cnmp delay calculation for user
%   wrsgpscnmpfun -   function to run for cnmp delay calculation for wrs   
%   wrsfile       -   file containing wrs position data
%   usrfile       -   file containing user position boundary polygon
%   igpfile       -   file containing IGP mask points
%   svfile        -   Yuma almanac file if TStep not zero
%                     Static satellite position file if TStep is zero
%   geodata       -   matrix with geo information
%   tstart        -   start time of simulation (for almanac option)
%   tend          -   end time of simulation (for almanac option)
%   tstep         -   time step of simulation (for almanac option),
%                     should be 0 for static satellite position option
%   userlatstep   -   latitude spacing of user grid
%   userlonstep   -   longitude spacing of user grid
%   outputs       -   array of ON-OFF flags (1 for ON, 0 for OFF) for 
%                     output options, corresponding to:
%                     1) availability     2) udre map     3) give map
%                     4) udre histogram   5) give histogram 6) V/HPL
%                     7) coverage/availability
%   percent       -   percent value to use for calculating availability, 
%                     give map, and/or V/HPL
%   vhal          -   VAL / HAL to use calculating availability
%   pa_mode       -   Whether to calulate vertical and horizontal or 
%                     horizontal only
%   dual_freq     -   Whether or not to calculate GIVE or have a dual
%                     frequency user

%Modified Todd Walter June 28, 2007 to include VAL, HAL and PA vs. NPA mode
%Modified Juan Blanch August 14, 2013 to adapt it ARAIM simulations.
%Modified Juan Blanch June 6, 2024

global CONST_H_IONO;
global COL_SAT_UDREI COL_SAT_PRN COL_SAT_XYZ COL_SAT_MINMON ...
%         COL_USR_UID COL_IGP_LL COL_IGP_DELAY...
%         COL_U2S_UID COL_U2S_PRN COL_U2S_AZ COL_U2S_MAX COL_U2S_TTRACK0 COL_U2S_GENUB...
%         COL_U2S_IPPLL COL_U2S_EL COL_IGP_GIVEI COL_IGP_MINMON COL_IGP_BETA...
%         COL_IGP_CHI2RATIO COL_U2S_SIGACC COL_U2S_SIGNOM
global  COL_USR_XYZ COL_USR_LL COL_USR_LLH COL_USR_EHAT ...
        COL_USR_NHAT COL_USR_UHAT COL_USR_INBND COL_USR_MAX
% global HIST_UDRE_NBINS HIST_GIVE_NBINS HIST_UDRE_EDGES HIST_GIVE_EDGES
% global HIST_UDREI_NBINS HIST_GIVEI_NBINS HIST_UDREI_EDGES HIST_GIVEI_EDGES
global MOPS_SIN_USRMASK MOPS_SIN_WRSMASK MOPS_NOT_MONITORED
% global MOPS_UDREI_NM MOPS_GIVEI_NM MOPS_UDREI_MAX MOPS_GIVEI_MAX
global CNMP_TL3
% global GRAPH_MOVIE_FIGNO TRUTH_FLAG RTR_FLAG
global UDREI_CONST GEOUDREI_CONST GIVEI_CONST


if nargin==17
   load(outputs_com_filename,'vpl','hpl','emt','sig_acc'); 
end     
% GEO Position Menu

      geodata = [];



fprintf('initializing run\n');
alm_param = read_yuma(svfile);


% initialize sat, wrs, usr, igp matrices (structures)
% see documentation for format of SATDATA,WRSDATA,USRDATA,IGPDATA,
%   WRS2SATDATA, USR2SATDATA
satdata=[];
satdata = init_satdata(geodata, alm_param, satdata, 0, 0);
ngps = size(alm_param,1);
ngeo = size(geodata,1);
nsat=ngps+ngeo;
[usrdata,usrlatgrid,usrlongrid] = init_usrdata(usrfile,usrlatstep,usrlonstep);


%wrsdata = init_wrsdata(wrsfile);
%nwrs=size(wrsdata,1);
%[igpdata, inv_igp_mask] = init_igpdata(igpfile);
%wrs2satdata = init_usr2satdata(wrsdata,satdata);
usr2satdata = init_usr2satdata(usrdata,satdata);
%truth_data = [];

%wrstrpfun = 'af_trpmops';
usrtrpfun = 'af_trpmops';
if ~isempty(which('init_trop_osp'))
    init_trop_osp();
    %wrstrpfun = 'af_trpadd';
end

if tstep>0
    tend = tend - 1; % e.g. 86400 becomes 86399
    ntstep = floor((tend-tstart)/tstep)+1;
else
    ntstep = 1;
end
if nargin<17
vpl = repmat(Inf,size(usrdata,1),ntstep);
hpl = repmat(Inf,size(usrdata,1),ntstep);
emt = repmat(Inf,size(usrdata,1),ntstep);
sig_acc =  repmat(Inf,size(usrdata,1),ntstep);
end
%ncrit = vpl;


nusr = size(usrdata,1);

sat_xyz=[];
% if 0
% % find all los rise times for cnmp calculation, 
% % start from tstart-CNMP_TL3 (below this, cnmp is at floor value)
% if isempty(CNMP_TL3)
%     CNMP_TL3 = 12000;
% end
% wrs2sat_trise = find_trise(tstart-CNMP_TL3,tend,MOPS_SIN_WRSMASK,alm_param,...
%             wrsdata(:,COL_USR_XYZ),wrsdata(:,COL_USR_EHAT),...
%             wrsdata(:,COL_USR_NHAT),wrsdata(:,COL_USR_UHAT));
% %add blank rows for geos
% 
% nrise=size(wrs2sat_trise,2);
% wrs2sat_trise = reshape(wrs2sat_trise, ngps, nwrs, nrise);
% wrs2sat_trise(ngps+1:ngps+ngeo,:,:)=zeros(ngeo, nwrs, nrise);
% wrs2sat_trise = reshape(wrs2sat_trise, nsat*nwrs, nrise);
% end

itstep = 1;
tcurr = tstart;
while tcurr<=tend
    tic
    % get current satellite positions
    satdata = init_satdata(geodata, alm_param, satdata, tcurr, 0);
    
   
    sat_xyz = [sat_xyz; satdata(:,COL_SAT_XYZ)];
    vhpl_comp(:,1) = vpl(:,itstep); 
    vhpl_comp(:,2) = emt(:,itstep); 
    vhpl_comp(:,3) = sig_acc(:,itstep); 
    vhpl_comp(:,4) = hpl(:,itstep); 
    % USER processing
    [vhpl, usr2satdata] = usrprocess_araim(satdata,usrdata,...
                        usr2satdata,usrtrpfun,...
                        usrcnmpfun,...
                        alm_param,tcurr,pa_mode,dual_freq, prm, isd, vhpl_comp);
    vpl(:,itstep) = vhpl(:,1);
    emt(:,itstep) = vhpl(:,2);
    sig_acc(:,itstep) = vhpl(:,3);
    hpl(:,itstep) = vhpl(:,4);
   
        


 
    % update time
    if tstep==0
        break;
    else
        fprintf('Time: %d / %d done\n',itstep,ntstep);
        tcurr = tcurr+tstep;
        itstep = itstep+1;
    end
    toc
end


save 'outputs' satdata usrdata sat_xyz vpl hpl emt sig_acc usrlatgrid usrlongrid svfile ;


plot_no_gui_araim(usrdata,vpl,hpl,usrlatgrid,usrlongrid,outputs,percent,vhal,pa_mode,...
            emt,sig_acc, latmax)




