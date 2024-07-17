%% 
%*************************************************************************
%*     Copyright c 2024 The board of trustees of the Leland Stanford     *
%*                      Junior University. All rights reserved.          *
%*     This script file may be distributed and used freely, provided     *
%*     this copyright notice is always kept with it.                     *
%*                                                                       *
%*     Questions and comments should be directed to Juan Blanch at:      *
%*     blanch@stanford.edu                                               *
%*************************************************************************
%
%MAIN_ARAIM_MAAST runs maast for ARAIM.
% all settings are edited in this file (RAIM ISM parameters, the user grid, the constellation
%configuration, and the receiver settings)
%allows for repeatable and recorded runs
%The parameters are described in:	Blanch, J., Walter, T., Enge, P., Lee, Y., Pervan, B., Rippl, M., Spletter, A., Kropp, V.,
%"Baseline Advanced RAIM User Algorithm and Possible Improvements," IEEE Transactions on Aerospace and Electronic Systems, 
%Volume 51,  No. 1, January 2015.

%In the position optimization option, the code attempts to lower the
%protection levels by applying the method described in: 
% Blanch, J., Walter, T., Enge, P., Kropp, V.,”A Simple Position Estimator that Improves Advanced RAIM Performance,” 
% IEEE Transactions on Aerospace and Electronic Systems Vol. 51, No. 3, July 2015.
%Some of the constants might be set to a different value than in the paper

%The protection levels are computed based on the algorithm described in
%the ARAIM Airborne Algorithm Algorithm Description Document v4.2

%Created 2024 June 14 by Juan Blanch


global ARAIM_USRMASK_GPS ARAIM_USRMASK_GLO ARAIM_USRMASK_GAL ARAIM_USRMASK_BDU 
global ARAIM_SIN_USRMASK_GPS ARAIM_SIN_USRMASK_GLO ARAIM_SIN_USRMASK_GAL
global ARAIM_SIN_USRMASK_BDU 



%Mask angles

ARAIM_USRMASK_GPS = 5;
ARAIM_USRMASK_GLO = 5;
ARAIM_USRMASK_GAL = 5;
ARAIM_USRMASK_BDU = 5;

ARAIM_SIN_USRMASK_GPS = sin(ARAIM_USRMASK_GPS*pi/180);
ARAIM_SIN_USRMASK_GLO = sin(ARAIM_USRMASK_GLO*pi/180);
ARAIM_SIN_USRMASK_GAL = sin(ARAIM_USRMASK_GAL*pi/180);
ARAIM_SIN_USRMASK_BDU = sin(ARAIM_USRMASK_BDU*pi/180);

%%%%% Integrity Support Data parameters (example based on current draft of ICAO Annex 10 default ISD (June 2024) and example broadcast URA for GPS)

isd.ura_gps  = 2.4;
isd.ura_gal  = 6;
isd.ura_glo  = 9;
isd.ura_bdu  = Inf;

isd.bias_gps = 0;% 0.75;%
isd.bias_gal = 0;%0.75;%
isd.bias_glo = 0;%0.75;
isd.bias_bdu = 0;%0.75

isd.ure_gps = 2.4;%0.66*isd.ura_gps%
isd.ure_gal = 4;%0.66*isd.ura_gal%
isd.ure_glo  = 8;
isd.ure_bdu  = Inf;

isd.bias_cont_gps = 0;
isd.bias_cont_gal = 0;
isd.bias_cont_glo = 0;
isd.bias_cont_bdu = 0;


isd.rsat_gps = 1e-5;%1e-5;
isd.rsat_gal = 2e-5;
isd.rsat_glo = 3.4e-5;
isd.rsat_bdu = 1e-5;


isd.psat_gal = 3e-5;%Galileo performance standards specify a state probability

isd.mfd_sat_gps = 1;%(in hours)
isd.mfd_sat_gal = isd.psat_gal/isd.rsat_gal; %Note that Galileo does not guarantee an mean fault duration 9MFD). Here the MFD is a derived parameter that is
%used to account for the specified rate and probability
isd.mfd_sat_glo = 3;
isd.mfd_sat_bdu = 1;

isd.rconst_gps = 1e-8; %For H-ARAIM, PCONST_GPS can be set to 10^-8
isd.rconst_gal = 1e-4;
isd.rconst_glo = 1e-5;
isd.rconst_bdu = 1e-4;

isd.pconst_gal =2e-4; %Galileo performance standards specify a probability and a rate

isd.mfd_const_gps = 1;
isd.mfd_const_gal = isd.pconst_gal/isd.rconst_gal;
isd.mfd_const_glo = 10;
isd.mfd_const_bdu = 1;

%%%%% ARAIM algorithm parameters %%%% 

prm.fde_flag    = 1;%Computes worst case PLs,EMT, sig_acc under a fault or outage scenario. In this version, only single outages and faults are taken into account for the continuity assessment.
%Probability of failed exclusion is set at PFA
prm.fde_wf_flag = 0;%Flag that sets whether constellation out is considered in the exclusion options
prm.pl0_fde     = 0;%When FDE_FLAG is off,and this flag is on the PLS are computed assuming that some of the integrity budget is reserved for the exclusion function 


% Vertical-ARAIM settings (example)

prm.phmi_vert = 9.8e-8;
prm.phmi_hor  = 2e-9;
prm.p_thres   = 6e-8;
prm.p_exc_thres = 1;%For FD, set above the maximum value of any fault mode
prm.pfa_vert  = 3.9e-6;
prm.pfa_hor   = 0.9e-6;
prm.p_emt = 1e-5;
prm.pl_tol= 1e-2;
prm.fc_thres = .01;

prm.sig_acc_max_vert = 1.86; %used in the estimator optimization
prm.sig_acc_max_hor1 = 3;
prm.sig_acc_max_hor2 = 3;

prm.N_es_int  = 25;
prm.N_es_cont = 25;  
prm.t_exp     = 1/24;

prm.vpl_target = 35; %these parameters are used in the fast mode and also in the optimization option
prm.hpl_target = 40;
prm.emt_target = 15;
prm.sig_acc_target = 1.87;

prm.hpl_variant = 0;

% HPL_VARIANT = 0: ARAIM ADD v4.2
% HPL_VARIANT = 1: Blanch ION GNSS+ 2023
% HPL_VARIANT = 2: Racelis ION GNSS+ 2022
prm.max_pmd_flag = 1;%if this flag is on, the integrity risk computation exploits the fact that the integrity risk over
% the exposure interval cannot exceed the probability that the fault is
% present in the first place.

prm.opt = 1;%determines whether we attempt to optimize the position estimation:
% OPT = 0: no optimization
% OPT = 1: optimization

% Horizontal-ARAIM settings (example)
if 1

prm.phmi_vert = 0.1e-8;
prm.phmi_hor  = 9.9e-8;
prm.p_thres   = 9e-8;
prm.p_exc_thres = 9e-7; %if FD results only, set to one
prm.pfa_vert  = 0.01e-7;
prm.pfa_hor   = 5e-7;
prm.p_emt = 1e-5;
prm.pl_tol= 1e-2;
prm.fc_thres = .05;

prm.sig_acc_max_vert = 10;
prm.sig_acc_max_hor1 = 60;
prm.sig_acc_max_hor2 = 60;


prm.N_es_int  = 450;
prm.N_es_cont = 450;  
prm.t_exp     = 1;

prm.vpl_target = Inf; %these parameters are used in the fast mode and also in the optimization option
prm.hpl_target = 185;
prm.emt_target = Inf;
prm.sig_acc_target = Inf;
prm.max_pmd_flag = 0;%if this flag is on, the integrity risk computation exploits the fact that the integrity risk over
% the exposure interval cannot exceed the probability that the fault is
% present in the first place. (Described in Blanch, J. and Walter, T. "Approaches to
%Improve Advanced RAIM Protection Levels" ION GNSS+2023)


end

prm.fast = 1; %This mode uses techniques the computational time of the MAAST run, like the bound on dual faults by one term in the PL equation, the consolidation of fault modes 

%%%%%%%%
% SV Menu
%WGC scenarios

svfile = {'almmops.txt','almanac Galileo 24 Week 703.alm.txt'};
%svfile = {'almmops-1.txt','almanac Galileo 24-1 Week 703.alm.txt'};
%svfile = {'almmops.txt','almanac Galileo 24 Week 703.alm.txt','almglonass.txt'};
%svfile = {'almmops.txt','almanac Galileo 24 + 3 Spare Week 703.alm.txt'};
%svfile ={'almgps24+3.txt','almanac Galileo 24 + 3 Spare Week 703.alm.txt'};
%svfile ={'almgps24+3.txt','almanac Galileo 24 + 8 Spare Week 703.alm.txt'};
%svfile = {'almmops-1.txt'};
%svfile = {'almmops.txt'};
%svfile ={'almgps24+3.txt'};
%svfile = {'almmops_22.txt'};
%svfile ={'alm050418.txt'};
%svfile = 'gps06262020.txt';
%svfile = {'current.txt'};  
%svfile =  {'almmops.txt','almgalileo.txt','almglonass.txt'};
% svfile =  {'almmops.txt','almgalileo24.txt'};

init_const;      % global physical and gps constants
init_col_labels_araim; % column indices 
init_mops;       % MOPS constants

close all;


%User signals:

% dual_freq = 0 : L1 only
% dual_freq = 1 : L1-L5
% dual_freq = 2 : L5 only

 dual_freq = 1;            
      

% USER CNMP Menu

     
      %select ED259A model
      usrcnmpfun = 'af_cnmp_259a';      
      init_cnmp_mops;
      
      %select AAD-A model
%      usrcnmpfun = 'af_cnmpaad';
%      init_aada;
      
      %select AAD-B model
%      usrcnmpfun = 'af_cnmpaad';      
%      init_aadb;
      
      %select AAD-B model
      %usrcnmpfun = 'af_cnmp_mops';  
      

% USER Menu
%select the world as the user area
      usrpolyfile = 'usrworld.dat';
      usrlatstep = 10;
      usrlonstep = 10;
      
%Start time for simulation
      TStart = 0;%20*600;
      
      %End time for simulation
      %TEnd = 862200;
      TEnd = 86400;
      
      % Size of time step
      TStep = 300;
      
%select CONUS as the user area
%      usrpolyfile = 'usrconus.dat';
      
      %select Alaska as the user area
%      usrpolyfile = 'usralaska.dat';
      
      %select Canada as the user area
%      usrpolyfile = 'usrcanada.dat';
      
      %select Mexico as the user area
%      usrpolyfile = 'usrmexico.dat';
      
      %select North America as the user area
      %usrpolyfile = 'usrn_america.dat';
      
      %select Europe as the user area
%      usrpolyfile = 'usreurope.dat';
      
      %select Japan as the user area
%      usrpolyfile = 'usrmsas.dat';
      
      %select Brazil as the user area
%      usrpolyfile = 'usrbrazil.dat';
        % check if file(s) exist
        i=1;
        while i<=size(svfile,2)
          if iscell(svfile)
            fid=fopen(svfile{i});
          else
            fid=fopen(svfile);
            i = size(svfile,2);
          end
          if fid==-1
              fprintf('Almanac file not found.  Please try again.\n');
              return;
          else
              fclose(fid);
          end 
          i=i+1;
        end
        
     
      


% Mode / Alert limit

      %choose PA mode vs NPA  
      pa_mode = 1;
      
      %choose VAL, HAL, EMT Threshold, and sigma accuracy threshold
      %vhal = [35, 40, 15, 1.87];
      vhal = [Inf, 556, Inf, Inf];
      
% OUTPUT Menu

      %initialize histograms
      %init_hist;
        
      
 
      % turn on or off output options
      %1: Availability  2: V/HPL  3: EMT  4: sig_acc  
      
      outputs = [1 1 0 0];

      % Assign percentage
      percent = 0.999; % 1 = 100%
      
      % Coverage values computed for users with |lat|<lamax
      latmax = 90/90*pi/2;
      


% RUN Simulation

profile on
 svmrun_araim(usrcnmpfun, usrpolyfile, svfile, TStart, TEnd, ...
             TStep, usrlatstep, usrlonstep, outputs, percent, vhal, ...
             pa_mode, dual_freq, latmax, prm, isd );%,'outputs.mat');%Adding an input results file will only recompute the hpl, vpl, emt,
 % sig-acc if they don't already meet the target requirements.  This can be
 % used in parametric studies when it is known that there is a monotonic
 % dependence on a parameter.
profile off
