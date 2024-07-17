%*************************************************************************
%*     Copyright c 2024 The board of trustees of the Leland Stanford     *
%*                      Junior University. All rights reserved.          *
%*     This script file may be distributed and used freely, provided     *
%*     this copyright notice is always kept with it.                     *
%*                                                                       *
%*     Questions and comments should be directed to Juan Blanch at:      *
%*     blanch@stanford.edu                                               *
%*************************************************************************
%uses plot_no_gui

clear; close; clc;
init_graph();
init_col_labels();

filename ='outputs';
load(filename)

%choose PA mode vs NPA  
pa_mode = 1;
      
%choose VAL, HAL, EMT Threshold, and sigma accuracy threshold
vhal = [35, 40, 15, 1.8700]; %LPV-200
%vhal = [50, 40, Inf, Inf]; %LPV
%vhal = [Inf, 185, Inf, Inf]; %RNP0.1  Note that the VPL and HPL computed
%depend on the choice of H-ARAIM and V-ARAIM.
%vhal = [Inf, 556, Inf, Inf]; %RNP0.3

% turn on or off output options
%1: Availability  2: V/HPL  3: EMT  4: sig_acc  
      
outputs = [1 1 0 0];

% Assign percentage
percent = 0.999; % 1 = 100%
      
% Coverage values computed for users with |lat|<latmax
latmax = 90/90*pi/2;
      
plot_no_gui_araim(usrdata,vpl,hpl,usrlatgrid,usrlongrid,outputs,percent,vhal,pa_mode,...
            emt,sig_acc, latmax)      


