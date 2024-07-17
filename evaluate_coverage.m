%*************************************************************************
%*     Copyright c 2013 The board of trustees of the Leland Stanford     *
%*                      Junior University. All rights reserved.          *
%*     This script file may be distributed and used freely, provided     *
%*     this copyright notice is always kept with it.                     *
%*                                                                       *
%*     Questions and comments should be directed to Juan Blanch at:      *
%*     blanch@stanford.edu                                               *
%*************************************************************************
%EVALUATE_COVERAGE applies the function COVERAGE for the ARAIM users included in
%OUTPUTS.MAT under different thresholds.

%Created by Juan Blanch August 2013

% load outputs_HARAIM_1em5 % _1em3_grps%outputs_baseline1
load outputs
%load outputs_1em3_baseline
val = 35;
hal = 40;
emt_th= 15;
acc_th = 1.87;
avail =.995;% .995;

% val = 50;
% hal = 40;
% emt_th= Inf;
% acc_th = Inf;
% avail = .995;

% val = Inf;
% hal = 185;%;
% emt_th= Inf;
% acc_th = Inf;
% avail = 1;%.995;
% 
% 
latmax = 90/90*pi/2;

covrg_RAIM995all  = coverage(vpl, hpl, emt, sig_acc, val, hal, emt_th, acc_th, avail, usrlatgrid, usrlongrid,latmax);

covrg_RAIM995_vpl = coverage(vpl, hpl, emt, sig_acc, val, hal, Inf,    Inf,    avail, usrlatgrid, usrlongrid,latmax);

covrg_RAIM995_emt = coverage(vpl, hpl, emt, sig_acc, Inf, hal, emt_th, Inf,    avail, usrlatgrid, usrlongrid,latmax);

covrg_RAIM995_acc = coverage(vpl, hpl, emt, sig_acc, Inf, hal, Inf,    acc_th, avail, usrlatgrid, usrlongrid,latmax);

[covrg_RAIM995all covrg_RAIM995_vpl covrg_RAIM995_emt covrg_RAIM995_acc]