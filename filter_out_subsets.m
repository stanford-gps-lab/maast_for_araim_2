function [sigma, sigma_ss, bias, bias_ss,...
     pfault_inst, pfault_exp, p_not_monitored, idx, x, chi2] = filter_out_subsets(sigma, sigma_ss, bias,...
     bias_ss, pfault_inst, pfault_exp, p_not_monitored, x, chi2, nmodes)
 %*************************************************************************
%*     Copyright c 2012 The board of trustees of the Leland Stanford     *
%*                      Junior University. All rights reserved.          *
%*     This script file may be distributed and used freely, provided     *
%*     this copyright notice is always kept with it.                     *
%*                                                                       *
%*     Questions and comments should be directed to Juan Blanch at:      *
%*     blanch@stanford.edu                                               *
%*************************************************************************

%Created 14 August 2012 by Juan Blanch
 
 
 
 idx = find(min(sigma(:,1:3),[],2)<Inf);
 sigma = sigma(idx,:);
 sigma_ss = sigma_ss(idx,:);      
 bias = bias(idx,:);
 bias_ss = bias_ss(idx,:);
 p_not_monitored = p_not_monitored +sum(pfault_exp)-sum(pfault_exp(idx,:));
 pfault_inst =pfault_inst(idx,:);
 pfault_exp =pfault_exp(idx,:);
 x = x(idx,:);
 chi2 = chi2(idx,:);
 nmodes =nmodes(idx,:);