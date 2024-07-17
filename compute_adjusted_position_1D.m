function [ s0 ] = compute_adjusted_position_1D(pap_subset,sigma,sigma_ss, svec, sigpr2_acc,sig_acc_max)
%*************************************************************************
%*     Copyright c 2015 The board of trustees of the Leland Stanford     *
%*                      Junior University. All rights reserved.          *
%*     This script file may be distributed and used freely, provided     *
%*     this copyright notice is always kept with it.                     *
%*                                                                       *
%*     Questions and comments should be directed to Juan Blanch at:      *
%*     blanch@stanford.edu                                               *
%*************************************************************************


%Computes estimator according to paper: Blanch, J., Walter, T., Enge, P., Kropp, V.,”A Simple Position Estimator
% that Improves Advanced RAIM Performance,” Accepted for publication in IEEE Transactions in Aerospace Electronic Systems.

%Created by Juan Blanch March 2015

idx_large_ap = find((pap_subset>=1e-11)&(sigma<Inf));
[~, idx_max_sigma] = max(sigma(idx_large_ap));         
idx_weak = idx_large_ap(idx_max_sigma);        

sigma_acc_best = sqrt((svec.*svec)*sigpr2_acc);
sig_acc_best = sigma_acc_best(1);

a = sigma_ss(idx_weak)^2;
b = sum(svec(1,:).*(svec(idx_weak,:)-svec(1,:)).*sigpr2_acc');
c = sig_acc_best^2 - (sig_acc_max-1e-5)^2;

if c<0
    t = (-b+sqrt(b^2 - a*c))/a;
    t = min(t,1);
else
    t = 0;
end

s0 = svec(1,:)+t*(svec(idx_weak,:)-svec(1,:));
end

