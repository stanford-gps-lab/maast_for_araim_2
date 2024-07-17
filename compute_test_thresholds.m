function [T1, T2, T3] = compute_test_thresholds(sigma_ss, bias_ss, pfa_vert, pfa_hor,nmodes)

%*************************************************************************
%*     Copyright c 2012 The board of trustees of the Leland Stanford     *
%*                      Junior University. All rights reserved.          *
%*     This script file may be distributed and used freely, provided     *
%*     this copyright notice is always kept with it.                     *
%*                                                                       *
%*     Questions and comments should be directed to Juan Blanch at:      *
%*     blanch@stanford.edu                                               *
%*************************************************************************

%This function computes the test thresholds
%See MHSS_RAIM_BASELINE

nsets = size(sigma_ss,1);
% nsets1 = max(1,sum(sigma_ss(:,1)<Inf));
% nsets2 = max(1,sum(sigma_ss(:,1)<Inf));
% nsets3 = max(1,sum(sigma_ss(:,1)<Inf));

nsets1 = max(1,sum(nmodes(sigma_ss(:,1)<Inf)));
nsets2 = max(1,sum(nmodes(sigma_ss(:,1)<Inf)));
nsets3 = max(1,sum(nmodes(sigma_ss(:,1)<Inf)));


if nsets>1

Kfa_1 = - norminv(.25*pfa_hor/(nsets1-1));
Kfa_2 = - norminv(.25*pfa_hor/(nsets2-1));
Kfa_3 = - norminv(.5*pfa_vert/(nsets3-1));

T1 = Kfa_1 * sigma_ss(:,1) + bias_ss(:,1);
T2 = Kfa_2 * sigma_ss(:,2) + bias_ss(:,2);
T3 = Kfa_3 * sigma_ss(:,3) + bias_ss(:,3);
else
T1 = 0;
T2 = 0;
T3 = 0;
end
    
%add delta to prevent numerical indeterminations. This delta does not
%impact integrity

T1 = T1 + eps;
T2 = T2 + eps;
T3 = T3 + eps;


