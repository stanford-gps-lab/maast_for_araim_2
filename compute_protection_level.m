
function [vpl, alloc] = compute_protection_level(sigma, threshold_plus_bias, pfault, phmi, pl_tol, alloc_max)
%*************************************************************************
%*     Copyright c 2024 The board of trustees of the Leland Stanford     *
%*                      Junior University. All rights reserved.          *
%*     This script file may be distributed and used freely, provided     *
%*     this copyright notice is always kept with it.                     *
%*                                                                       *
%*     Questions and comments should be directed to Juan Blanch at:      *
%*     blanch@stanford.edu                                               *
%*************************************************************************

%Created 14 August 2012 by Juan Blanch
%Updated 24 February 2015 by Juan Blanch
%COMPUTE_PROTECTION_LEVEL_V2 solves the equation:   pfault.*normcdf((threshold_plus_bias - vpl)./sigma) = phmi
%alloc is the vector of corresponding allocations normcdf((threshold_plus_bias - vpl)./sigma) 
%alloc max is the maximum necessary allocation

if nargin<6
    alloc_max = ones(length(sigma),1);
end

%%%% Exclude sigmas that are infinite and evaluate their integrity contribution
index_Inf = find(sigma == Inf);
index_fin = setdiff(1:length(sigma),index_Inf);

p_not_monitorable = sum(pfault(index_Inf).*alloc_max(index_Inf));
alloc = zeros(length(sigma),1);
if p_not_monitorable>=phmi    
    vpl = Inf;
    alloc = zeros(length(sigma),1);
else
    sigma = sigma(index_fin);
    threshold_plus_bias = threshold_plus_bias(index_fin);
    pfault = pfault(index_fin);
    alloc_max = alloc_max(index_fin);
    phmi = phmi - p_not_monitorable;
    maxCount = 10;    %maximum number of iterations

    %determine lower bound on vpl
    %Klow=-norminv(phmi./pfault);
    Klow=-norminv(min(1,phmi./(pfault.*alloc_max)));
    vpl_low=max(threshold_plus_bias + Klow.*sigma);

    %determine upper bound on vpl
    Khigh = max(0,-norminv(phmi./(length(sigma)*pfault)));
    vpl_high=max(threshold_plus_bias + Khigh.*sigma);

    %compute logarithm of phmi
    log10phmi=log10(phmi);

    if vpl_low>vpl_high     
       warning('vpl_low>vpl_high, check code')
    end
    
    
    count=0;
    while ((vpl_high-vpl_low>pl_tol)&&(count<maxCount))
        count = count+1;
        vpl_half = 0.5*(vpl_low+vpl_high);
        cdfhalf = log10(sum(pfault.*min(modnormcdf((threshold_plus_bias-vpl_half)./sigma),alloc_max)));
        if cdfhalf>log10phmi
           vpl_low = vpl_half;
        else
           vpl_high = vpl_half;
        end
    end

   vpl = vpl_high;
   alloc(index_fin) = min(modnormcdf((threshold_plus_bias-vpl)./sigma),alloc_max);
end
%End