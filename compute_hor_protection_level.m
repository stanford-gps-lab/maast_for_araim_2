
function [vpl, alloc] = compute_hor_protection_level(sigma1,sigma2, t1,t2, pfault, phmi, pl_tol, alloc_max)
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
%Updated June 2023 by Juan Blanch
%COMPUTE_HOR_PROTECTION_LEVEL implements the HPL calculation described in Blanch, J. and Walter, T. "Approaches to
%Improve Advanced RAIM Protection Levels" ION GNSS+2023 

if nargin<6
    alloc_max = ones(length(sigma1),1);
end

%%%% Exclude sigmas that are infinite and evaluate their integrity contribution
index_Inf = find(sigma1.*sigma2 == Inf);
index_fin = setdiff(1:length(sigma1),index_Inf);
p_not_monitorable = sum(pfault(index_Inf).*alloc_max(index_Inf));
alloc = zeros(length(sigma1),1);
if p_not_monitorable>=phmi    
    vpl = Inf;
    alloc = zeros(length(sigma1),1);
else
    sigma1 = sigma1(index_fin);
    sigma2 = sigma2(index_fin);
    t1 = t1(index_fin);
    t2 = t2(index_fin);
    sigmah2 = sigma1.^2+sigma2.^2;
    sigmah =sqrt(sigmah2);
    a = (t1.*sigma1 + t2.*sigma2)./sigmah;
    c = (t1.*sigma2 - t2.*sigma1)./sigmah;

    pfault = pfault(index_fin);
    alloc_max = alloc_max(index_fin);
    phmi = phmi - p_not_monitorable;
    maxCount = 10;    %maximum number of iterations

    %determine lower bound on vpl
    
    Klow=-norminv(min(1,phmi./(pfault.*alloc_max)));  
   
    vpl_low =sqrt(max(max((max(a+Klow.*sigmah,0)).^2-c.^2),0));
    %vpl_low = 0;
    %determine upper bound on vpl
    Khigh = max(0,-norminv(.5*phmi./(length(sigma1)*pfault)));
    vpl_high =sqrt(max((a+Khigh.*sigmah).^2-c.^2));

    %compute logarithm of phmi
    log10phmi=log10(phmi);

    if vpl_low>vpl_high     
       warning('vpl_low>vpl_high, check code')
    end
    
    
    count=0;
    while ((vpl_high-vpl_low>pl_tol)&&(count<maxCount))
        count = count+1;
        vpl_half = (vpl_low+vpl_high)/2;
                 
         K = (((vpl_half.^2)+c.^2).^.5 - a)./sigmah;
         pl1 = sigma1.*K+t1;
         pl2 = sigma2.*K+t2;
         K1 = (pl1-t1(1))./sigma1(1);
         K2 = (pl2-t2(1))./sigma2(1);
        %K3 = min([K1 K2],[],2);
        cdfhalf = log10(sum(pfault.*min(modnormcdf(-K)+.5*modnormcdf(-K1)+.5*modnormcdf(-K2),alloc_max)));
       %cdfhalf = log10(sum(pfault.*min(modnormcdf(-K)+modnormcdf(-K3),alloc_max)));


        if cdfhalf>log10phmi
           vpl_low = vpl_half;
        else
           vpl_high = vpl_half;
        end
    end

   vpl = vpl_high;
   alloc(index_fin) = min(modnormcdf(K),alloc_max);
end
