
function [vpl] = compute_hor_protection_level_v5(sigma1,sigma2, t1,t2, pfault, phmi, pl_tol, alloc_max,vpl_target)
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
if nargin<7
    vpl_target = 0;
end


%find worst case subset that must be monitored
sigmah2 = sigma1.^2+sigma2.^2;
    sigmah =sqrt(sigmah2);
    a = (t1.*sigma1 + t2.*sigma2)./sigmah;
    c = (t1.*sigma2 - t2.*sigma1)./sigmah;
sigmah2_min = min(sigmah2);

%consolidate all modes that are within alpha of smaller sigma;
alpha = 1.3;
idcons =(sigmah2<sigmah2_min*alpha^2);
idlarge = ~idcons;
idcons(1)= boolean(0);

sigma1_cons = max(sigma1(idcons));
sigma2_cons = max(sigma2(idcons));
 
t1_cons = max(t1(idcons));
t2_cons = max(t2(idcons));

sigmah_cons =max(sigmah(idcons));
a_cons = max(a(idcons));
c_cons = max(abs(c(idcons)));

pfault_cons = sum(pfault(idcons));
alloc_max_cons = min(sum(alloc_max(idcons)),1);

%idlarge = find(sigma>sigma_max*alpha);

%consolidated inputs

 sigma1 = [sigma1(1);sigma1_cons;sigma1(idlarge)];
 sigma2 = [sigma2(1);sigma2_cons;sigma2(idlarge)];
 t1 = [t1(1);t1_cons;t1(idlarge)];
 t2 = [t2(1);t2_cons;t2(idlarge)];

 sigmah = [sigmah(1);sigmah_cons;sigmah(idlarge)];
 a      = [a(1);a_cons;a(idlarge)];
 c      = [c(1);c_cons;c(idlarge)];

pfault = [pfault(1); pfault_cons; pfault(idlarge)];
alloc_max = [alloc_max(1);alloc_max_cons; alloc_max(idlarge)];


%%%% Exclude sigmas that are infinite and evaluate their integrity contribution
% index_Inf = find(sigma1.*sigma2 == Inf);
% index_fin = setdiff(1:length(sigma1),index_Inf);
index_Inf = (sigmah==Inf);
index_fin = ~index_Inf;

p_not_monitorable = sum(pfault(index_Inf).*alloc_max(index_Inf));
%alloc = zeros(length(sigma1),1);
if p_not_monitorable>=phmi    
    vpl = Inf;
    %alloc = zeros(length(sigma1),1);
else
    sigmah = sigmah(index_fin);
    %sigma2 = sigma2(index_fin);
    a = a(index_fin);
    c = c(index_fin);
    sigma1 = sigma1(index_fin);
    sigma2 = sigma2(index_fin);
    t1 = t1(index_fin);
    t2 = t2(index_fin);
    % sigmah2 = sigma1.^2+sigma2.^2;
    % sigmah =sqrt(sigmah2);
    % a = (t1.*sigma1 + t2.*sigma2)./sigmah;
    % c = (t1.*sigma2 - t2.*sigma1)./sigmah;

    pfault = pfault(index_fin);
    alloc_max = alloc_max(index_fin);
    phmi = phmi - p_not_monitorable;
    maxCount = 10;    %maximum number of iterations

    %determine lower bound on vpl
    %Klow=-norminv(phmi./pfault);
    Klow=-norminv(min(1,phmi./(pfault.*alloc_max)));
    vpl_low =sqrt(max(max((max(a+Klow.*sigmah,0)).^2-c.^2),0));

    %determine upper bound on vpl
    Khigh = max(0,-norminv(.5*phmi./(length(sigma1)*pfault)));
    vpl_high =sqrt(max((a+Khigh.*sigmah).^2-c.^2));

    %compute logarithm of phmi
    %% 
    log10phmi=log10(phmi);

    if vpl_low>vpl_high     
       warning('vpl_low>vpl_high, check code')
    end
    
 
    
    count=0;
    while ((vpl_high-vpl_low>pl_tol)&&(count<maxCount))&&vpl_high>vpl_target
        count = count+1;
        vpl_half = (vpl_low+vpl_high)/2;
        K = (((vpl_half.^2)+c.^2).^.5 - a)./sigmah;
         
        pl1 = sigma1.*K+t1;
        pl2 = sigma2.*K+t2;
        K1 = (pl1-t1(1))./sigma1(1);
        K2 = (pl2-t2(1))./sigma2(1);
        
        cdfhalf = log10(sum(pfault.*min(modnormcdf(-K)+.5*modnormcdf(-K1)+.5*modnormcdf(-K2),alloc_max)));
                
        if cdfhalf>log10phmi
           vpl_low = vpl_half;
        else
           vpl_high = vpl_half;
        end
    end

   vpl = vpl_high;
  
   %alloc(index_fin) = min(modnormcdf(K),alloc_max);
end
%End