function emt = compute_emt(T3, pap_subset, p_emt)


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


%Effective Monitor Threshold formulation proposed at the WG-C ARAIM
%subgroup on June 2012.  p_emt is set at 10^-7

%Find all modes with a probability larger than p_emt 

index_emt = (pap_subset>=p_emt);
emt = max(T3(index_emt));
