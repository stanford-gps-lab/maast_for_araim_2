function [nmax, p_not_monitored]=determine_nmax(p,p_thres)
%*************************************************************************
%*     Copyright c 2013 The board of trustees of the Leland Stanford     *
%*                      Junior University. All rights reserved.          *
%*     This script file may be distributed and used freely, provided     *
%*     this copyright notice is always kept with it.                     *
%*                                                                       *
%*     Questions and comments should be directed to Juan Blanch at:      *
%*     blanch@stanford.edu                                               *
%*************************************************************************

%Created 14 August 2013 by Juan Blanch
%DETERMINE_NMAX determines the maximum size of the subsets to be monitored.

n = length(p);
p_not_monitored = 1;
r = 0;

while ((p_not_monitored>p_thres)&(r<(n+1)))
    r = r+1;
    p_not_monitored = compute_p_not_monitored(p,r);
end

nmax = r-1;

%End

