function p_not_monitored = compute_p_not_monitored(p_event,r);
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
%COMPUTE_P_NOT_MONITORED returns an upper bound of the probability of r or more
% of the events with probability p

%p_event is a vector containing the probability of n independent events

%Find events with probability one

idcrtn = find(p_event==1);
ncrtn = length(idcrtn);

idcmp =setdiff(1:length(p_event),idcrtn);

p_event = p_event(idcmp);

r = r-ncrtn;

%Compute probability of no fault
p_no_fault = prod(1-p_event);

if r<=0
    p_not_monitored = 1;
end

if r==1
    p_not_monitored = 1-p_no_fault;
end

if r==2
    p_not_monitored = 1-p_no_fault*(1+sum(p_event./(1-p_event)));
end

if r==3
    p_not_monitored = 1-p_no_fault*(1+sum(p_event./(1-p_event)))...
     - .5*p_no_fault*((sum(p_event./(1-p_event)))^2 - sum((p_event./(1-p_event)).^2))   ;
end

if r>=4 
    p_not_monitored = ((sum(p_event))^r)/prod(1:r);
end
    