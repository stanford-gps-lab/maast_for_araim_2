function covg = coverage(vpl, hpl, emt, acc, val, hal, emt_th, acc_th, avlb, usrlatgrid, usrlongrid,latmax)

%*************************************************************************
%*     Copyright c 2013 The board of trustees of the Leland Stanford     *
%*                      Junior University. All rights reserved.          *
%*     This script file may be distributed and used freely, provided     *
%*     this copyright notice is always kept with it.                     *
%*                                                                       *
%*     Questions and comments should be directed to Juan Blanch at:      *
%*     blanch@stanford.edu                                               *
%*************************************************************************
%COVERAGE calculates the percentage of the users on the grid USRLATGRID
%USRLONGRID between -LATMAX and +LATMAX that meet VPL<VAL, EMT<EMT_TH, ACC<ACC_TH
%with an AVLB availability. Users are weighed by the size of the cell 
%centered on the grid point.

%Created by Juan Blanch August 2013

if nargin<10
    latmax = pi/2;
end


A = size(vpl);
Nloc = A(1);
Nsteps = A(2);

nlat = length(usrlatgrid);
nlon = length(usrlongrid);

lat   = ones(nlon,nlat)* diag(usrlatgrid*pi/180);
lat   = reshape(lat,Nloc,1);

% [latmesh,lonmesh] = meshgrid(usrlatgrid,usrlongrid);
% lat   =latmesh(:)*pi/180;

clat = cos(lat);
i = clat<cos(latmax);
clat(i)=0;

av = sum(((vpl<=val)'.*(emt<=emt_th)'.*(acc<=acc_th)'.*(hpl<=hal)'))/Nsteps;
covg = 100*sum((av>=avlb).*clat')/sum(clat);
