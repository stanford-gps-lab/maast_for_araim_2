function sig2=af_cnmp_259a(del_t,el, freq, const)
%*************************************************************************
%*     Copyright c 2009 The board of trustees of the Leland Stanford     *
%*                      Junior University. All rights reserved.          *
%*     This script file may be distributed and used freely, provided     *
%*     this copyright notice is always kept with it.                     *
%*                                                                       *
%*     Questions and comments should be directed to Todd Walter at:      *
%*     twalter@stanford.edu                                              *
%*************************************************************************
%
%SIG2_AAD calculate airborne psueudorange confidence (variance) per
%DFMC SBAS MOPS ED259A
%SIG2=CNMP_MOPS(DEL_T, EL)
%   DEL_T is the track time in seconds since last cycle slip and is not used
%   EL is the elevation angle in radians
%   SIG2 is the psueudorange confidence (variance) in meters^2 
%
% SEE ALSO INIT_AADA INIT_AADB

%created 12 October, 2007 by Todd Walter
%Modified 22 March, 2022 by Juan Blanch to harmonize with 259A. Code noise
%formulas need to be double-checked

global CNMP_MOPS_A0 CNMP_MOPS_A1 CNMP_MOPS_THETA0
global CNMP_MOPS_B0 CNMP_MOPS_B1 CNMP_MOPS_PHI0 GAMMA_L1_L5
%freq = 0 : GPS L1, freq =1 : dual freq, freq = 2, GPS L5

if freq==1
    %GPS and Galileo L1-L5 iono free
    sig2 = 1/GAMMA_L1_L5*(0.4^2 + (0.34+0.40*exp(-el*180/(pi*14))).^2) ; %scaled by GAMMA_L1_L5 to be consistent with af_cnmp_mops use 

else
    if freq==2
            %GPS and Galileo L5
            sig2 = (0.11 + 0.18*exp(-el/(15*pi/180))).^2 +...
           (CNMP_MOPS_B0 + CNMP_MOPS_B1*exp(-el/CNMP_MOPS_PHI0)).^2;  
    else
      if const == 'GPS'
             %GPS L1
             sig2 = (CNMP_MOPS_A0 + CNMP_MOPS_A1*exp(-el/CNMP_MOPS_THETA0)).^2 +...
            (CNMP_MOPS_B0 + CNMP_MOPS_B1*exp(-el/CNMP_MOPS_PHI0)).^2;
      end

      if const == 'Gal'
             %Galileo L1
             sig2 = (0.13 + 0.17*exp(-el/(13*pi/180))).^2 +...
            (CNMP_MOPS_B0 + CNMP_MOPS_B1*exp(-el/CNMP_MOPS_PHI0)).^2;
      end
    end

end
