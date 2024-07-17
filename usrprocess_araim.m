function [vhpl,usr2satdata] = usrprocess_araim(satdata,usrdata,usr2satdata,usrtrpfun,...                                                     usrcnmpfun, alm_param,time, pa_mode,dual_freq, prm, isd, vhpl_comp)%*************************************************************************%*     Copyright c 2024 The board of trustees of the Leland Stanford     *%*                      Junior University. All rights reserved.          *%*     This script file may be distributed and used freely, provided     *%*     this copyright notice is always kept with it.                     *%*                                                                       *%*     Questions and comments should be directed to Juan Blanch at:      *%*     blanch@stanford.edu                                               *%*************************************************************************%%% User Processing% Created by Juan Blanch 14 August 2013 global MOPS_SIN_USRMASK CONST_H_IONO global COL_SAT_PRN COL_SAT_XYZ COL_USR_XYZ COL_USR_LL COL_SAT_UDREI ...        COL_USR_EHAT COL_USR_NHAT COL_USR_UHAT  ...        COL_U2S_UID COL_U2S_PRN COL_U2S_LOSXYZ COL_U2S_GXYZB ...        COL_U2S_LOSENU COL_U2S_GENUB COL_U2S_EL COL_U2S_AZ ...        COL_U2S_SIG2MON COL_U2S_IPPLL COL_U2S_IPPXYZ COL_U2S_TTRACK0 COL_U2S_INITNAN...global MOPS_MIN_GPSPRN MOPS_MAX_GPSPRN MOPS_MIN_GLOPRN MOPS_MAX_GLOPRN ...       MOPS_MIN_GALPRN MOPS_MAX_GALPRN MOPS_MIN_GEOPRN MOPS_MAX_GEOPRN ...       MOPS_MIN_BDUPRN MOPS_MAX_BDUPRNglobal ARAIM_SIN_USRMASK_GPS ARAIM_SIN_USRMASK_GLO ARAIM_SIN_USRMASK_GAL ...       ARAIM_SIN_USRMASK_BDUglobal GAMMA_L1_L5 CONST_GAMMA_L1L5global COL_U2S_SIGACC COL_U2S_SIGNOM    nsat = size(satdata,1);nusr = size(usrdata,1);nlos = nsat*nusr;% initialize some values to NaNusr2satdata(:,COL_U2S_INITNAN) = NaN;% form los data from usr to satellitesusr2satdata(:,COL_U2S_GXYZB) = find_los_xyzb(usrdata(:,COL_USR_XYZ), ...                                            satdata(:,COL_SAT_XYZ));usr2satdata(:,COL_U2S_GENUB) = find_los_enub(usr2satdata(:,COL_U2S_GXYZB),...   usrdata(:,COL_USR_EHAT),usrdata(:,COL_USR_NHAT),usrdata(:,COL_USR_UHAT));%Find satellites above GPS and Galileo mask angleidusrgps=find(usr2satdata(:,COL_U2S_PRN) >= MOPS_MIN_GPSPRN & ...          usr2satdata(:,COL_U2S_PRN) <= MOPS_MAX_GPSPRN);idusrglo=find(usr2satdata(:,COL_U2S_PRN) >= MOPS_MIN_GLOPRN & ...          usr2satdata(:,COL_U2S_PRN) <= MOPS_MAX_GLOPRN);idusrgal=find(usr2satdata(:,COL_U2S_PRN) >= MOPS_MIN_GALPRN & ...          usr2satdata(:,COL_U2S_PRN) <= MOPS_MAX_GALPRN);idusrbdu=find(usr2satdata(:,COL_U2S_PRN) >= MOPS_MIN_BDUPRN & ...          usr2satdata(:,COL_U2S_PRN) <= MOPS_MAX_BDUPRN);abv_mask_gps = find(-usr2satdata(idusrgps,COL_U2S_GENUB(3)) >= ARAIM_SIN_USRMASK_GPS);abv_mask_glo = find(-usr2satdata(idusrglo,COL_U2S_GENUB(3)) >= ARAIM_SIN_USRMASK_GLO);abv_mask_gal = find(-usr2satdata(idusrgal,COL_U2S_GENUB(3)) >= ARAIM_SIN_USRMASK_GAL);abv_mask_bdu = find(-usr2satdata(idusrbdu,COL_U2S_GENUB(3)) >= ARAIM_SIN_USRMASK_BDU);abv_mask = [idusrgps(abv_mask_gps);idusrglo(abv_mask_glo);idusrgal(abv_mask_gal);idusrbdu(abv_mask_bdu)];%abv_mask = find(-usr2satdata(:,COL_U2S_GENUB(3)) >= MOPS_SIN_USRMASK);if(~isempty(abv_mask))  [usr2satdata(abv_mask,COL_U2S_EL),usr2satdata(abv_mask,COL_U2S_AZ)] = ...        find_elaz(usr2satdata(abv_mask,COL_U2S_LOSENU));  usr2satdata(abv_mask,COL_U2S_IPPLL) = find_ll_ipp(usrdata(:,COL_USR_LL),...                                usr2satdata(:,COL_U2S_EL),...                                usr2satdata(:,COL_U2S_AZ), abv_mask);endidxold = find(~isnan(usr2satdata(:,COL_U2S_TTRACK0)));idxnew = setdiff(abv_mask,idxold);% set start time of track for lost los's to NaNusr2satdata(setdiff(idxold,abv_mask),COL_U2S_TTRACK0) = NaN; % lost losusr2satdata(idxnew,COL_U2S_TTRACK0) = time; % new los% tropo    sig2_trop = feval(usrtrpfun,usr2satdata(:,COL_U2S_EL));% cnmp% if ~isempty(usrcnmpfun)%     sig2_cnmp = feval(usrcnmpfun,time-usr2satdata(:,COL_U2S_TTRACK0),...%                                 usr2satdata(:,COL_U2S_EL));% end% initialize outputsvhpl = repmat(NaN,nusr,4);[t1 t2]=meshgrid(1:nusr,1:nsat);usridx=reshape(t1,nlos,1);satidx=reshape(t2,nlos,1);los_xyzb = usr2satdata(:,COL_U2S_GXYZB);los_enub = usr2satdata(:,COL_U2S_GENUB);el = usr2satdata(:,COL_U2S_EL);if(~isempty(abv_mask))        el = el(abv_mask);    prn = usr2satdata(abv_mask,COL_U2S_PRN);    indgps=find(prn >= MOPS_MIN_GPSPRN & prn <= MOPS_MAX_GPSPRN);    indglo=find(prn >= MOPS_MIN_GLOPRN & prn <= MOPS_MAX_GLOPRN);    indgal=find(prn >= MOPS_MIN_GALPRN & prn <= MOPS_MAX_GALPRN);    indbdu=find(prn >= MOPS_MIN_BDUPRN & prn <= MOPS_MAX_BDUPRN);      % Dual frequency measurement noise for integrity       ura     = ones(length(abv_mask),1)*Inf;       sig2_if = ones(length(abv_mask),1)*Inf;       bnom_i  = ones(length(abv_mask),1)*Inf;       rsat_i  = ones(length(abv_mask),1);       msat_i  = ones(length(abv_mask),1);              rsat_i(indgps) = isd.rsat_gps;        rsat_i(indglo) = isd.rsat_glo;       rsat_i(indgal) = isd.rsat_gal;       rsat_i(indbdu) = isd.rsat_bdu;              msat_i(indgps) = isd.mfd_sat_gps;       msat_i(indglo) = isd.mfd_sat_glo;       msat_i(indgal) = isd.mfd_sat_gal;       msat_i(indbdu) = isd.mfd_sat_bdu;              ura(indgps) = isd.ura_gps;       ura(indglo) = isd.ura_glo;       ura(indgal) = isd.ura_gal;       ura(indbdu) = isd.ura_bdu;              bnom_i(indgps) = isd.bias_gps;       bnom_i(indglo) = isd.bias_glo;       bnom_i(indgal) = isd.bias_gal;       bnom_i(indbdu) = isd.bias_bdu;             sig_trv = 0.12*1.001./sqrt(0.002001+sin(el).^2);       sig_uire_df = 40./(261.0+(180/pi*el).^2)+0.018;  if dual_freq ==1            sig2_if(indgps) = GAMMA_L1_L5*feval(usrcnmpfun,NaN,el(indgps), dual_freq);            sig2_if(indglo) = GAMMA_L1_L5*feval(usrcnmpfun,NaN,el(indglo), dual_freq); %use GPS cnmp characterization             sig2_if(indbdu) = GAMMA_L1_L5*feval(usrcnmpfun,NaN,el(indbdu), dual_freq);            sig2_if(indgal) = GAMMA_L1_L5*feval(usrcnmpfun,NaN, el(indgal),dual_freq);            sig2_if = sig2_if + sig_uire_df.^2;         else            mag_lat = usr2satdata(abv_mask,COL_U2S_IPPLL(1)) + ...               0.064*180*cos((usr2satdata(abv_mask,COL_U2S_IPPLL(2))/180-1.617)*pi);                       %mid-latitude klobuchar confidence            sig2_uire(abv_mask) = 20.25*obliquity2(el);            %low-latitude klobuchar confidence            idx = find(abs(mag_lat) < 20);            if(~isempty(idx))                sig2_uire(abv_mask(idx)) = 81*obliquity2(el(idx));            end            %high-latitude klobuchar confidence            idx = find(abs(mag_lat) > 55);            if(~isempty(idx))                sig2_uire(abv_mask(idx)) = 36*obliquity2(el(idx));            end                        sig2_uire = sig2_uire(abv_mask)';            if dual_freq ==0                 sig2_if(indgps) = sig2_uire(indgps) + feval(usrcnmpfun,NaN,el(indgps),dual_freq,'GPS');%af_cnmp_mops(NaN,el(indgps));       sig2_if(indgal) = sig2_uire(indgal) + feval(usrcnmpfun,NaN, el(indgal),dual_freq,'Gal');%af_cnmp_galileo(el(indgal));       sig2_if(indglo) = sig2_uire(indglo) + feval(usrcnmpfun,NaN,el(indglo),dual_freq,'GPS');%af_cnmp_mops(NaN,el(indglo)); %use GPS cnmp characterization        sig2_if(indbdu) = sig2_uire(indbdu) + feval(usrcnmpfun,NaN,el(indbdu),dual_freq,'GPS');%af_cnmp_mops(NaN,el(indbdu));            else       sig2_if(indgps) = CONST_GAMMA_L1L5^2*sig2_uire(indgps) + feval(usrcnmpfun,NaN,el(indgps), dual_freq);%af_cnmp_mops(NaN,el(indgps));       sig2_if(indgal) = CONST_GAMMA_L1L5^2*sig2_uire(indgal) + feval(usrcnmpfun,NaN, el(indgal), dual_freq);%af_cnmp_galileo(el(indgal));       sig2_if(indglo) = CONST_GAMMA_L1L5^2*sig2_uire(indglo) + feval(usrcnmpfun,NaN,el(indglo), dual_freq);%af_cnmp_mops(NaN,el(indglo)); %use GPS cnmp characterization        sig2_if(indbdu) = CONST_GAMMA_L1L5^2*sig2_uire(indbdu) + feval(usrcnmpfun,NaN,el(indbdu), dual_freq);%af_cnmp_mops(NaN,el(indbdu)); %use GPS cnmp characterization       end  end                     sig2    = ura.^2+ sig_trv.^2+ sig2_if;        usr2satdata(abv_mask,COL_U2S_SIGNOM) = sqrt(sig2);       % Dual frequency measurement noise for accuracy              ure = ones(length(abv_mask),1)*Inf;       bcont_i = ones(length(abv_mask),1)*Inf;       ure(indgps) = isd.ure_gps;       ure(indglo) = isd.ure_glo;       ure(indgal) = isd.ure_gal;       ure(indbdu) = isd.ure_bdu;              bcont_i(indgps) = isd.bias_cont_gps;       bcont_i(indglo) = isd.bias_cont_glo;       bcont_i(indgal) = isd.bias_cont_gal;       bcont_i(indbdu) = isd.bias_cont_bdu;       sig2acc = ure.^2+ sig_trv.^2+ sig2_if;              usr2satdata(abv_mask,COL_U2S_SIGACC) = sqrt(sig2acc);             vhpl(1:max(usridx(abv_mask)),:)=usr_vhpl_araim(los_enub(abv_mask,:), ...                             usridx(abv_mask), sig2, ...                             usr2satdata(abv_mask,COL_U2S_PRN),...                             sig2acc, bnom_i, bcont_i, rsat_i, msat_i,prm, isd, vhpl_comp(1:max(usridx(abv_mask)),:));      bad_usr=find(vhpl(:,1) <= 0 | vhpl(:,4) <= 0);      if(~isempty(bad_usr))        vhpl(bad_usr,:)=NaN;      end  end