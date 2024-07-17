function  [subsets, pfault_inst, pfault_exp, p_not_monitored] = determine_subsets_v5(G, p_sat_inst, p_const_inst, tau_sat, tau_const, p_thres, fc_thres,p_exc_thres, phmi)
%*************************************************************************
%*     Copyright c 2024 The board of trustees of the Leland Stanford     *
%*                      Junior University. All rights reserved.          *
%*     This script file may be distributed and used freely, provided     *
%*     this copyright notice is always kept with it.                     *
%*                                                                       *
%*     Questions and comments should be directed to Juan Blanch at:      *
%*     blanch@stanford.edu                                               *
%*************************************************************************

%Created september 2014 by Juan Blanch
%Updated april 2016 by Juan Blanch (integrates changes proposed in WGC ARAIM briefing
%from April 2016)
%Modified November 2019 to improve treatment of temporal exposure effects
%tau_sat and tau_const is equal to the ratio of the exposure window to the mean time to
%detect
%updated March 2024 by Juan Blanch

if nargin<8
    p_exc_thres = 1;
end
if nargin<9
    phmi = 1e-7;
end

Nsat = size(G,1);
Nconst = size(G,2)-3;

if max(size(p_sat_inst))==1              %if p_sat is a scalar, apply to all satellites
   p_sat_inst = ones(Nsat,1)*p_sat_inst;
end
if max(size(p_const_inst))==1
   p_const_inst = ones(Nconst,1)*p_const_inst;
end

%Remove constellations for which there are no satellites
idc = find(sum(abs(G(:,4:end))>0));
G = G(:,[1:3 3+idc]);
Nconst = size(G,2)-3;
p_const_inst = p_const_inst(idc);

%Compute primary fault probabilities over exposure window
p_sat_exp = p_sat_inst.*(1+tau_sat);
p_const_exp = p_const_inst.*(1+tau_const);

%Compute probability of no fault over exposure window
pnofault_exp = prod(1-p_sat_exp)*prod(1-p_const_exp);
pnofault_inst = prod(1-p_sat_inst)*prod(1-p_const_inst);

% %Initialize pnotmonitored
% p_not_monitored = 1-pnofault_exp;

%Determine upper bound of subset size
[nsatconstmax, ~] = determine_nmax([p_sat_exp;p_const_exp],p_thres);
subsetsize =0;
nsatconstmax =max(nsatconstmax,2);
%Determine upperbound of number of subsets
for j=0:nsatconstmax
    subsetsize = subsetsize + nchoosek(Nconst+Nsat,j);
end

%Allocate memory for subsets
subsets_ex = zeros(subsetsize,Nsat+Nconst);
subsets_ex(1,:) = zeros(1,Nsat+Nconst); %subset corresponding to the all-in-view

%Initialize probabilities
pfault_inst = zeros(subsetsize,1);
pfault_inst(1) = 1;

pfault_exp = zeros(subsetsize,1);
pfault_exp(1) = pnofault_exp;

%Form initial list of subsets

j = 1;

for k=1:nsatconstmax
    subsets_k = determine_k_subsets(Nsat+Nconst,k);

    pfault_inst_k = prod((subsets_k*diag([p_sat_inst./(1-p_sat_inst);p_const_inst./(1-p_const_inst)])+ (1 - subsets_k)),2)*pnofault_inst;
    pfault_exp_k = prod((subsets_k*diag([p_sat_exp./(1-p_sat_exp);p_const_exp./(1-p_const_exp)])+ (1 - subsets_k)),2)*pnofault_exp;

    subsets_ex(j+1:j+size(subsets_k,1),:) = subsets_k;
    
    pfault_inst(j+1:j+size(subsets_k,1))= pfault_inst_k;
    pfault_exp(j+1:j+size(subsets_k,1))= pfault_exp_k;

    j = j+size(subsets_k,1);
end

subsets_const = (G(:,4:(3+Nconst))*subsets_ex(:,Nsat+1:Nsat+Nconst)')';
subsets_sat  = min(subsets_ex(:,1:Nsat) + subsets_const ,1);

idconst = find(sum(subsets_ex(:,end-length(idc)+1:end))>0);


%Consolidate low probability faults in high probability constellation faults
for  jj=1:Nconst
    %find subset corresponding to constellation wide fault first
    id_const_c=find((sum(subsets_ex(:,1:Nsat),2)==0).*...
        (sum(subsets_ex,2)==1).*(subsets_ex(:,Nsat+idconst(jj))==1));
   
    %find subsets that include a constellation fault and satellites whithin
    %that constellation
    
    index_cs =find((((subsets_sat*subsets_sat(id_const_c,:)'))>0).*...
    (((subsets_sat*(1-subsets_sat(id_const_c,:))'))<=0));
   
    idremove = find(pfault_inst(index_cs)<fc_thres*pfault_inst(id_const_c));
    idrmv    = index_cs(idremove);

    %remove from the list and add probability to constellation wide fault
    pfault_inst(id_const_c) =pfault_inst(id_const_c) + sum(pfault_inst(idrmv));
    pfault_exp(id_const_c) =pfault_exp(id_const_c) + sum(pfault_exp(idrmv));
    idnew = setdiff(1:length(pfault_inst),idrmv);
    
    pfault_inst = pfault_inst(idnew);
    pfault_exp  = pfault_exp(idnew);
    subsets_sat = subsets_sat(idnew,:);
    subsets_ex  = subsets_ex(idnew,:);
end
nsubsets = size(subsets_sat,1);
%From this initial list, remove potentially weaker subsets 
sats_out = sum(subsets_sat,2);

%order list by the number of satellites out
[sats_out_s, index_s]=sort(sats_out,1,'ascend');

% [~, index_s]=sort(Nsat_min,2,'descend');
 pfault_exp_s = pfault_exp(index_s);
 pfault_inst_s = pfault_inst(index_s);
 subsets_sat_s = subsets_sat(index_s,:);

%Identify subsets corresponding to valid exclusion options
%sats_in = Nsat - sats_out_s;
exclusion_subsets = subsets_sat(pfault_inst>p_exc_thres,:);
Nsat_min = min(((1-subsets_sat_s)*(1-exclusion_subsets)')')';
 
%Find all modes that exceed total integrity budget.  These must be
%monitored
idlarge = find(pfault_exp_s>phmi);
p_not_monitored = 1-sum(pfault_exp_s(idlarge));

idexc = find((Nsat_min==0)); %second condition ensures that we exclude from the list modes that cannot be monitored under an exclusion option 


idother = setdiff(1:nsubsets,[idlarge;idexc]);
h = 0;
n_other = length(idother);

while (h<n_other)&&(p_not_monitored>p_thres)
        h=h+1;
        if pfault_exp_s(idother(h))>0 
         p_not_monitored = p_not_monitored - pfault_exp_s(idother(h));
        end
end

idother = idother(1:h);
%pfault_inst, pfault_exp

subsets =1 - subsets_sat_s([idlarge' idother],:);
pfault_inst = pfault_inst_s([idlarge' idother]);
pfault_exp = pfault_exp_s([idlarge' idother]);

if p_not_monitored>1e-7
    
    p_not_monitored
    
end
