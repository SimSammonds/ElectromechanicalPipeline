function [spike] = ContSpkSpatial(CONdata,CONclst,SPKdata,SPKclst,oldElecConfig,elecConfig,flipConfig)
% ContSpkSpatial: Take a set of temporally correlated contractions and
% spikes, and validate them based on location. Spikes which have an overlap
% of >=15 mm are considered valid. Due to the nature of contractile
% measurements with the EndoFLIP, only vertical SPK position matters.
% This function will also export important SPK metrics.

% CONdata and SPK data must be in a specific structure as generated from
% ContractionMetrics & DataOrganiser, respectively.
% elecConfig & flipConfig are generated by RecentreElecConfigs.m

% Outputs:
% 'spikes' idx: 1=SPKclst, 2=propLength(mm), 3=duration(s), 4=1stATdiff(s),
% 5=overlap(mm), 6=minStartAT(s), 7=maxFinalAT(s)

% Author: Sam Simmonds
% Date: 8th Novemember 2022

spike = []; % Leave blank if not valid

% Contraction propagation distances:
CON_idx = find(CONdata(:,1) == CONclst);
CON_elecs = CONdata(CON_idx,2);
CON_dist = flipConfig(CON_elecs);

% Times for export metrics
CON_AT = min(CONdata(CON_idx,3));

% Spike propagation distances:
SPK_idx = find(SPKdata(:,1) == SPKclst);

SPK_AT = min(SPKdata(SPK_idx,3));
SPK_elecs = SPKdata(SPK_idx,2);
for i = 1:length(SPK_elecs)
    [elec_r, elec_c] = find(oldElecConfig == SPK_elecs(i));
    SPK_dist(i) = elecConfig(elec_r,elec_c);
end

% Find overlap (mm)
min_list = [min(SPK_dist), min(CON_dist)];
max_list = [max(SPK_dist), max(CON_dist)];
overlap = (min(max_list) - max(min_list));

% Does it surpass threshold?
thresh = 10; % (mm)
if overlap >= thresh
    % Calculate additional metrics
    propLen = max(SPK_dist)-min(SPK_dist);
    minStartAT = min(SPKdata(SPK_idx,3));
    maxFinalAT = max(SPKdata(SPK_idx,4));
    dur = maxFinalAT - minStartAT;
    startATdiff = CON_AT - SPK_AT;
    % Output
    spike = [SPKclst,propLen,dur,startATdiff,overlap,minStartAT,maxFinalAT];
end

end
