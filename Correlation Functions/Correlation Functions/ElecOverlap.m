function [numElecs] = ElecOverlap(SWdata,SWclst,SPKdata,SPKclst)
% ElecOverlap: Finds the number of overlapping electrodes for specified slow wave and
% spike clusters.
% SW and SPK data must be in a specific structure as generated from
% DataOrganiser function.

% Author: Sam Simmonds
% Date: 8th Novemember 2022

% Find function to gather indices of interest
SW_idx = find(SWdata(:,1) == SWclst);
SPK_idx = find(SPKdata(:,1) == SPKclst);

% Use indices to create a list of electrodes
SW_elec_list = SWdata(SW_idx,2);
SPK_elec_list = SPKdata(SPK_idx,2);

numElecs = numel(intersect(SW_elec_list,SPK_elec_list));

end

