function [valid_spikes] = ATOverlap(SWdata,SWclst,SPKdata)
% ATOverlap: Create a list of valid spikes (overlapping ATs) for specified
% slow wave cluster.
% SW and SPK data must be in a specific structure as generated from
% DataOrganiser function. Spike are temporally correlated when their 
% activation time is within 2 seconds of SW AT.
% valid_spikes has the form: 1=SPKclst, 2=meanAT, 3=meanATdiff.

% Author: Sam Simmonds
% Date: 8th Novemember 2022

% Find function to gather indices of interest
SW_idx = find(SWdata(:,1) == SWclst);

% Use indices to create a list of electrodes & find mean SW AT
SW_AT_list = SWdata(SW_idx,3);
SW_AT = mean(SW_AT_list);

% Initialise variables:
valid_spikes = [];
temp = [];
% Loop through each spike
for i = 1:(length(SPKdata)-1)
    % Gather cluster and starting AT info
    cur_clst = SPKdata(i,1);
    nxt_clst = SPKdata(i+1,1);
    cur_AT = SPKdata(i,3); % grab current AT
    temp = [temp, cur_AT]; % Save to array for averaging
    
    if (cur_clst ~= nxt_clst)
        % Check mean AT diff is <= 2 seconds
        if (abs(mean(temp)-SW_AT) <= 2)
            valid_spikes(end+1,1:3) = [cur_clst,mean(temp),mean(temp)-SW_AT];
        end
        temp = []; % Reset temp array
    elseif i == (length(SPKdata)-1) % Last entry, manually save
        nxt_AT = SPKdata(i+1,3); % grab next AT
        temp = [temp, nxt_AT]; % Save to array for averaging
        if (abs(mean(temp)-SW_AT) <= 2)
            valid_spikes(end+1,1:3) = [cur_clst,mean(temp),mean(temp)-SW_AT];
        end
        temp = [];
    end % ELSE; do nothing and continue loop
end

end