function [newElecConfig,newFlipConfig] = RecentreElecConfigs(elecConfig,pyl_elec,pyl_flip)
% RecentreElecConfigs takes the GEMS elecConfig file, and the EndoFLIP
% catheter sensors and recentres around the pylorus. This will allow for
% spatial correlations between the two. Electrode configurations should be
% exported in mm.

% Author: Sam Simmonds
% Date: 8th Novemember 2022

% Initialise:
newElecConfig = [];
newFlipConfig = [];

% Start with GEMS electrodes:
[r,c] = size(elecConfig);
for i = 1:r
    % Check that row has at least one real electrode (not 257)
    if (numel(unique(elecConfig(i,:))) > 1) % At least one valid elec
        for j = 1:c
            % Wipe value if = 257
            if elecConfig(i,j) == 257
                newElecConfig(i,j) = NaN;
            else
            % For each element, assign a distance value from pyl
            % Positive == gastric   % Negative == duodenal
                newElecConfig(i,j) = pyl_elec - i;
            end
        end
    end % Else ignore data
end
newElecConfig = newElecConfig*5;

% EndoFLIP electrodes:
newFlipConfig = [1:16] - pyl_flip;
newFlipConfig = -newFlipConfig';
newFlipConfig = newFlipConfig*10;

end

