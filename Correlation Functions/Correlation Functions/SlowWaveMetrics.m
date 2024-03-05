function [metrics] = SlowWaveMetrics(SWclst,SWdata,oldElecConfig,elecConfig)
% SlowWaveMetrics is a function to extract useful data for a specific slow wave.
% metrics idx:
% 1=SWclst, 2=duration, 3=VRTpropLength
% Also manually find SW amplitude - GEMS export

% Author: Sam Simmonds
% Date: 14th Novemember 2022

% Find indices
idx = find(SWdata(:,1) == SWclst);

% Temp arrays
AT_temp = [];
dist_temp = [];

for i = idx' % For each elec this SW appears in
    
    % Spatial
    elec = SWdata(i,2);
    [r,c] = find(oldElecConfig == elec);
    dist = elecConfig(r,c);
    dist_temp = [dist_temp, dist];
    
    % Temporal
    AT = SWdata(i,3);
    AT_temp = [AT_temp, AT];
    
end

% Export
duration = max(AT_temp) - min(AT_temp);
prop = abs(max(dist_temp) - min(dist_temp));
metrics = [SWclst, duration, prop];

end

