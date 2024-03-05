function [elecs] = PassProkElecs(toappData,i,first,last)
% Quick function to pull useful elecs from previous PROKINETICS file
% i == 1 for SWs
% i == 2 for SPKs

% Author: Sam Simmonds
% Date: 31st May 2023

if i == 1 % SW
    Data = toappData.toapp.TimeAmplCluster(1,1:end-1);
    [~, Data] = DataOrganiser(Data,1);
    elecs = unique(Data(:,2));
    
elseif i == 2 % SPK
    Data = toappData.toapp.TimeAmplSpikeCluster(first:last,1);
    [~, Data] = DataOrganiser(Data,2);
    elecs = unique(Data(:,2));
end

end