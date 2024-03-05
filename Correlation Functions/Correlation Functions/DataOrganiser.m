function [numClst,OrgData] = DataOrganiser(data,datatype)
% DataOrganiser is a function to organise messy data as derived from GEMS'
% 'toapp' data. The output will be organised as follows:
% For slow wave data:
    % cols: 1=SWclst, 2=elec, 3=AT
% For spike data:
    % cols: 1=SPKclst, 2=elec, 3=StartAT, 4=EndAT
% datatype: 1=SWs, 2=SPKs

% Author: Sam Simmonds
% Date: 8th Novemember 2022

% Initialise parameters:
count = 0;
len = length(data);
if datatype == 2 % SPIKES
    OrgData = zeros(10,4);
    for i = 1:len
        [entries,~] = size(data{i,1});
        if entries ~= 0
            for j = 1:entries
                count = count + 1;
                OrgData(count,1) = i; % Cluster
                OrgData(count,2) = data{i,1}(j,4); % Elec
                OrgData(count,3) = data{i,1}(j,3); % startTime
                OrgData(count,4) = data{i,1}(j,8); % endTime
            end
        end
    end
    
elseif datatype == 1 % SWs
    OrgData = zeros(100,3);
    for i = 1:len
        [entries,~] = size(data{1,i});
        if entries ~= 0
            for j = 1:entries
                count = count + 1;
                OrgData(count,1) = i; % Cluster
                OrgData(count,2) = data{1,i}(j,4); % Elec
                OrgData(count,3) = data{1,i}(j,3); % AT
            end
        end
    end
end

numClst = numel(unique(OrgData(:,1)));

end

