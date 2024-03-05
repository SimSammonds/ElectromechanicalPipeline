function [gSPKs, iSPKs, idpSPKs, dpSPKs] = independentSPKs(SPKdata,SW_SPK_ML,g_clsts)
% independentSPKs is a function to extract the clusters, and label them as
% GASTRIC, INTESTINAL, or INDEPENDENT

% Author: Sam Simmonds
% Date: 1st September 2023

% Inputs:
% SPK data, processed; SPK_data; 1=SPKclst, 2=elec, 3=StartAT, 4=EndAT
% Associated list; SW_SPK_ML; 1=SWclst, 2=SPKclst, 3=meanStartATdiff, 4=overlap(mm)
    % contains all associated SWs/SPKs, find the inverse to locate the independent list

% Internal functionality:
    % Find the independent SPKs and output the metrics of those SPKs

% Output:
% [category]SPKs; IDX: 1=SPKclst, 2=elec, 3=StartAT, 4=EndAT

% Indentify the independent SPK bursts
idpSPKs = [];
gSPKs = [];
iSPKs = [];
dpSPKs = [];
allSPKs = unique(SPKdata(:,1));

% Loop through each row of SPKdata and allocate the respective row to one
% of the above matrices:
for i = 1:length(allSPKs)
    if ~ismember(SPKdata(i,1),SW_SPK_ML(:,2)) % independent
        idpSPKs = [idpSPKs; SPKdata(i,:)]; % Append
    elseif ismember(SPKdata(i,1),g_clsts) % gastric
        gSPKs = [gSPKs; SPKdata(i,:)]; % Append
        dpSPKs = [dpSPKs; SPKdata(i,:)]; % Append
    else % duodenal
        iSPKs = [iSPKs; SPKdata(i,:)]; % Append
        dpSPKs = [dpSPKs; SPKdata(i,:)]; % Append
    end
end

end