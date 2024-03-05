function [DataOut] = mergeMetrics(SW_SPK,SW_CON,CON_SPK)
% mergeMetrics is a function to format only the useful (correlatable)
% metrics from each script.
% Date: 23rd March 2023
% Inputs: from formatMetrics

DataOut = []; % Establish matrix

% SW_SPKs:
if ~isempty(SW_SPK)
    DataOut(:,1) = SW_SPK(:,2);     % SW dur
    DataOut(:,2) = SW_SPK(:,4);     % SW amp
    DataOut(:,3) = SW_SPK(:,10);    % SPK dur
    DataOut(:,4) = SW_SPK(:,11);    % SPK amp
    DataOut(:,5) = SW_SPK(:,12);    % SPK nrg
end
% SW_CONs:
if ~isempty(SW_CON)
    [r,~] = size(SW_CON);           % Size
    DataOut(1:r,6) = SW_CON(:,2);   % SW dur
    DataOut(1:r,7) = SW_CON(:,4);     % SW amp
    DataOut(1:r,8) = SW_CON(:,10);    % CON dur
    DataOut(1:r,9) = SW_CON(:,11);    % CON amp
    DataOut(1:r,10) = SW_CON(:,12);   % CON nrg
end
% CON_SPKs:
if ~isempty(CON_SPK)
    [r,~] = size(CON_SPK);          % Size
    DataOut(1:r,11) = CON_SPK(:,3);   % CON amp
    DataOut(1:r,12) = CON_SPK(:,4);   % CON nrg
    DataOut(1:r,13) = CON_SPK(:,5);   % CON dur
    DataOut(1:r,14) = CON_SPK(:,12);  % SPK dur
    DataOut(1:r,15) = CON_SPK(:,13);  % SPK amp
    DataOut(1:r,16) = CON_SPK(:,14);  % SPK nrg
end
end







