function [firstATstart,firstATend] = CONtiming(CONdata,FLIPconfig)
% CONtiming() is a function to discern if contraction timings are
% co-ordinated between the stomach, pylorus, and duodenum.

% Kolmogorov–Smirnov test can check against a distribution - uniform may be
% the easiest
% Shapiro-wilks + Q-Q plot can check for normality

% Hypothesis: If KS test can reject uniformly distributed intervals
% & SW / QQ cannot reject normality, the data is likely temporally
% correlated around the peak of the histogram

% Null-hypothesis: If we cannot reject KS test, and can reject SW/QQ, then
% we cannot say that the timings are correlated
 
% FLIPconfig is zeroed around pylorus
% Create a list of all the GAS/PYL elecs, to be ismember() compared. 
GP_list = FLIPconfig(FLIPconfig>=0);
GP_CONs = [];
D_CONs = [];

% For each gastric/pyloric contraction:
clsts = unique(CONdata(:,1));
for clst = clsts'
    % Find idx
    idx = find(CONdata(:,1) == clst);
    elecs = CONdata(idx,2);
    
    if numel(intersect(elecs,GP_list)) > 0  % Gastric/Pyl contraction
        startAT = min(CONdata(idx,3));  % Timings
        endAT = max(CONdata(idx,5));
        GP_CONs = [GP_CONs; clst, startAT, endAT]; % Append

    elseif numel(intersect(elecs,GP_list)) < 1 % Duodenal only contraction
        startAT = min(CONdata(idx,3));  % Timings
        endAT = max(CONdata(idx,5));
        D_CONs = [D_CONs; clst, startAT, endAT]; % Append

    end
end

for i = 1:length(GP_CONs(:,1)) % For each gastric/pyl contraction
    GP_CON = GP_CONs(i,1); % Find label

    % Finding timing differences for both start and end times
    startATdiffs = D_CONs(:,2) - GP_CONs(i,2);
    endATdiffs = D_CONs(:,2) - GP_CONs(i,3);
    startATdiffs2 = sort(startATdiffs(startATdiffs>=0),'ascend');
    endATdiffs2 = sort(endATdiffs(endATdiffs>=0),'ascend');

    if ~isempty(startATdiffs2)
        firstATstart(i) = startATdiffs2(1); % Closest in each group
    end
    if ~isempty(endATdiffs2)
        firstATend(i) = endATdiffs2(1);
    end

end

firstATstart = firstATstart';
firstATend = firstATend';

% (defined as any contraction which touches a single gas/pyl elec)
% Extract the timing (start & finish)
% Find all duodenal contractions (start) (do not touch gas/pyl elec)
% Find interval
% Plot as QQ
% Perform SW test
% Craft a uniform distribution across the average interval

end