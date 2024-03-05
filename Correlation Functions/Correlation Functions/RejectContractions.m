function [CONdata4] = RejectContractions(CONdata,EXPdata,time)
%% OUTDATED - USE ContractionMetrics.m INSTEAD
% RejectContractions is a function to confirm that contractions manually
% marked are suitable. Valid contractions must have at least 2 electrodes
% which are >= 1 mm amplitude threshold.
% CONdata4 has the same form as CONdata, but omits location labels.

CONdata2 = [];
thresh = 1; % Threshold for valid contraction (mm)
for i = 1:length(CONdata(:,1))
    % Check each electrode for amplitude - reject if below threshold:
    % Collect indices
    startTime = find(time == CONdata(i,3));
    endTime = find(time == CONdata(i,5));
    elec = CONdata(i,2);
    
    % Extract range of relevant data for this electrode
    cont = EXPdata(startTime:endTime,elec);
    ampl = max(cont) - min(cont);
    
    if ampl >= thresh
        % Save details to new array
        CONdata2(end+1,:) = CONdata(i,1:6); % Omit location
    end
    
end

% Check that contraction propagates at least 2 electrodes
CONdata3 = [];
clsts = unique(CONdata2(:,1));
for i = 1:length(clsts)
    cur_clst = clsts(i);
    clstIdx = find(cur_clst == CONdata2(:,1));
    if numel(clstIdx) >= 2
        % Save clsts to a 3rd version of the array
        CONdata3(end+1:end+numel(clstIdx),1) = CONdata2(clstIdx,1);
    end
end

CONdata3 = unique(CONdata3);

% Take list of valid contractions and resave data
CONdata4 = [];
for i = 1:length(CONdata)
    if ismember(CONdata(i,1),CONdata3)
        % If valid, save
        CONdata4(end+1,:) = CONdata(i,1:6);
    end
end


end

