function [OUT] = ProContProc(CCL_P,conData,FLIP_tStart)
%ProContProc : Prokinetic Contraction Processing
% This script will take in CCL_P.xlsx files and determine:
% Mean +/- SD for AMP, DUR, & Propagation Length

% INPUTS: ContractionClusters, EndoFLIP Data, time offset
% OUTPUTS: OUT = [DUR(s), PROP(mm), AMP(MEAN; mm), AMP(SD; mm)]
%                > Each row is a new contraction

% Author: Sam Simmonds
% Date: 31st May 2023

% IDX: [clst, elec, startTime, ~, endTime, ~]
CONs = unique(CCL_P(:,1));
OUT = [];
count = 1;

for CON = CONs'                     % For each contraction
    idx = find(CCL_P(:,1) == CON);  % idx
    
    earliest = min(CCL_P(idx,3));   % DUR
    latest = max(CCL_P(idx,5));
    DUR = latest - earliest;        % s
    
    PROP = 10*length(idx);          % PROP; mm
    
    AMP = [];                       % AMP
    % Finding BDF zero:
    if CON <= FLIP_tStart(count,1)
        BDF_zero = FLIP_tStart(count,2);
    else
        count = count + 1;
        BDF_zero = FLIP_tStart(count,2);
    end
    
    Fs_FLIP = 10;                   % Hz
    for i = idx
        CON_st = BDF_zero + (Fs_FLIP*CCL_P(i,3)); % start time
        CON_et = BDF_zero + (Fs_FLIP*CCL_P(i,5)); % end time
        elec = CCL_P(i,2);
        minCONdia = min(conData(CON_st:CON_et,3+elec));
        maxCONdia = max(conData(CON_st:CON_et,3+elec));
        AMP = [AMP, maxCONdia-minCONdia]; % mm
    end
    % For each contraction, average over all electrodes
    MEAN = mean(AMP);
    SD = std(AMP);
    
    OUT = [OUT; DUR, PROP, MEAN, SD];
end


end

