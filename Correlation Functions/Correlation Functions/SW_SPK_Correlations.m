function [ML] = SW_SPK_Correlations(SWdata,SPKdata,threshold)
%SW_SPK_Correlations is a function which takes in a slow wave and spike
%data, and determine which clusters are temporally and spatially
%correlated.

% Outputs: ML (master list)
% IDX: [SWclst, SPKclst, meanStartATdiff(s), overlap(mm), tally]

% Author: Sam Simmonds
% Date: 14th November 2022

ML = [];
if isempty(SWdata) || isempty(SPKdata)
    return
end

SWclsts = unique(SWdata(:,1));
SPKclsts = unique(SPKdata(:,1));

for SW = SWclsts' % For each slow wave

    % Find idx, elec, AT data
    SWidx = find(SWdata(:,1) == SW);
    SWelecs = SWdata(SWidx,2);
    SW_ATs = SWdata(SWidx,3);
    %[r,c] = find(oldElecConfig == SWdata(SWidx,2));
    
    for SPK = SPKclsts' % For each spike
        
        % Establish or reset tally;
        tally = 0;
        % Establish temp array for averaging AT differences
        AT_temp = [];
        
        % Find SPK idx, elecs, AT data
        SPKidx = find(SPKdata(:,1) == SPK);
        SPKelecs = SPKdata(SPKidx,2);
        SPK_startATs = SPKdata(SPKidx,3);
        SPK_endATs = SPKdata(SPKidx,4);
        
        % Find number of shared electrodes
        sharedElec = ismember(SWelecs,SPKelecs);
        numSharedElec = sum(sharedElec);
        
        if (5*numSharedElec) >= threshold(2) % Check overlap threshold
            
            for elec = SWelecs(sharedElec)' % for each shared electrode
                
                % Find AT difference at this electrode
                curSWidx = find(SWdata(SWidx,2) == elec);
                curSPKidx = find(SPKdata(SPKidx,2) == elec);
                SW_AT = SWdata(SWidx(curSWidx),3);
                SPK_startAT = SPKdata(SPKidx(curSPKidx),3);
                SPK_endAT = SPKdata(SPKidx(curSPKidx),4);
                startATdiff = SPK_startAT - SW_AT; % Pos = SW 1st
                
                if ((startATdiff) <= threshold(1) && (startATdiff >= -0.5))
                    tally = tally + 1; % Add to tally
                    AT_temp = [AT_temp, startATdiff];
                end % if
                
            end % for
            
            % Check tally against threshold
            if tally >= threshold(3)
                % Save to master list
                ML(end+1,:) = [SW,SPK,mean(AT_temp),5*numSharedElec,tally];
            end % if
            
        end % if
    end % for
end % for

end % function

