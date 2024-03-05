function [ML] = CON_SPK_Correlations(SPKdata,CONdata,elecConfig,oldElecConfig,flipConfig,threshold,BDFstartTime)
% CON_SPK_Correlations is a function to find all valid correlating
% contractions and spikes.

% Outputs: ML (master list)
% IDX: [CONclst, SPKclst, meanStartATdiff(s), overlap(mm), tally]

% Author: Sam Simmonds
% Date: 14th November 2022

ML = [];
if isempty(CONdata) || isempty(SPKdata)
    return
end

SPKclsts = unique(SPKdata(:,1));
CONclsts = unique(CONdata(:,1));

for SPK = SPKclsts' % For each spike

    % Find idx, elec, AT data
    SPKidx = find(SPKdata(:,1) == SPK);
    SPKelecs = SPKdata(SPKidx,2);
    
    % SPK distances
    for j = 1:length(SPKelecs) % For each electrode
        [elec_r, elec_c] = find(oldElecConfig == SPKelecs(j)); % Find co-ords
        SPKdist(j) = elecConfig(elec_r,elec_c); % Get distance
    end % for
    
    for CON = CONclsts' % For each contraction
        
        % Establish or reset tally;
        tally = 0;
        % Establish temp array for averaging AT differences
        AT_temp = [];
        
        % Find CON idx, elecs, AT data
        CONidx = find(CONdata(:,1) == CON);
        CONelecs = CONdata(CONidx,2);
        CONdist = flipConfig(CONelecs);
        
        % Find overlap (mm)
        min_list = [min(SPKdist), min(CONdist)];
        max_list = [max(SPKdist), max(CONdist)];
        overlap = (min(max_list) - max(min_list));
        
        if overlap >= threshold(8) % Check overlap threshold

            sharedElec = []; % IDX for overlapping SPK elecs
            % Which electrodes overlap?
            for i = 1:length(SPKelecs) % For each electrode
                if (SPKdist(i) <= max(CONdist)) && (SPKdist(i) >= min(CONdist))
                    sharedElec = [sharedElec, i];
                end % if
            end % for
            
            for elec = SPKelecs(sharedElec)' % for each overlapping SPK elec
                
                % Find AT difference at this electrode
                curSPKidx = find(SPKdata(SPKidx,2) == elec);
                SPK_AT = SPKdata(SPKidx(curSPKidx),3);
                [elec_r, elec_c] = find(oldElecConfig == elec); % Find co-ords
                curSPKdist = elecConfig(elec_r,elec_c); % Get distance
                remainder = rem(curSPKdist,10); % Check remainder
                
                if remainder ~= 0 % Then we must take an average of adjacent CON_ATs
                    % Proximal
                    CONelec = find(flipConfig == curSPKdist+5);
                    CONdistIdx = find(CONdata(CONidx,2) == CONelec);
                    CON_proxAT = CONdata(CONidx(CONdistIdx),3);
                    
                    % Distal
                    CONelec = find(flipConfig == curSPKdist-5);
                    CONdistIdx = find(CONdata(CONidx,2) == CONelec);
                    CON_distAT = CONdata(CONidx(CONdistIdx),3);
                    
                    % Average ATs
                    CON_AT = mean([CON_proxAT, CON_distAT]);
                    
                else % Take exact AT at this verticle level
                    CONelec = find(flipConfig == curSPKdist);
                    CONdistIdx = find(CONdata(CONidx,2) == CONelec);
                    CON_AT = CONdata(CONidx(CONdistIdx),3);
                end
                
                startATdiff = CON_AT - (SPK_AT - BDFstartTime); % Pos = SPK 1st
                
                if (startATdiff) <= threshold(7) && (startATdiff >= -0.5)
                    tally = tally + 1; % Add to tally
                    AT_temp = [AT_temp, startATdiff];
                end % if
                
            end % for
            
            % Check tally against threshold
            if tally >= threshold(9)
                % Save to master list
                ML(end+1,:) = [SPK,CON,mean(AT_temp),overlap,tally];
            end % if
            
        end % if
    end % for
end % for

% Sort masterlist so contractions are on the leftmost column
if ~isempty(ML)
    ML = sortrows(ML,2,'ascend');
    CON_list = ML(:,2);
    ML(:,2) = ML(:,1);
    ML(:,1) = CON_list;
end

end

