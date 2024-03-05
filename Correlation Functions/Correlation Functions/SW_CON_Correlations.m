function [ML] = SW_CON_Correlations(SWdata,CONdata,elecConfig,oldElecConfig,flipConfig,threshold,BDFstartTime)
% SW_CON_Correlations is a function to find valid pairs of slow waves and
% contractions. 

% Outputs: ML (master list)
% IDX: [SWclst, CONclst, meanStartATdiff(s), overlap(mm), tally]

% Author: Sam Simmonds
% Date: 14th November 2022

ML = [];
SWclsts = unique(SWdata(:,1));
CONclsts = unique(CONdata(:,1));

for SW = SWclsts' % For each slow wave

    % Find idx, elec, AT data
    SWidx = find(SWdata(:,1) == SW);
    SWelecs = SWdata(SWidx,2);
    
    % SW distances
    for j = 1:length(SWelecs) % For each electrode
        [elec_r, elec_c] = find(oldElecConfig == SWelecs(j)); % Find co-ords
        SWdist(j) = elecConfig(elec_r,elec_c); % Get distance
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
        min_list = [min(SWdist), min(CONdist)];
        max_list = [max(SWdist), max(CONdist)];
        overlap = (min(max_list) - max(min_list));
        
        if overlap >= threshold(5) % Check overlap threshold
            
            sharedElec = []; % IDX for overlapping SW elecs
            % Which electrodes overlap?
            for i = 1:length(SWelecs) % For each electrode
                if (SWdist(i) <= max(CONdist)) && (SWdist(i) >= min(CONdist))
                    sharedElec = [sharedElec, i];
                end % if
            end % for
            
            for elec = SWelecs(sharedElec)' % for each overlapping slow wave elec
                
                % Find AT difference at this electrode
                curSWidx = find(SWdata(SWidx,2) == elec);
                SW_AT = SWdata(SWidx(curSWidx),3);
                [elec_r, elec_c] = find(oldElecConfig == elec); % Find co-ords
                curSWdist = elecConfig(elec_r,elec_c); % Get distance
                remainder = rem(curSWdist,10); % Check remainder
                
                if remainder ~= 0 % Then we must take an average of adjacent CON_ATs
                    % Proximal
                    CONelec = find(flipConfig == curSWdist+5);
                    CONdistIdx = find(CONdata(CONidx,2) == CONelec);
                    CON_proxAT = CONdata(CONidx(CONdistIdx),3);
                    
                    % Distal
                    CONelec = find(flipConfig == curSWdist-5);
                    CONdistIdx = find(CONdata(CONidx,2) == CONelec);
                    CON_distAT = CONdata(CONidx(CONdistIdx),3);
                    
                    % Average ATs
                    CON_AT = mean([CON_proxAT, CON_distAT]);
                    
                else % Take exact AT at this verticle level
                    CONelec = find(flipConfig == curSWdist);
                    CONdistIdx = find(CONdata(CONidx,2) == CONelec);
                    CON_AT = CONdata(CONidx(CONdistIdx),3);
                end

                startATdiff = CON_AT - (SW_AT - BDFstartTime); % Pos = SW 1st

                if (startATdiff <= threshold(4)) && (startATdiff >= -0.5)
                    tally = tally + 1; % Add to tally
                    AT_temp = [AT_temp, startATdiff];
                end % if
                
            end % for
            
            % Check tally against threshold
            if tally >= threshold(6)
                % Save to master list
                ML(end+1,:) = [SW,CON,mean(AT_temp),overlap,tally];
            end % if
        end % if
    end % for
end % for


end % function