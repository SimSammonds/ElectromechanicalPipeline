function [OrgData,gOrgData,iOrgData] = FormatMetrics(Type1,Type2,metrics1,metrics2,ML,SW_AMPs,g_clsts)
% FormatMetrics is a function to take correlations and metrics, and to
% organise them in such a way that they can be easily exported and
% immediately plotted in excel - save time reorganising data later.
% Type1 should be either 'SW' or 'CON', and Type2 should be 'SPK' or 'CON'

% Outputs: OrgData
% OrgData: [Type1clst, Type2clst, meanStartATdiff(s), overlap(mm), tally]

% Author: Sam Simmonds
% Date: 15th November 2022

OrgData = [];
gOrgData = [];
iOrgData = [];
TempData = [];

if Type1 == "SW" % For primary data type SWs
    
    SWclsts = unique(ML(:,1)); % Get slow wave clusters
    [r,~] = size(SWclsts);
    
    if Type2 == "SPK" % For SW-SPK correlations
        for i = 1:r % For each slow wave

            SW_AMP_idx = find(SW_AMPs(:,1) == SWclsts(i)); % Slow wave amplitude data
            SW_AMP = SW_AMPs(SW_AMP_idx,2);
            
            MLidx = find(ML(:,1) == SWclsts(i)); % Indexing for ML
            SW_Metidx = find(metrics1(:,1) == SWclsts(i)); % Indexing for Metrics
            SPKclsts = (ML(MLidx,2)); % Get SPK clusters for this SW
            [r2,~] = size(SPKclsts);
            
            % First 4 columns dedicated to SWs (4th is ampl - manual)
            TempData(1,1:3) = metrics1(SW_Metidx,:);
            TempData(1,4) = SW_AMP;
            
            for j = 1:r2 % For each spike correlate with this SW
                SPK_Metidx = find(metrics2(:,1) == SPKclsts(j)); % Indexing for Metrics
                % Columns 5-12 dedicated to SPKs
                TempData(j,5:11) = metrics2(SPK_Metidx,:);
            end
            
            % Averaged / Summed metrics
            dur = max(TempData(:,11)) - min(TempData(:,10));
            TempData(:,10:11) = 0;
            meanAmp = mean(TempData(:,6));
            totalEnergy = sum(TempData(:,7));
            TempData(1,10:12) = [dur, meanAmp, totalEnergy];
            OrgData = [OrgData; TempData]; % Save to output
            
            % Organ specific outputs
            if ismember(SWclsts(i),g_clsts)
                gOrgData = [gOrgData; TempData];
            else
                iOrgData = [iOrgData; TempData];
            end
            TempData = []; % Reset temp array
        end
        
    elseif Type2 == "CON" % For SW-CON correlations
        for i = 1:r % For each slow wave
            
            SW_AMP_idx = find(SW_AMPs(:,1) == SWclsts(i)); % Slow wave amplitude data
            SW_AMP = SW_AMPs(SW_AMP_idx,2);
            
            MLidx = find(ML(:,1) == SWclsts(i)); % Indexing for ML
            SW_Metidx = find(metrics1(:,1) == SWclsts(i)); % Indexing for Metrics
            CONclsts = (ML(MLidx,2)); % Get CON clusters for this SW
            [r2,~] = size(CONclsts);
            
            % First 4 columns dedicated to SWs (4th is ampl - manual)
            TempData(1,1:3) = metrics1(SW_Metidx,:);
            TempData(1,4) = SW_AMP;
            
            for j = 1:r2 % For each spike correlate with this SW
                CON_Metidx = find(metrics2(:,1) == CONclsts(j)); % Indexing for Metrics
                % Columns 5-12 dedicated to CONs
                TempData(j,5:11) = metrics2(CON_Metidx,:);
            end
            
            % Averaged / Summed metrics
            dur = max(TempData(:,11)) - min(TempData(:,10));
            TempData(:,10:11) = 0;
            sumAmp = sum(TempData(:,7));
            totalEnergy = sum(TempData(:,8));
            TempData(1,10:12) = [dur, sumAmp, totalEnergy];
            OrgData = [OrgData; TempData]; % Save to output
            
            % Organ specific outputs
            if ismember(SWclsts(i),g_clsts)
                gOrgData = [gOrgData; TempData];
            else
                iOrgData = [iOrgData; TempData];
            end
            TempData = []; % Reset temp array
            
        end
    else
        % Error
        disp('Non-valid type inputs for FormatMetrics function')
    end
end

if Type1 == 'CON' % For primary data type CONs
    
    CONclsts = unique(ML(:,1)); % Get contraction clusters
    [r,~] = size(CONclsts);
    
    if Type2 == 'SPK' % For CON - SPK correlations
        for i = 1:r % For each contraction
            
            MLidx = find(ML(:,1) == CONclsts(i)); % Indexing for ML
            CON_Metidx = find(metrics1(:,1) == CONclsts(i)); % Indexing for Metrics
            SPKclsts = (ML(MLidx,2)); % Get SPK clusters for this CON
            [r2,~] = size(SPKclsts);
            
            % First 6 columns dedicated to CONs
            TempData(1,1:6) = metrics1(CON_Metidx,1:6);
            
            for j = 1:r2 % For each spike correlate with this SW
                SPK_Metidx = find(metrics2(:,1) == SPKclsts(j)); % Indexing for Metrics
                % Columns 7-14 dedicated to SPKs
                TempData(j,7:13) = metrics2(SPK_Metidx,:);
            end
            
            % Averaged / Summed metrics
            dur = max(TempData(:,13)) - min(TempData(:,12));
            TempData(:,12:13) = 0;
            meanAmp = mean(TempData(:,8));
            totalEnergy = sum(TempData(:,9));
            TempData(1,12:14) = [dur, meanAmp, totalEnergy];
            OrgData = [OrgData; TempData]; % Save to output
            TempData = []; % Reset temp array
            
        end
        
    else
        % Error
        disp('Non-valid type inputs for FormatMetrics function')
    end
end

end

