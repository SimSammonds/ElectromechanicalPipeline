function [perWaveVals] = Metrics(type,clsts,data,SW_AMPs,toappSPKdata,BDFstartTime)
% METRICS is a function to extract valuable data from marked events, which
% can be SWs, SPKs, or CONs.
% This is a merged and streamlined version of SlowWaveMetrics,
% SpikeMetrics, ContractionMetrics, MergeMetrics, and FormatMetrics.

% Outputs:
% perWaveVals: -> IDX: [clst, DUR, AMP, NRG]

% Author: Sam Simmonds
% Date: 1st September 2023

perWaveVals = [];
tol = 0.01;

% Determine which event has been fed:
if type == "SW"

    % perWaveVals:
    for clst = clsts' % For each slow wave that has been input
        idx = find(data(:,1) == clst);  % Find indices
        % Find temporal data:
        ATs = data(idx,3);
        DUR = max(ATs) - min(ATs);
        % Find amplitude data:
        SW_AMP_idx = find(SW_AMPs(:,1) == clst); % Slow wave amplitude data
        AMP = SW_AMPs(SW_AMP_idx,2);
        NRG = AMP*DUR*numel(ATs)*5; % Amp(mV) * Dur(s) * Area(mm^2)
        perWaveVals = [perWaveVals; clst, DUR, AMP, NRG];
    end

elseif type == "SPK"

    data(:,3:4) = data(:,3:4) - BDFstartTime; % Shift data

    % perWaveVals:
    for clst = clsts' % For each spike that has been input
        idx = find(data(:,1) == clst);  % Find indices
        % Find temporal data:
        startATs = data(idx,3);
        endATs = data(idx,4);
        DUR = max(endATs) - min(startATs);

        % Initialise primary variables:
        filtdata = toappSPKdata.toapp.filtdata;
        [r,c] = size(filtdata); % Determine size
        t = (0:c-1)*(1/512); % Create time vector

        % Initialise blank variables
        AMP = [];
        start_temp = [];
        end_temp = [];

        for i = idx' % For each spike marking
            elec = data(i,2);    % Find electrode

            % For when < 256 channels are filtered:
            % Need to ensure that FC swapped arrays are being treated properly
            if (r ~= length(toappSPKdata.toapp.showchans)) && (elec >= 97) % Check if FC swapped has been passed
                filt_elec_idx = elec - 32; % Adjust
            else
                filt_elec_idx = elec; % Keep the same
            end

            % Find start and end times indices
            ST_idx = find(abs(t-data(i,3))<tol,1,'last');
            ET_idx = find(abs(t-data(i,4))<tol,1,'last');
            start_temp = [start_temp, t(ST_idx)];
            end_temp = [end_temp, t(ET_idx)];

            % Max amplitude (to be averaged)
            sig = filtdata(filt_elec_idx,ST_idx:ET_idx);
            amp_temp = max(sig) - min(sig);
            AMP = [AMP, amp_temp];
        end

        AMP = mean(AMP)/1000;
        NRG = AMP*DUR*numel(idx)*5;
        perWaveVals = [perWaveVals; clst, DUR, AMP, NRG];
    end

elseif type == "CON"
    % perWaveVals:
    for clst = clsts' % For each slow wave that has been input
        idx = find(data(:,1) == clst);  % Find indices
        elecs = data(idx,2);         % Elecs
        AMP = mean(data(idx,7)); % AMP
        DUR = max(data(idx,5)) - min(data(idx,3));
        NRG = AMP*DUR*numel(elecs)*10;
        perWaveVals = [perWaveVals; clst, DUR, AMP, NRG];
    end
end
end























