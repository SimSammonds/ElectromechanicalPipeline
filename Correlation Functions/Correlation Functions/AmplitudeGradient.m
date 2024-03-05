function [OUT] = AmplitudeGradient(SWdata,toappSWdata,SPKdata,toappSPKdata,CONdata,elecConfig,oldElecConfig,flipConfig,BDFstart)
% AmplitudeGradient is a function to determine how SW, SPK, and CON AMPs
% change as a function of spatial location. SWs, SPKs, and CONs should be
% averaged per electrode, over the entire analysed period.

% Inputs: SWdata,toappSWdata,SPKdata,toappSPKdata,CONdata,ALL_SW_AMPs,elecConfig,oldElecConfig
% Outputs: matrix containing 1D amplitude changes for each of the event

% Author: Sam Simmonds
% Date: 20th June 2023

%% SWs
SWlabels = unique(SWdata(:,1));
numSWs = length(SWlabels);
SW_AMPL = [];

for i = 1:numSWs % For each slow wave that was clustered
   WaveData = toappSWdata.toapp.TimeAmplCluster{1,SWlabels(i)};
   
   for j = 1:height(WaveData) % For each electrode
       
       curElec = WaveData(j,5); % Find current elec
       [idxR,idxC] = find(oldElecConfig == curElec); % Find distance
       dist = elecConfig(idxR,idxC);
       
       amp = WaveData(j,7); % Amplitude pulled from toapp data
       
       if ~isnan(amp)
           temp_SW_AMPL = [curElec, dist, amp];
           SW_AMPL = [SW_AMPL; temp_SW_AMPL]; % Save to master list
       end
   end
end

SW_OUT = [];

% Average SWs at each distance and save for output
for d = max(max(elecConfig)):-5:min(min(elecConfig))
    idx = find(SW_AMPL(:,2) == d);
    temp_SW_OUT = [d,mean(SW_AMPL(idx,3))];
    SW_OUT = [SW_OUT; temp_SW_OUT];
end


%% SPKs
SPKlabels = unique(SPKdata(:,1));
numSPKs = length(SPKlabels);
SPK_AMPL = [];
tol = 0.01;

% Filtered Data
filtdata = toappSPKdata.toapp.filtdata;
[r,c] = size(filtdata);
Fs = 512; % (Hz)
t = (0:c-1)*(1/Fs); % Create time vector

for i = 1:numSPKs % For each spike that was clustered
   WaveData = toappSPKdata.toapp.TimeAmplSpikeCluster{SPKlabels(i),1};
   
   for j = 1:height(WaveData) % For each electrode
       
       curElec = WaveData(j,5); % Find current elec
       [idxR,idxC] = find(oldElecConfig == curElec); % Find distance
       dist = elecConfig(idxR,idxC);
       
       % For when < 256 channels are filtered:
       % Need to ensure that FC swapped arrays are being treated properly
       if (r ~= length(toappSPKdata.toapp.showchans)) && (curElec >= 97) % Check if FC swapped has been passed
           filt_elec_idx = curElec - 32; % Adjust
       else
           filt_elec_idx = curElec; % Keep the same
       end
       
       % Amplitude will need to be manually determined:
       startT = WaveData(j,3) - BDFstart;
       endT = WaveData(j,8) - BDFstart;
       
       startTime = find(abs(t-startT)<tol,1,'last');
       endTime = find(abs(t-endT)<tol,1,'last');
       sig = filtdata(filt_elec_idx,startTime:endTime);   % Relevant time points
       amp = max(sig) - min(sig);
       
       if ~isnan(amp)
           temp_SPK_AMPL = [curElec, dist, amp];
           SPK_AMPL = [SPK_AMPL; temp_SPK_AMPL]; % Save to master list
       end
   end
end

SPK_OUT = [];

% Average SPKs at each distance and save for output
for d = max(max(elecConfig)):-5:min(min(elecConfig))
    idx = find(SPK_AMPL(:,2) == d);
    temp_SPK_OUT = [d,mean(SPK_AMPL(idx,3))];
    SPK_OUT = [SPK_OUT; temp_SPK_OUT];
end


%% CONs

% [clsts, elec, startTime, nan, endTime, nan, AMPL, dur, location]
% [1,      2,      3,             5,            7,               ]

CON_AMPL = [];
CON_OUT = [];
[r,~] = size(CONdata);

for i = 1:r % Each row in CONdata
    elec = CONdata(i,2);
    dist = flipConfig(elec);
    amp = CONdata(i,7);
    
    temp_CON_AMPL = [elec, dist, amp];
    CON_AMPL = [CON_AMPL; temp_CON_AMPL];
end

% Average CONs at each distance and save for output
for d = max(max(flipConfig)):-5:min(min(flipConfig))
    idx = find(CON_AMPL(:,2) == d);
    temp_CON_OUT = [d,mean(CON_AMPL(idx,3))];
    CON_OUT = [CON_OUT; temp_CON_OUT];
end

%% Combine
OUT = [];

maximum1 = max(SW_OUT(:,1));
maximum2 = max(CON_OUT(:,1));
minimum1 = min(SW_OUT(:,1));
minimum2 = min(CON_OUT(:,1));

for d = max(maximum1,maximum2):-5:min(minimum1,minimum2)
    % SWs
    idx = find(SW_AMPL(:,2) == d);
    temp_SW_OUT = [mean(SW_AMPL(idx,3))];
    % SPKs
    idx = find(SPK_AMPL(:,2) == d);
    temp_SPK_OUT = [mean(SPK_AMPL(idx,3))];
    % CONs
    CON_idx = find(CON_AMPL(:,2) == d);
    temp_CON_OUT = [mean(CON_AMPL(CON_idx,3))];
    % Combined
    OUT = [OUT; d, temp_SW_OUT/1000, temp_SPK_OUT/1000, temp_CON_OUT];
end




end
























