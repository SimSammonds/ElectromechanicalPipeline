function [SPKmetrics] = ProSpkProc(inputArg1,inputArg2)
%ProSpkProc: Prokinetic Spike Processing
%   This function will take in spike data and output useful parameters
%   including:
    % Amplitude, duration, propagation.
% This script is essentially a mirror of SpikeMetrics.m, but will need to
% work for reduced channel numbers.

% SpikeMetrics is a function to determine average amplitude and spike
% signal energy from original toapp data, for a single spike cluster.
% SPKmetrics IDX: 1=SPKclst, 2=aveAmp(mV), 3=totalEnergy(mV^2), 4=duration(s), 
% 5=vrtproplength(mm), 6=startAT(s), 7=endAT(s)

% Author: Sam Simmonds
% Date: 31st May 2023

% Pull filtered data from toapp
filtdata = toappSPKdata.toapp.filtdata;
% 256 rows = channels;  % > 60,000 cols = time points (512 Hz)

% Shift spike data
SPKdata(:,3:4) = SPKdata(:,3:4) - BDFstartTime;

% Determine size
[~,c] = size(filtdata);

% Create time vector
t = (0:c-1)*(1/Fs);
% As we are working with shifted SPK data, need to shift time vector too?

% Find SPK indices
idx = find(SPKdata(:,1) == SPKclst);

% Initialise variables
tol = 0.01;
amp_temp = [];
energy = 0;
start_temp = [];
end_temp = [];

for i = idx' % For each spike marking
    
    % Find electrode
    elec = SPKdata(i,2);
    
    % SPECIAL CONDITION FOR SWAPPED F & C ELEC CONFIGS:
    F_C_Swapped = 1;
    if F_C_Swapped == 1 && elec > 96
            elec = elec - 32;
    end
    
    % Find start and end times indices
    startTime = find(abs(t-SPKdata(i,3))<tol,1,'last');
    endTime = find(abs(t-SPKdata(i,4))<tol,1,'last');
    start_temp = [start_temp, t(startTime)];
    end_temp = [end_temp, t(endTime)];
    
    % Segment signal and time vectors
    sig = filtdata(elec,startTime:endTime);
    sig2 = sig.^2;
    time = t(startTime:endTime);
    
    % Max amplitude (to be averaged)
    amp = max(sig) - min(sig);
    amp_temp = [amp_temp, amp];
    
    % Perform energy integration
    energy = trapz(time,sig2) + energy;
    
end

startTimeOut = min(start_temp);
endTimeOut = max(end_temp);
dur = (max(end_temp) - min(start_temp));

SPK_elecs = SPKdata(idx,2);
for i = 1:length(SPK_elecs)
    [elec_r, elec_c] = find(oldElecConfig == SPK_elecs(i));
    SPK_dist(i) = elecConfig(elec_r,elec_c);
end
propLen = abs(max(SPK_dist)-min(SPK_dist));

ave_amp = mean(amp_temp);
SPKmetrics = [SPKclst,ave_amp,energy,dur,propLen,startTimeOut,endTimeOut];

end
