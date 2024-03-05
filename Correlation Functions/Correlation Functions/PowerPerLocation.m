function [OUT] = PowerPerLocation(toappdata,elecConfig,oldElecConfig,gasFreq,intFreq)
% AmplitudeGradient is a function to determine how SW and CON Powers
% change as a function of spatial location. Averaging the power spectra of
% each electrode (RAW data; at a particular distance) containing gastric or intestinal slow wave signals will
% remove phasing effects. Then, taking the frequency bins of interest, we
% can determine the average power at each distance from the pylorus.

% Inputs: SWdata,toappSWdata,CONdata,ALL_SW_AMPs,elecConfig,oldElecConfig
% Outputs: matrix containing 1D amplitude changes for each of the event

% Author: Sam Simmonds
% Date: 4th July 2023

%% Disable warnings

warning('off','all')

%% Electrical Waveforms

all_unfilt = toappdata.toapp.edata(1:256,:);
[~,L] = size(all_unfilt);
tol = 0.01;
delta = 0.05; % Freq spread on each side of the fundamental and harmonics

% PSD parameters - 30 Hz
SegLength = L/2;          % Cannot be greater than number of samples
SampleOverlap = L/4;      % Increase this number to improve/reduce accuracy?
DFTpoints = L*10;        % Optimise this so ~3 cpm falls directly on a DFT bin
Fs = 30;                    % Sampling frequency
t = (0:L-1)*(1/Fs);      % Time vector (for plotting)

% Index each electrode based upon it's distance from the pylorus
chan_list = toappdata.toapp.showchans;

for i = 1:length(chan_list) % For each channel
    
    curElec = chan_list(i);
    
    % Find distance
    [r,c] = find(oldElecConfig == curElec);
    DIST(curElec) = elecConfig(r,c);
    
    % Extract unfiltered data
    unfilt = all_unfilt(curElec,:);
    unfilt = unfilt./1000; % Change units from uV to mV
    
    % Perform PSD
    [PSD_temp, W_temp] = pwelch(unfilt,SegLength,SampleOverlap,DFTpoints,Fs);
    
    % Save with Distance label for averaging later
    PSD(curElec,:) = PSD_temp';
    
    if i == 1
        HZ = W_temp';
        CPM = 60*W_temp';    % Convert to cpm (Hz*60 = cpm)
    end
end

MAX = max(DIST);
MIN = min(DIST);

PSD2 = [];
DIST2 = [];

for d = MAX:-5:MIN % For each distance
    
    idx = find(DIST == d); % Find appropriate elecs
    AVE_PSD_temp = mean(PSD(idx,:)); % Average each elec at each time point
    
    % Save to lists
    PSD2 = [PSD2; AVE_PSD_temp];
    DIST2 = [DIST2; d];

end

% Find gastric bins
GF_IDX_1 = find(abs(CPM-(gasFreq-delta))<tol,1,'last');
GF_IDX_2 = find(abs(CPM-(gasFreq+delta))<tol,1,'last');
GH1_IDX_1 = find(abs(CPM-(gasFreq*2-delta))<tol,1,'last');
GH1_IDX_2 = find(abs(CPM-(gasFreq*2+delta))<tol,1,'last');
GH2_IDX_1 = find(abs(CPM-(gasFreq*3-delta))<tol,1,'last');
GH2_IDX_2 = find(abs(CPM-(gasFreq*3+delta))<tol,1,'last');

% Find intestinal bins
IF_IDX_1 = find(abs(CPM-(intFreq-delta))<tol,1,'last');
IF_IDX_2 = find(abs(CPM-(intFreq+delta))<tol,1,'last');
IH1_IDX_1 = find(abs(CPM-(intFreq*2-delta))<tol,1,'last');
IH1_IDX_2 = find(abs(CPM-(intFreq*2+delta))<tol,1,'last');
IH2_IDX_1 = find(abs(CPM-(intFreq*3-delta))<tol,1,'last');
IH2_IDX_2 = find(abs(CPM-(intFreq*3+delta))<tol,1,'last');

% Ventilator bins
vf = 22.8; % Ventilator frequency
VF_IDX_1 = find(abs(CPM-(vf-delta))<tol,1,'last');
VF_IDX_2 = find(abs(CPM-(vf+delta))<tol,1,'last');
VH1_IDX_1 = find(abs(CPM-(vf*2-delta))<tol,1,'last');
VH1_IDX_2 = find(abs(CPM-(vf*2+delta))<tol,1,'last');
VH2_IDX_1 = find(abs(CPM-(vf*3-delta))<tol,1,'last');
VH2_IDX_2 = find(abs(CPM-(vf*3+delta))<tol,1,'last');


for i = 1:length(DIST2) % For each location
    
    % Perform integration over defined limits
    PWR(i,1) = trapz(CPM(GF_IDX_1:GF_IDX_2),PSD2(i,GF_IDX_1:GF_IDX_2));
    PWR(i,2) = trapz(CPM(GH1_IDX_1:GH1_IDX_2),PSD2(i,GH1_IDX_1:GH1_IDX_2));
    PWR(i,3) = trapz(CPM(GH2_IDX_1:GH2_IDX_2),PSD2(i,GH2_IDX_1:GH2_IDX_2));
    
    PWR(i,4) = trapz(CPM(IF_IDX_1:IF_IDX_2),PSD2(i,IF_IDX_1:IF_IDX_2));
    PWR(i,5) = trapz(CPM(IH1_IDX_1:IH1_IDX_2),PSD2(i,IH1_IDX_1:IH1_IDX_2));
    PWR(i,6) = trapz(CPM(IH2_IDX_1:IH2_IDX_2),PSD2(i,IH2_IDX_1:IH2_IDX_2));
    
    PWR(i,7) = trapz(CPM(VF_IDX_1:VF_IDX_2),PSD2(i,VF_IDX_1:VF_IDX_2));
    PWR(i,8) = trapz(CPM(VH1_IDX_1:VH1_IDX_2),PSD2(i,VH1_IDX_1:VH1_IDX_2));
    PWR(i,9) = trapz(CPM(VH2_IDX_1:VH2_IDX_2),PSD2(i,VH2_IDX_1:VH2_IDX_2));
    
    PWR(i,10) = sum(PWR(i,1:3)); % Sum gastric powers
    PWR(i,11) = sum(PWR(i,4:6)); % Sum intestinal powers
    PWR(i,12) = sum(PWR(i,7:9)); % Sum ventilator powers
    
%     % Bandpower() - is there any difference between manual and auto?
%     BP(i,1) = bandpower(PSD2(i,:),CPM,[CPM(GF_IDX_1),CPM(GF_IDX_2)],"psd");
%     BP(i,2) = bandpower(PSD2(i,:),CPM,[CPM(GH1_IDX_1),CPM(GH1_IDX_2)],"psd");
%     BP(i,3) = bandpower(PSD2(i,:),CPM,[CPM(GH2_IDX_1),CPM(GH2_IDX_2)],"psd");
%     
%     BP(i,4) = bandpower(PSD2(i,:),CPM,[CPM(IF_IDX_1),CPM(IF_IDX_2)],"psd");
%     BP(i,5) = bandpower(PSD2(i,:),CPM,[CPM(IH1_IDX_1),CPM(IH1_IDX_2)],"psd");
%     BP(i,6) = bandpower(PSD2(i,:),CPM,[CPM(IH2_IDX_1),CPM(IH2_IDX_2)],"psd");
%     
%     BP(i,7) = sum(BP(i,1:3)); % Sum gastric powers
%     BP(i,8) = sum(BP(i,4:6)); % Sum intestinal powers
    
end

% Output gastric and intestinal powers
OUT(:,1) = DIST2;
OUT(:,2) = PWR(:,10);
OUT(:,3) = PWR(:,11);
OUT(:,4) = PWR(:,12);


end
























