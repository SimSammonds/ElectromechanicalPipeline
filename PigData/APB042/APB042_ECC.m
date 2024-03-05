% EndoFLIP data analysis - Electrocontractile Coupling for APB042.
% This script will read-in and extract appropriate data points which were
% recorded simultaneously with electrical activity.

% Author: Sam Simmonds
% Date: 21st November 2022

%% Importing and processing data into usable format:

import = 1;
if import == 1
    filename1 = 'APB042_LIVE_DATA_1.TXT';
    data_t1 = readtable(filename1);
    data_m1 = readmatrix(filename1);
    [nRows1, nCols1] = size(data_m1);
    filename2 = 'APB042_LIVE_DATA_2.TXT';
    data_t2 = readtable(filename2);
    data_m2 = readmatrix(filename2);
    [nRows2, nCols2] = size(data_m2);
    filename3 = 'APB042_LIVE_DATA_3.TXT';
    data_t3 = readtable(filename3);
    data_m3 = readmatrix(filename3);
    [nRows3, nCols3] = size(data_m3);
    filename4 = 'APB042_LIVE_DATA_4.TXT';
    data_t4 = readtable(filename4);
    data_m4 = readmatrix(filename4);
    [nRows4, nCols4] = size(data_m4);
    CON_data_t = [data_t1;data_t2;data_t3;data_t4];
    CON_data_m = [data_m1;data_m2;data_m3;data_m4];
    disp('Finished Importing Contaction Data')
    Fs_flip = 10;
end


%% Begin with probing - experiment 3

% PROBING IRRELEVANT BECAUSE WE HAVE ACTUAL BDF TIMINGS IN THIS STUDY!

% Start time: 10:40:35 am - End time: 10:41:35 am
% Entry - Data 3: 3863 (start of 10:40:35) - 4471 (end of 10:41:35)

Fs_flip = 10; % Sampling frequency (Hz)
% Extract and save to matrix:
probe = data_m3(3863:4471,4:20); % 1-16: Dest, 17: Pressure
[nRowsP,nColsP] = size(probe);

% Plot pressure over time:
subplot(5,1,1);
t_probe = (0:nRowsP-1)*(1/Fs_flip); % Establish time vector
plot(t_probe,probe(:,17),'k');
%xlim([150,450]);

% Plot Dest over time across the entire catheter:
subplot(5,1,2);
plot(t_probe,probe(:,1:8));
%xlim([160,210]);

subplot(5,1,3);
plot(t_probe,probe(:,3:10));
%xlim([210,260]);

subplot(5,1,4);
plot(t_probe,probe(:,7:14));
%xlim([260,330]);

subplot(5,1,5);
plot(t_probe,probe(:,9:16));
%xlim([320,440]);

% Superimpose BDF time stamps onto plot and shift until they line up with
% pressure peaks
shifter = 0; % The amount to shift the points in time (s)
start_time = 0; % Below times are subtracted to find relative time
points_s = [12,24,36,50,83,94,108,120,134,166,179,191]; % (s) start points
points_e = [18,30,40,55,87,98,113,126,139,172,184,196]; % (s) end points
points_s = points_s - start_time + shifter;
points_e = points_e - start_time + shifter;

for j = 1:5
    subplot(5,1,j);
    for col = 1:length(points_s)
        xline(points_s(col),'r');
        xline(points_e(col),'b');
    end
end

% Probing is NOT in the same bdf file as the baseline recordings.
% How can we correlate the two temporally?
% vent pauses? might see a change in Dest activity?


%% Experiment 3 - Contractions

% Start time: 10:40:35 am - End time: 10:41:35 am
% Entry - Data 3: 3863 (start of 10:40:35) - 4471 (end of 10:41:35)

% Timing workspace:
% 10:37:00 == BDF 305 == 73721
% Therefore BDF 0 == 73721 - 3050 = 70671
% check ... yeep

% Extract contractile data:
expDATA = data_m2(3863:4471,4:20); % 1-16: Dest, 17: Pressure
[nRows,~] = size(expDATA);

% Find pylorus (most common minimum CSA):
for row = 1:nRows % Pylorus == minimum CSA on EndoFLIP
    minimum = min(expDATA(row,1:16)); % Save as index 18
    expDATA(row,19) = find(expDATA(row,1:16)==minimum,1);
    expDATA(row,18) = minimum;
end
pyl_flip = mode(expDATA(:,19)); % Find most common sensor with minCSA
pyl_flip = 10; % Manual override - 
% 10 is the 2nd most common MinCSA and is more expected than sensor 15

% USER SETTINGS:
wind = [0,60]; % seconds
viewElecs = 1:16;

% Create a series of subplots for visualisation:
t_exp = (0:nRows-1)*(1/Fs_flip); % Establish time vector
for plt = 1:length(viewElecs)
    subplot(length(viewElecs),1,plt)
    plot(t_exp, expDATA(:,viewElecs(plt)))
    xlim(wind);
    ylabel(num2str(viewElecs(plt)))
    set(gca,'Xtick',[]);
    set(gca,'Ytick',[]);
end
set(gcf,'color','w');

% Manually save contractile time stamps to contraction list excel file

%% Import contraction, SW, and spike data :::::::::::::::::::::::::::::::::
import2 = 1;
if import2 == 1
    addpath 'C:\Users\ssim749\Documents\PhD Backup March 2022\Data Analysis\Pig Data\APB040'
    toappSWdata = load('APB042_exp3_520-580_SW_elec3','toapp');
    toappSPKdata = load('APB042_exp3_520-580_SPK_elec3','toapp');
    excelCONdata = readmatrix('APB042_ContractionClusterList.xlsx');
    SW_AMPs = readmatrix('APB042_SW_AMPs.xlsx');
    % excelCONdata IDX: [Cluster, Electrode, StartTime, EndTime, StartAmp, EndAmp, Location]
    oldElecConfig = readmatrix('Custom_Config_3.txt'); % Import electrode config
    close all
    
    % Organise slow wave and spike data:
    SWdata = toappSWdata.toapp.TimeAmplCluster(1,1:end-1);
    [numSWs, SWdata] = DataOrganiser(SWdata,1);
    SPKdata = toappSPKdata.toapp.TimeAmplSpikeCluster(1:200,1);
    [numSPKs, SPKdata] = DataOrganiser(SPKdata,2);
    
    % Manually remove SW = 7
%     SWdata(346:end-1,:) = SWdata(347:end,:);
%     SWdata(end,:) = [];
    
    Fs_SPK = 512; % (Hz)
    
    % Change everything so that it is threshold dependent
    Thrs_SW_SPK_temp = 1.5; % seconds %1
    Thrs_SW_SPK_sptl = 10; % shared electrode area (mm^2) %15
    Thrs_SW_SPK_tally = 1; % valid electrodes %3
    Thrs_SW_CON_temp = 3; % seconds %2
    Thrs_SW_CON_sptl = 10; % electrode overlap (mm) %15
    Thrs_SW_CON_tally = 1; %3
    Thrs_CON_SPK_temp = 2; % seconds %2
    Thrs_CON_SPK_sptl = 10; % mm %15
    Thrs_CON_SPK_tally = 1; %3
    
    thresholds = [Thrs_SW_SPK_temp,...
        Thrs_SW_SPK_sptl,...
        Thrs_SW_SPK_tally,...
        Thrs_SW_CON_temp,...
        Thrs_SW_CON_sptl,...
        Thrs_SW_CON_tally,...
        Thrs_CON_SPK_temp,...
        Thrs_CON_SPK_sptl,...
        Thrs_CON_SPK_tally];
    
end

%% NEW CORRELATIONS

% GEMS parameters and outputs
BDFstartTime = 520; % (s)
gasFreq = 3.81;
intFreq = 17.41;

% Reconfigure electrode structures
pyl_elec = 13; % Manually determined using GEMS - num of elecs from top of array?
g_clsts = [1:4];
i_clsts = [11:13;15:17;20:22;24:26];
[qz_AVE,qz_SD,pyl_elec] = QuiescentZoneFunc(SWdata,oldElecConfig,g_clsts,i_clsts);
[elecConfig,flipConfig] = RecentreElecConfigs(oldElecConfig,pyl_elec,pyl_flip);
qz_metrics = [qz_AVE,qz_SD];
i_clsts = 10:27; % Changing the list to incorporate all intestinal waves

% Extract valid contractions & metrics
[CONdata, ~] = ValidateContractions(flipConfig,excelCONdata,expDATA,t_exp);
CONmetricsForHistology = CONmetricsForHisto(CONdata);

% Find organ specific SWs:
[gSWdata,iSWdata] = organSWs(SWdata,g_clsts);
% Find directional SPKs:
[isotSPKdata,longSPKdata,circSPKdata] = directionalSPKs(SPKdata,oldElecConfig,1);

% For each combination, find associations -> check for overlap!
% IDX: [SWclst, SPKclst, meanStartATdiff, overlap(mm)]
isot_g_SW_SPK_ML = SW_SPK_Correlations(gSWdata,isotSPKdata,thresholds);
long_g_SW_SPK_ML = SW_SPK_Correlations(gSWdata,longSPKdata,thresholds);
circ_g_SW_SPK_ML = SW_SPK_Correlations(gSWdata,circSPKdata,thresholds);

isot_i_SW_SPK_ML = SW_SPK_Correlations(iSWdata,isotSPKdata,thresholds);
long_i_SW_SPK_ML = SW_SPK_Correlations(iSWdata,longSPKdata,thresholds);
circ_i_SW_SPK_ML = SW_SPK_Correlations(iSWdata,circSPKdata,thresholds);

% IDX: [CONclst, SPKclst, meanStartATdiff(s), overlap(mm), tally]
isot_CON_SPK_ML = CON_SPK_Correlations(isotSPKdata,CONdata,elecConfig,oldElecConfig,flipConfig,thresholds,BDFstartTime);
long_CON_SPK_ML = CON_SPK_Correlations(longSPKdata,CONdata,elecConfig,oldElecConfig,flipConfig,thresholds,BDFstartTime);
circ_CON_SPK_ML = CON_SPK_Correlations(circSPKdata,CONdata,elecConfig,oldElecConfig,flipConfig,thresholds,BDFstartTime);

% Metrics:
% Gastric SW-SPK associations:
% Isotopic
if ~isempty(isot_g_SW_SPK_ML)
    isot_g_SW_Met = Metrics('SW',isot_g_SW_SPK_ML(:,1),gSWdata,SW_AMPs,toappSPKdata,BDFstartTime);
    isot_g_SPK_Met = Metrics('SPK',isot_g_SW_SPK_ML(:,2),isotSPKdata,SW_AMPs,toappSPKdata,BDFstartTime);
    isot_g_SW_SPK_ForMet = FormatMetrics2('SW','SPK',isot_g_SW_Met,isot_g_SPK_Met,isot_g_SW_SPK_ML);
else; isot_g_SW_SPK_ForMet = NaN(1,6);
end
% Longitudinal:
if ~isempty(long_g_SW_SPK_ML)
    long_g_SW_Met = Metrics('SW',long_g_SW_SPK_ML(:,1),gSWdata,SW_AMPs,toappSPKdata,BDFstartTime);
    long_g_SPK_Met = Metrics('SPK',long_g_SW_SPK_ML(:,2),longSPKdata,SW_AMPs,toappSPKdata,BDFstartTime);
    long_g_SW_SPK_ForMet = FormatMetrics2('SW','SPK',long_g_SW_Met,long_g_SPK_Met,long_g_SW_SPK_ML);
else; long_g_SW_SPK_ForMet = NaN(1,6);
end
% Circular:
if ~isempty(circ_g_SW_SPK_ML)
    circ_g_SW_Met = Metrics('SW',circ_g_SW_SPK_ML(:,1),gSWdata,SW_AMPs,toappSPKdata,BDFstartTime);
    circ_g_SPK_Met = Metrics('SPK',circ_g_SW_SPK_ML(:,2),circSPKdata,SW_AMPs,toappSPKdata,BDFstartTime);
    circ_g_SW_SPK_ForMet = FormatMetrics2('SW','SPK',circ_g_SW_Met,circ_g_SPK_Met,circ_g_SW_SPK_ML);
else; circ_g_SW_SPK_ForMet = NaN(1,6);
end

% Intestinal SW-SPK associations:
% Isotopic
if ~isempty(isot_i_SW_SPK_ML)
    isot_i_SW_Met = Metrics('SW',isot_i_SW_SPK_ML(:,1),iSWdata,SW_AMPs,toappSPKdata,BDFstartTime);
    isot_i_SPK_Met = Metrics('SPK',isot_i_SW_SPK_ML(:,2),isotSPKdata,SW_AMPs,toappSPKdata,BDFstartTime);
    isot_i_SW_SPK_ForMet = FormatMetrics2('SW','SPK',isot_i_SW_Met,isot_i_SPK_Met,isot_i_SW_SPK_ML);
else; isot_i_SW_SPK_ForMet = NaN(1,6);
end
% Longitudinal:
if ~isempty(long_i_SW_SPK_ML)
    long_i_SW_Met = Metrics('SW',long_i_SW_SPK_ML(:,1),iSWdata,SW_AMPs,toappSPKdata,BDFstartTime);
    long_i_SPK_Met = Metrics('SPK',long_i_SW_SPK_ML(:,2),longSPKdata,SW_AMPs,toappSPKdata,BDFstartTime);
    long_i_SW_SPK_ForMet = FormatMetrics2('SW','SPK',long_i_SW_Met,long_i_SPK_Met,long_i_SW_SPK_ML);
else; long_i_SW_SPK_ForMet = NaN(1,6);
end
% Circular:
if ~isempty(circ_i_SW_SPK_ML)
    circ_i_SW_Met = Metrics('SW',circ_i_SW_SPK_ML(:,1),iSWdata,SW_AMPs,toappSPKdata,BDFstartTime);
    circ_i_SPK_Met = Metrics('SPK',circ_i_SW_SPK_ML(:,2),circSPKdata,SW_AMPs,toappSPKdata,BDFstartTime);
    circ_i_SW_SPK_ForMet = FormatMetrics2('SW','SPK',circ_i_SW_Met,circ_i_SPK_Met,circ_i_SW_SPK_ML);
else; circ_i_SW_SPK_ForMet = NaN(1,6);
end

% CON - SPK associations
% Isotopic:
if ~isempty(isot_CON_SPK_ML)
    isot_CON_Met = Metrics('CON',isot_CON_SPK_ML(:,1),CONdata,SW_AMPs,toappSPKdata,BDFstartTime);
    isot_SPK_Met = Metrics('SPK',isot_CON_SPK_ML(:,2),isotSPKdata,SW_AMPs,toappSPKdata,BDFstartTime);
    isot_CON_SPK_ForMet = FormatMetrics2('CON','SPK',isot_CON_Met,isot_SPK_Met,isot_CON_SPK_ML);
else; isot_CON_SPK_ForMet = NaN(1,6);
end
% Longitudinal:
if ~isempty(long_CON_SPK_ML)
    long_CON_Met = Metrics('CON',long_CON_SPK_ML(:,1),CONdata,SW_AMPs,toappSPKdata,BDFstartTime);
    long_SPK_Met = Metrics('SPK',long_CON_SPK_ML(:,2),longSPKdata,SW_AMPs,toappSPKdata,BDFstartTime);
    long_CON_SPK_ForMet = FormatMetrics2('CON','SPK',long_CON_Met,long_SPK_Met,long_CON_SPK_ML);
else; long_CON_SPK_ForMet = NaN(1,6);
end
% Circular:
if ~isempty(circ_CON_SPK_ML)
    circ_CON_Met = Metrics('CON',circ_CON_SPK_ML(:,1),CONdata,SW_AMPs,toappSPKdata,BDFstartTime);
    circ_SPK_Met = Metrics('SPK',circ_CON_SPK_ML(:,2),circSPKdata,SW_AMPs,toappSPKdata,BDFstartTime);
    circ_CON_SPK_ForMet = FormatMetrics2('CON','SPK',circ_CON_Met,circ_SPK_Met,circ_CON_SPK_ML);
else; circ_CON_SPK_ForMet = NaN(1,6);
end

% Merge all formated metrics into one large variable for extraction:
All_ForMet = mergeMetrics2(isot_g_SW_SPK_ForMet,long_g_SW_SPK_ForMet,...
    circ_g_SW_SPK_ForMet,isot_i_SW_SPK_ForMet,long_i_SW_SPK_ForMet,...
    circ_i_SW_SPK_ForMet,long_CON_SPK_ForMet,circ_CON_SPK_ForMet,isot_CON_SPK_ForMet);

VeloDirectionMetrics = Propagation3(gSWdata,iSWdata,SPKdata,oldElecConfig);
PropMetrics = Propagation2(gSWdata,iSWdata,isotSPKdata,longSPKdata,circSPKdata,oldElecConfig);

disp('DONE!')

%% Slow Wave - Spike Correlations

% Master List
SW_SPK_ML = [];
% IDX: [SWclst, SPKclst, meanStartATdiff, overlap(mm)]

SW_SPK_ML = SW_SPK_Correlations(SWdata,SPKdata,thresholds);

% Slow wave - Contraction Correlations
BDFstartTime = 520; % (s)
gasFreq = 3.81;
intFreq = 17.41;

% Reconfigure electrode structures
pyl_elec = 13; % Manually determined using GEMS - num of elecs from top of array?
g_clsts = [1:4];
i_clsts = [11:13;15:17;20:22;24:26];
[qz_AVE,qz_SD,pyl_elec] = QuiescentZoneFunc(SWdata,oldElecConfig,g_clsts,i_clsts);
[elecConfig,flipConfig] = RecentreElecConfigs(oldElecConfig,pyl_elec,pyl_flip);
qz_metrics = [qz_AVE,qz_SD];
% Extract valid contractions & metrics
[CONdata, ~] = ValidateContractions(flipConfig,excelCONdata,expDATA,t_exp);

i_clsts = 10:27; % Changing the list to incorporate all intestinal waves

% Master List
SW_CON_ML = [];
% IDX: [SWclst, CONclst, meanStartATdiff, overlap(mm)]

SW_CON_ML = SW_CON_Correlations(SWdata,CONdata,elecConfig,oldElecConfig,flipConfig,thresholds,BDFstartTime);

% Contraction - Spike Correlations

% Master List
CON_SPK_ML = [];
% IDX: [CONclst, SPKclst, meanStartATdiff, overlap(mm)]

CON_SPK_ML = CON_SPK_Correlations(SPKdata,CONdata,elecConfig,oldElecConfig,flipConfig,thresholds,BDFstartTime);

% Slow waves correlated with spikes
SWs = unique(SW_SPK_ML(:,1));
metrics_temp = [];
SW_SPK_metrics = [];
for SWclst = SWs'
    metrics_temp = SlowWaveMetrics(SWclst,SWdata,oldElecConfig,elecConfig);
    SW_SPK_metrics = [SW_SPK_metrics; metrics_temp];
end

% Slow waves correlated with contractions
SWs = unique(SW_CON_ML(:,1));
metrics_temp = [];
SW_CON_metrics = [];
for SWclst = SWs'
    metrics_temp = SlowWaveMetrics(SWclst,SWdata,oldElecConfig,elecConfig);
    SW_CON_metrics = [SW_CON_metrics; metrics_temp];
end

% Contractions correlated with slow waves
CONs = unique(SW_CON_ML(:,2));
metrics_temp = [];
CON_SW_metrics = [];
for CONclst = CONs'
    metrics_temp = ContractionMetrics(CONclst,flipConfig,CONdata,expDATA,t_exp);
    CON_SW_metrics = [CON_SW_metrics; metrics_temp];
end

% Contractions correlated with spikes
CONs = unique(CON_SPK_ML(:,1));
metrics_temp = [];
CON_SPK_metrics = [];
for CONclst = CONs'
    metrics_temp = ContractionMetrics(CONclst,flipConfig,CONdata,expDATA,t_exp);
    CON_SPK_metrics = [CON_SPK_metrics; metrics_temp];
end

% Spikes associated with slow waves
SPKs = unique(SW_SPK_ML(:,2));
metrics_temp = [];
SPK_SW_metrics = [];
for SPKclst = SPKs'
    metrics_temp = SpikeMetrics(SPKdata,toappSPKdata,SPKclst,Fs_SPK,oldElecConfig,elecConfig,BDFstartTime);
    SPK_SW_metrics = [SPK_SW_metrics; metrics_temp];
end

% Spikes associated with contractions
SPKs = unique(CON_SPK_ML(:,2));
metrics_temp = [];
SPK_CON_metrics = [];
for SPKclst = SPKs'
    metrics_temp = SpikeMetrics(SPKdata,toappSPKdata,SPKclst,Fs_SPK,oldElecConfig,elecConfig,BDFstartTime);
    SPK_CON_metrics = [SPK_CON_metrics; metrics_temp];
end

% Numerical metrics for SWs and SPKs
[firstATstart,firstATend] = CONtiming(CONdata,flipConfig);
assoc_SWs = unique(SW_SPK_ML(:,1));
g_assoc = [];
i_assoc = [];
for i = 1:length(assoc_SWs)
    idx = find(SW_SPK_ML(:,1) == assoc_SWs(i));
    assoc_SPKs = unique(SW_SPK_ML(idx,2));
    if ismember(assoc_SWs(i),g_clsts)
        g_assoc = [g_assoc, numel(assoc_SPKs)];
    else
        i_assoc = [i_assoc, numel(assoc_SPKs)];
    end
end
disp(numel(unique(g_clsts)))
disp(numel(unique(i_clsts)))
disp(numel(unique(SPKdata(:,1))))
disp(mean(g_assoc))
disp(mean(i_assoc))

% Formatting metrics
[SW_SPK_FM,gSW_SPK_FM,iSW_SPK_FM] = FormatMetrics("SW","SPK",SW_SPK_metrics,SPK_SW_metrics,SW_SPK_ML,SW_AMPs,g_clsts);
[SW_CON_FM,gSW_CON_FM,iSW_CON_FM] = FormatMetrics("SW","CON",SW_CON_metrics,CON_SW_metrics,SW_CON_ML,SW_AMPs,g_clsts);
[CON_SPK_FM,~,~] = FormatMetrics("CON","SPK",CON_SPK_metrics,SPK_CON_metrics,CON_SPK_ML,SW_AMPs,g_clsts);  
[MergeData] = mergeMetrics(SW_SPK_FM,SW_CON_FM,CON_SPK_FM);
[gMergeData] = mergeMetrics(gSW_SPK_FM,gSW_CON_FM,[]);
[iMergeData] = mergeMetrics(iSW_SPK_FM,iSW_CON_FM,[]);

[SPKprop,totals] = SpikePropagation(SPKdata,oldElecConfig,elecConfig);
RgnlPWR = PowerPerLocation(toappSWdata,elecConfig,oldElecConfig,gasFreq,intFreq);
AmpGrad = AmplitudeGradient(SWdata,toappSWdata,SPKdata,toappSPKdata,CONdata,elecConfig,oldElecConfig,flipConfig,BDFstartTime);
[ap_list,ap_metrics,dp_list,dp_metrics] = Pacemaker(toappSWdata,elecConfig,oldElecConfig,g_clsts,i_clsts,qz_AVE);
[PropOut] = Propagation(SWdata,toappSWdata,SPKdata,elecConfig,oldElecConfig,g_clsts,i_clsts,qz_AVE);

% MasterList means and SDs:
SW_SPK_smry = [mean(SW_SPK_ML(:,3)),std(SW_SPK_ML(:,3)),mean(SW_SPK_ML(:,4)),std(SW_SPK_ML(:,4))];
SW_CON_smry = [mean(SW_CON_ML(:,3)),std(SW_CON_ML(:,3)),mean(SW_CON_ML(:,4)),std(SW_CON_ML(:,4))];
CON_SPK_smry = [mean(CON_SPK_ML(:,3)),std(CON_SPK_ML(:,3)),mean(CON_SPK_ML(:,4)),std(CON_SPK_ML(:,4))];
all_smry = [SW_SPK_smry,SW_CON_smry,CON_SPK_smry];

disp('DONE!');

%% Quantify values from each region for comparison

% Define region boundaries (electrical clusters)
gas_SWs = (1:8);
gas_SPKs = (1:50);
pyl_SPKs = (51:100);
duo_SWs = (10:100);
duo_SPKs = (101:200);

gas_SW_metrics = [];
for gas_SW = gas_SWs
    temp = [];
    if ismember(gas_SW,SWdata(:,1)) % Check it exist
        % Get metrics
        temp = SlowWaveMetrics(gas_SW,SWdata,oldElecConfig,elecConfig);
        % Save to permanent list
        gas_SW_metrics = [gas_SW_metrics; temp];
    end
end

duo_SW_metrics = [];
for duo_SW = duo_SWs
    temp = [];
    if ismember(duo_SW,SWdata(:,1)) % Check it exist
        % Get metrics
        temp = SlowWaveMetrics(duo_SW,SWdata,oldElecConfig,elecConfig);
        % Save to permanent list
        duo_SW_metrics = [duo_SW_metrics; temp];
    end
end

gas_SPK_metrics = [];
for gas_SPK = gas_SPKs
    temp = [];
    if ismember(gas_SPK,SPKdata(:,1)) % Check it exist
        % Get metrics
        temp = SpikeMetrics(SPKdata,toappSPKdata,gas_SPK,Fs_SPK,oldElecConfig,elecConfig,BDFstartTime);
        % Save to permanent list
        gas_SPK_metrics = [gas_SPK_metrics; temp];
    end
end

duo_SPK_metrics = [];
for duo_SPK = duo_SPKs
    temp = [];
    if ismember(duo_SPK,SPKdata(:,1)) % Check it exist
        % Get metrics
        temp = SpikeMetrics(SPKdata,toappSPKdata,duo_SPK,Fs_SPK,oldElecConfig,elecConfig,BDFstartTime);
        % Save to permanent list
        duo_SPK_metrics = [duo_SPK_metrics; temp];
    end
end

pyl_SPK_metrics = [];
for pyl_SPK = pyl_SPKs
    temp = [];
    if ismember(pyl_SPK,SPKdata(:,1)) % Check it exist
        % Get metrics
        temp = SpikeMetrics(SPKdata,toappSPKdata,pyl_SPK,Fs_SPK,oldElecConfig,elecConfig,BDFstartTime);
        % Save to permanent list
        pyl_SPK_metrics = [pyl_SPK_metrics; temp];
    end
end

disp('DONE!');

%% Spikes & SWs not associated with contractions

% list of spikes
all_SPKs = unique(SPKdata(:,1));

tonic_SPK_metrics = [];

% for each spike
for SPK = all_SPKs'
    
    % if not a associated with a contraction
    if ~ismember(SPK,CON_SPK_ML(:,2))
        
        % get spike metrics
        temp = SpikeMetrics(SPKdata,toappSPKdata,SPK,Fs_SPK,oldElecConfig,elecConfig);
        tonic_SPK_metrics = [tonic_SPK_metrics; temp];
    end

end


% list of spikes
all_SWs = unique(SWdata(:,1));

tonic_SW_metrics = [];

% for each spike
for SW = all_SWs'
    
    % if not a associated with a contraction
    if ~ismember(SW,SW_CON_ML(:,1))
        
        % get spike metrics
        temp = SlowWaveMetrics(SW,SWdata,oldElecConfig,elecConfig);
        % Amp
        SW_AMP_idx = find(SW_AMPs(:,1) == SW); % Slow wave amplitude data
        SW_AMP = SW_AMPs(SW_AMP_idx,2);
        temp = [temp, SW_AMP];
        tonic_SW_metrics = [tonic_SW_metrics; temp];
    end

end


%% Spike propagation per organ

% Find all spikes which touch pylorus and categorise their organ-level
% propagation

% List of quiescent electrodes
g_idx = find(SWdata(:,1) <= 8);
g_elecs = unique(SWdata(g_idx,2))';

i_idx = find(SWdata(:,1) >= 9);
i_elecs = unique(SWdata(i_idx,2))';

sw_elecs = [g_elecs, i_elecs];
p_elecs = setxor(1:256,sw_elecs);

p_elecs2 = [129 125 121 117 113 109 105 101 99 130 126 122 118 114 110 106 102 100 131 127 123 119 115 111 107 103 128 124 120 116 108 104];

% Initialise tallies
sd_tally = 0;
d_tally = 0;
s_tally = 0;
iso_tally = 0;
qz_tally = 0;

for SPK = unique(SPKdata(:,1))' % For each spike
    
    SPKidx = find(SPKdata(:,1) == SPK);
    elecs = SPKdata(SPKidx,2)';
    
    
    if (sum(ismember(elecs,p_elecs2)) > 0) % If it touches quiescent zone
        qz_tally = qz_tally + 1;
        
        if (sum(ismember(elecs,g_elecs))>0) && (sum(ismember(elecs,i_elecs))>0) % If touches both
            sd_tally = sd_tally + 1;
            
        elseif (sum(ismember(elecs,i_elecs))>0) % If touches duodenum
            d_tally = d_tally + 1;
            
        elseif (sum(ismember(elecs,g_elecs))>0) % If touches stomach
            s_tally = s_tally + 1;
            
        else % Else isolated to pylorus
            iso_tally = iso_tally + 1;
        end
        
    end
    
end

% Convert to percentages
sd_perc = sd_tally / qz_tally * 100;
d_perc = d_tally / qz_tally * 100;
s_perc = s_tally / qz_tally * 100;
iso_perc = iso_tally / qz_tally * 100;

sd_tally
d_tally
s_tally
iso_tally

disp('DONE!');




