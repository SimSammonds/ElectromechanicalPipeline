% EndoFLIP data analysis - Electrocontractile Coupling for APC005.
% This script will read-in and extract appropriate data points which were
% recorded simultaneously with electrical activity.

% Author: Sam Simmonds
% Date: 9th June 2023

%% Importing and processing data into usable format:

import = 1;
if import == 1
    Fs_flip = 10; % Sampling frequency (Hz)
    filename1 = 'APC005_LIVE_DATA_1.TXT';
    data_t1 = readtable(filename1);
    data_m1 = readmatrix(filename1);
    filename2 = 'APC005_LIVE_DATA_2.TXT';
    data_t2 = readtable(filename2);
    data_m2 = readmatrix(filename2);
    filename3 = 'APC005_LIVE_DATA_3.TXT';
    data_t3 = readtable(filename3);
    data_m3 = readmatrix(filename3);
    CON_data_t = [data_t1;data_t2;data_t3];
    CON_data_m = [data_m1;data_m2;data_m3];
    disp('Finished Importing Contaction Data')
end

%% Experiment 2 - Baseline Contractions

% Time corre workspace:
% 330 bdf == 10:41:18 == entry 25438
% bdf 0 == 22138

% Start time: 630s (BDF) / 28438 (flip entry)
% End time: 750s (BDF) / 29638 (flip entry)

CON_BL = CON_data_m(28438:29638,4:20); % 1-16: Dest, 17: Pressure
[r,c] = size(CON_BL);

% Find pylorus (most common minimum CSA):
for row = 1:r % Pylorus == minimum CSA on EndoFLIP
    minimum = min(CON_BL(row,1:16)); % Save as index 18
    CON_BL(row,19) = find(CON_BL(row,1:16)==minimum,1);
    CON_BL(row,18) = minimum;
end
pyl_flip = mode(CON_BL(:,19)); % Find most common sensor with minCSA

% USER SETTINGS:
wind = [0,120]; % seconds
viewElecs = 1:16;

% Create a series of subplots for visualisation:
t_CON_BL = (0:r-1)*(1/Fs_flip); % Establish time vector

% Pressure
% subplot(length(viewElecs)+1,1,1)
% plot(t_CON_BL,CON_BL(:,17),'k');
% xlim(wind);
% ylabel('Pressure')
% set(gca,'Xtick',[]);
% set(gca,'Ytick',[]);

% Dest
for plt = 1:length(viewElecs)
    subplot(length(viewElecs),1,plt)
    plot(t_CON_BL, CON_BL(:,viewElecs(plt)))
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
    addpath 'C:\Users\ssim749\Desktop\PhD - New\Data Analysis\Pig Data\APC007'
    toappSWdata = load('APC005_exp3_630_750_SW.mat','toapp');
    toappSPKdata = load('APC005_exp3_630_750_SPK.mat','toapp');
    excelCONdata = readmatrix('APC005_ContractionClusterList.xlsx');
    SW_AMPs = readmatrix('APC005_SW_AMPs.xlsx');
    %SW_AMPs = SW_AMPs(49:end,:);
    % excelCONdata IDX: [Cluster, Electrode, StartTime, EndTime, StartAmp, EndAmp, Location]
    oldElecConfig = readmatrix('Config_Swapped_FC_2.txt'); % Import electrode config
    close all
    
    % Organise slow wave and spike data:
    SWdata = toappSWdata.toapp.TimeAmplCluster(1,1:end-1);
    [numSWs, SWdata] = DataOrganiser(SWdata,1);
    SPKdata = toappSPKdata.toapp.TimeAmplSpikeCluster(1:200,1);
    [numSPKs, SPKdata] = DataOrganiser(SPKdata,2);
    
    % Swap start and end time of SW X in electrode Y
    
    Fs_SPK = 512; % (Hz)
    Fs_flip = 10; % (Hz)
    
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
    disp('Finished importing SW and SPK data')
    
end

%% NEW CORRELATIONS

% GEMS parameters and outputs
BDFstartTime = 630; % (s)
gasFreq = 2.49;
intFreq = 17.44;

% Reconfigure electrode structures
pyl_elec = 14; % Determined in the QuiescentZone script
g_clsts = [40:44];
i_clsts = [2:4;9:11;16:18;24:26;31:33];
[qz_AVE,qz_SD,pyl_elec] = QuiescentZoneFunc(SWdata,oldElecConfig,g_clsts,i_clsts);
qz_metrics = [qz_AVE,qz_SD];
[elecConfig,flipConfig] = RecentreElecConfigs(oldElecConfig,pyl_elec,pyl_flip);
i_clsts = 2:35;

% Extract valid contractions & metrics
[CONdata, ~] = ValidateContractions(flipConfig,excelCONdata,CON_BL,t_CON_BL);

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

%% Correlations

% Master List
SW_SPK_ML = [];
% IDX: [SWclst, SPKclst, meanStartATdiff, overlap(mm)]

SW_SPK_ML = SW_SPK_Correlations(SWdata,SPKdata,thresholds);

% Slow wave - Contraction Correlations
BDFstartTime = 630; % (s)
gasFreq = 2.49;
intFreq = 17.44;

% Reconfigure electrode structures
pyl_elec = 14; % Determined in the QuiescentZone script
g_clsts = [40:44];
i_clsts = [2:4;9:11;16:18;24:26;31:33];
[qz_AVE,qz_SD,pyl_elec] = QuiescentZoneFunc(SWdata,oldElecConfig,g_clsts,i_clsts);
qz_metrics = [qz_AVE,qz_SD];
[elecConfig,flipConfig] = RecentreElecConfigs(oldElecConfig,pyl_elec,pyl_flip);

i_clsts = 2:35;

% Extract valid contractions & metrics
[CONdata, ~] = ValidateContractions(flipConfig,excelCONdata,CON_BL,t_CON_BL);
CONmetricsForHistology = CONmetricsForHisto(CONdata);

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
    metrics_temp = ContractionMetrics(CONclst,flipConfig,CONdata,CON_BL,t_CON_BL);
    CON_SW_metrics = [CON_SW_metrics; metrics_temp];
end

% Contractions correlated with spikes
CONs = unique(CON_SPK_ML(:,1));
metrics_temp = [];
CON_SPK_metrics = [];
for CONclst = CONs'
    metrics_temp = ContractionMetrics(CONclst,flipConfig,CONdata,CON_BL,t_CON_BL);
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


%% Experiment 4 - Prokinetic Contractions COPIED FROM APC007 SCRIPT NEEDS UPDATING!

% Manaully identify contractions for: (extend window as needed)
% [X] exp3: 10:50:10 - 10:51:30 // 50 - 130 BDF // 37837 - 38693 FLIP
% [X] exp3: 10:55:40 - 10:57:00 // 350 - 430 BDF // 41125 - 41981 FLIP
% [X] exp4: 11:20:05 - 11:21:30 // 5 - 90 BDF // 55969 - 56819 FLIP
% [X] exp4: 11:24:55 - 11:27:15 // 295 - 375 BDF // 58609 - 60012 FLIP
% [X] exp5: 11:50:07 - 11:51:27 // 0 - 80 BDF // 73669 - 74472 FLIP
% [X] exp5: 11:55:07 - 11:56:27 // 300 - 380 BDF // 76655 - 77761 FLIP
% [X] exp6: 12:21:40 - 12:23:00 // 6 - 86 BDF // 92525 - 93631 FLIP
% [X] exp6: 12:26:40 - 12:28:00 // 306 - 386 BDF // 95514 - 96320 FLIP

% Temporal Correlation Workspace:
% 

% FLIP start times (These == 0 on marked BDF)
% idx == [end wave, BDF zero; 'next row == next SW file']
FLIP_tStart = [6,37837; 12,41125; 19,55969; 25,58609; 30,73669; 36,76655; 42,92525; 48,95514];

% Extract and save to matrix:
expP = CON_data_m(95514:96320,4:20); % 1-16: Dest, 17: Pressure
[nRows,~] = size(expP);

% USER SETTINGS:
wind = [0,90]; % seconds
viewElecs = 1:16;

% Create a series of subplots for visualisation:
t_CON_BL = (0:nRows-1)*(1/Fs_flip); % Establish time vector
subplot(length(viewElecs)+1,1,1)
plot(t_CON_BL,expP(:,17),'k');
xlim(wind);
ylabel('Pressure')
set(gca,'Xtick',[]);
set(gca,'Ytick',[]);

for plt = 1:length(viewElecs)
    subplot(length(viewElecs)+1,1,plt+1)
    plot(t_CON_BL, expP(:,viewElecs(plt)))
    xlim(wind);
    ylabel(num2str(viewElecs(plt)))
    set(gca,'Xtick',[]);
    set(gca,'Ytick',[]);
end
set(gcf,'color','w');

% Manually save contractile time stamps to contraction list excel file

%% Prokinetic Contraction Processing
% Determine Mean +/- SD for:
    % Amplitude, duration, propagation

CCL_P = readmatrix('APC007_CCL_P.xlsx');
ProConMet = ProContProc(CCL_P,CON_data_m,FLIP_tStart);
disp('Prokinetic contractions successfully processed!')

%% Prokinetic Spike Processing
% Determine Mean +/- SD for:
    % Amplitude, duration, propagation

% USER SETTINGS:
firstClst = 1; % clst
lastClst = 6; % clst
BDFstartTime = 306;  % (s)
toappSPKdata = load('APC007_exp5_prok3_2_SPK','toapp');

[elecs] = PassProkElecs(toappSPKdata,2,firstClst,lastClst);
SPKdata = toappSPKdata.toapp.TimeAmplSpikeCluster(firstClst:lastClst,1);
[numSPKs, SPKdata] = DataOrganiser(SPKdata,2);
Fs_SPK = 512; % (Hz)

% Im pretty sure that these numbers don't matter for SPK metrics:
pyl_elec = 14;
pyl_flip = 10;
oldElecConfig = readmatrix('Custom_Config_3.txt'); % Import electrode config
[elecConfig,~] = RecentreElecConfigs(oldElecConfig,pyl_elec,pyl_flip);
close all

% Spikes associated with slow waves
SPKs = unique(SPKdata(:,1));
metrics_temp = [];
ProSpkMet = [];
for SPKclst = 1:numSPKs
    metrics_temp = SpikeMetrics(SPKdata,toappSPKdata,SPKclst,Fs_SPK,oldElecConfig,elecConfig,BDFstartTime);
    ProSpkMet = [ProSpkMet; metrics_temp];
end

disp('Prokinetic spikes successfully processed!')


