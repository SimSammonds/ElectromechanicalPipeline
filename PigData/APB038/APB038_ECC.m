%% APB038 - Correlations Script 2 %%
% Rewriting this script in order to reassess how I am making my
% correlations

% Change everything so that it is threshold dependent
Thrs_SW_SPK_temp = 1.5; % seconds %1 [%15]
Thrs_SW_SPK_sptl = 10; % shared electrode area (mm^2) %15 [%10]
Thrs_SW_SPK_tally = 1; % valid electrodes %3
Thrs_SW_CON_temp = 3; % seconds %2 [%3]
Thrs_SW_CON_sptl = 10; % electrode overlap (mm) %15 [%10]
Thrs_SW_CON_tally = 1; %3
Thrs_CON_SPK_temp = 2; % seconds %2
Thrs_CON_SPK_sptl = 10; % mm %15 [%10]
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

% Thresholds tightened 1st September 2023, to attempt to segregate
% SW-dependent and SW-independent SPKs. values prior are in square
% parentheses.

%% Import contraction, SW, and spike data :::::::::::::::::::::::::::::::::

%cd 'C:\Users\ssim749\Desktop\PhD - New\Data Analysis\Pig Data\APB038'
%addpath 'C:\Users\ssim749\Desktop\PhD - New\Data Analysis\MATLAB Scripts\Correlation Functions'
%addpath 'C:\Users\ssim749\Desktop\PhD - New\Data Analysis\Pig Data\APB038'
toappSWdata = load('APB038_exp4_280-400_SWs_elec3','toapp');
toappSPKdata = load('APB038_exp4_280-400_SPK_elec3','toapp');
excelCONdata = readmatrix('APB038_ContractionClusterList.xlsx');
SW_AMPs = readmatrix('APB038_SW_AMPs.xlsx');
% excelCONdata IDX: [Cluster, Electrode, StartTime, EndTime, StartAmp, EndAmp, Location]
oldElecConfig = readmatrix('Custom_Config_3.txt'); % Import electrode config
close all

% Organise slow wave and spike data:
SWdata = toappSWdata.toapp.TimeAmplCluster(1,1:end-1);
[numSWs, SWdata] = DataOrganiser(SWdata,1);
SPKdata = toappSPKdata.toapp.TimeAmplSpikeCluster(1:200,1);
[numSPKs, SPKdata] = DataOrganiser(SPKdata,2);

% Manually remove outlier from APB038, SPKclst == 152
[r,c] = size(SPKdata);
if (r == 1527) && (c == 4)% If APB038
    SPKdata(1477:end-1,:) = SPKdata(1478:end,:);
    SPKdata(end,:) = [];
end

Fs_SPK = 512; % (Hz)

filename1 = 'APB038_LIVE_DATA_1.TXT';
data_m1 = readmatrix(filename1);
filename2 = 'APB038_LIVE_DATA_2.TXT';
data_m2 = readmatrix(filename2);
filename3 = 'APB038_LIVE_DATA_3.TXT';
data_m3 = readmatrix(filename3);

% Extract contractile data:
Fs_flip = 10; % (Hz)
exp4 = data_m2(17004:18797,4:20); % 1-16: Dest, 17: Pressure
[nRowsExp4,nColsExp4] = size(exp4);

% Find pylorus (most common minimum CSA):
for row = 1:nRowsExp4 % Pylorus == minimum CSA on EndoFLIP
    minimum = min(exp4(row,1:16)); % Save as index 18
    exp4(row,19) = find(exp4(row,1:16)==minimum,1);
    exp4(row,18) = minimum;
end
pyl_flip = mode(exp4(:,19)); % Find most common sensor with minCSA

% Shift by amount determined from probing:
shifter = 40; % The amount to shift the points in time (s)
exp4_shft = exp4(1+(shifter*Fs_flip):end,:);
[nRows, nCols] = size(exp4_shft);
t_exp4 = (0:nRows-1)*(1/Fs_flip); % Establish time vector

%% Slow Wave - Spike Correlations

BDFstartTime = 280; % (s)
gasFreq = 4.44; % Taken from collated data spreadsheet
intFreq = 19.49;

% Master List
SW_SPK_ML = [];
% IDX: [SWclst, SPKclst, meanStartATdiff, overlap(mm)]

SW_SPK_ML = SW_SPK_Correlations(SWdata,SPKdata,thresholds);

% Slow wave - Contraction Correlations

% Reconfigure electrode structures
pyl_elec = 13; % Manually determined using GEMS
g_clsts = [1,3,4,5,6,7,8,9];
i_clsts = [11:13;15:17;21:23;25:27;30:32;34:36;39:41;42:44];
[qz_AVE,qz_SD,pyl_elec] = QuiescentZoneFunc(SWdata,oldElecConfig,g_clsts,i_clsts);
qz_metrics = [qz_AVE,qz_SD];
[elecConfig,flipConfig] = RecentreElecConfigs(oldElecConfig,pyl_elec,pyl_flip);

i_clsts = 10:47; % Changing the list to incorporate all intestinal waves

% Extract valid contractions & metrics
[CONdata, ~] = ValidateContractions(flipConfig,excelCONdata,exp4_shft,t_exp4);
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
    metrics_temp = ContractionMetrics(CONclst,flipConfig,CONdata,exp4_shft,t_exp4);
    CON_SW_metrics = [CON_SW_metrics; metrics_temp];
end

% Contractions correlated with spikes
CONs = unique(CON_SPK_ML(:,1));
metrics_temp = [];
CON_SPK_metrics = [];
for CONclst = CONs'
    metrics_temp = ContractionMetrics(CONclst,flipConfig,CONdata,exp4_shft,t_exp4);
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

% FUNCTION TO INDENTIFY CON TIMINGS
[SPKprop,totals] = SpikePropagation(SPKdata,oldElecConfig,elecConfig);
RgnlPWR = PowerPerLocation(toappSWdata,elecConfig,oldElecConfig,gasFreq,intFreq);
AmpGrad = AmplitudeGradient(SWdata,toappSWdata,SPKdata,toappSPKdata,CONdata,elecConfig,oldElecConfig,flipConfig,BDFstartTime);
[ap_list,ap_AVE,dp_list,dp_AVE] = Pacemaker(toappSWdata,elecConfig,oldElecConfig,g_clsts,i_clsts,qz_AVE);
[PropOut] = Propagation(SWdata,toappSWdata,SPKdata,elecConfig,oldElecConfig,g_clsts,i_clsts,qz_AVE);


% Independent SPKs
% [idp_SPK_metrics] = independentSPKs(SPKdata,toappSPKdata,SW_SPK_ML,Fs_SPK,oldElecConfig,elecConfig,BDFstartTime);

% MasterList means and SDs:
SW_SPK_smry = [mean(SW_SPK_ML(:,3)),std(SW_SPK_ML(:,3)),mean(SW_SPK_ML(:,4)),std(SW_SPK_ML(:,4))];
SW_CON_smry = [mean(SW_CON_ML(:,3)),std(SW_CON_ML(:,3)),mean(SW_CON_ML(:,4)),std(SW_CON_ML(:,4))];
CON_SPK_smry = [mean(CON_SPK_ML(:,3)),std(CON_SPK_ML(:,3)),mean(CON_SPK_ML(:,4)),std(CON_SPK_ML(:,4))];
all_smry = [SW_SPK_smry,SW_CON_smry,CON_SPK_smry];

disp('DONE!');

%% Rewriting the ECC script to simplify the data analysis process

% GEMS parameters and outputs:
BDFstartTime = 280; % (s)
gasFreq = 4.44; % Taken from collated data spreadsheet
intFreq = 19.49;

% Reconfigure electrode structures:
pyl_elec = 13; % Manually determined using GEMS
g_clsts = [1,3,4,5,6,7,8,9];
i_clsts = [11:13;15:17;21:23;25:27;30:32;34:36;39:41;42:44];
[qz_AVE,qz_SD,pyl_elec] = QuiescentZoneFunc(SWdata,oldElecConfig,g_clsts,i_clsts);
qz_metrics = [qz_AVE,qz_SD];
[elecConfig,flipConfig] = RecentreElecConfigs(oldElecConfig,pyl_elec,pyl_flip);
i_clsts = 10:47; % Changing the list to incorporate all intestinal waves

% Extract valid contractions & metrics:
[CONdata, ~] = ValidateContractions(flipConfig,excelCONdata,exp4_shft,t_exp4);

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


%% OLD CODE

% %% Quantify values from each region for comparison
% 
% 
% 
% % Define region boundaries (electrical clusters)
% gas_SWs = (1:9);
% gas_SPKs = (1:50);
% pyl_SPKs = (51:100);
% duo_SWs = (10:100);
% duo_SPKs = (101:200);
% 
% 
% 
% gas_SW_metrics = [];
% for gas_SW = gas_SWs
%     temp = [];
%     if ismember(gas_SW,SWdata(:,1)) % Check it exist
%         % Get metrics
%         temp = SlowWaveMetrics(gas_SW,SWdata,oldElecConfig,elecConfig);
%         % Save to permanent list
%         gas_SW_metrics = [gas_SW_metrics; temp];
%     end
% end
% 
% duo_SW_metrics = [];
% for duo_SW = duo_SWs
%     temp = [];
%     if ismember(duo_SW,SWdata(:,1)) % Check it exist
%         % Get metrics
%         temp = SlowWaveMetrics(duo_SW,SWdata,oldElecConfig,elecConfig);
%         % Save to permanent list
%         duo_SW_metrics = [duo_SW_metrics; temp];
%     end
% end
% 
% gas_SPK_metrics = [];
% for gas_SPK = gas_SPKs
%     temp = [];
%     if ismember(gas_SPK,SPKdata(:,1)) % Check it exist
%         % Get metrics
%         temp = SpikeMetrics(SPKdata,toappSPKdata,gas_SPK,Fs_SPK,oldElecConfig,elecConfig);
%         % Save to permanent list
%         gas_SPK_metrics = [gas_SPK_metrics; temp];
%     end
% end
% 
% duo_SPK_metrics = [];
% for duo_SPK = duo_SPKs
%     temp = [];
%     if ismember(duo_SPK,SPKdata(:,1)) % Check it exist
%         % Get metrics
%         temp = SpikeMetrics(SPKdata,toappSPKdata,duo_SPK,Fs_SPK,oldElecConfig,elecConfig);
%         % Save to permanent list
%         duo_SPK_metrics = [duo_SPK_metrics; temp];
%     end
% end
% 
% pyl_SPK_metrics = [];
% for pyl_SPK = pyl_SPKs
%     temp = [];
%     if ismember(pyl_SPK,SPKdata(:,1)) % Check it exist
%         % Get metrics
%         temp = SpikeMetrics(SPKdata,toappSPKdata,pyl_SPK,Fs_SPK,oldElecConfig,elecConfig);
%         % Save to permanent list
%         pyl_SPK_metrics = [pyl_SPK_metrics; temp];
%     end
% end
% 
% disp('DONE!');

%% Spikes & SWs not associated with contractions

% % list of spikes
% all_SPKs = unique(SPKdata(:,1));
% 
% tonic_SPK_metrics = [];
% 
% % for each spike
% for SPK = all_SPKs'
%     
%     % if not a associated with a contraction
%     if ~ismember(SPK,CON_SPK_ML(:,2))
%         
%         % get spike metrics
%         temp = SpikeMetrics(SPKdata,toappSPKdata,SPK,Fs_SPK,oldElecConfig,elecConfig);
%         tonic_SPK_metrics = [tonic_SPK_metrics; temp];
%     end
% 
% end
% 
% 
% % list of spikes
% all_SWs = unique(SWdata(:,1));
% 
% tonic_SW_metrics = [];
% 
% % for each spike
% for SW = all_SWs'
%     
%     % if not a associated with a contraction
%     if ~ismember(SW,SW_CON_ML(:,1))
%         
%         % get spike metrics
%         temp = SlowWaveMetrics(SW,SWdata,oldElecConfig,elecConfig);
%         % Amp
%         SW_AMP_idx = find(SW_AMPs(:,1) == SW); % Slow wave amplitude data
%         SW_AMP = SW_AMPs(SW_AMP_idx,2);
%         temp = [temp, SW_AMP];
%         tonic_SW_metrics = [tonic_SW_metrics; temp];
%     end
% 
% end


%% Spike propagation per organ

% % Find all spikes which touch pylorus and categorise their organ-level
% % propagation
% 
% % List of quiescent electrodes
% g_idx = find(SWdata(:,1) <= 8);
% g_elecs = unique(SWdata(g_idx,2))';
% 
% i_idx = find(SWdata(:,1) >= 9);
% i_elecs = unique(SWdata(i_idx,2))';
% 
% sw_elecs = [g_elecs, i_elecs];
% p_elecs = setxor(1:256,sw_elecs);
% 
% p_elecs2 = [121,117,113 109 105 101 99 122 118 114 110 106 102 100 123 119 115 111 107 103 124 120 116 112 108 104 134 133 135 205 209 213 217 221 225 206 210 214 218 222 226 203 207 211 215 219 223 227 200 204 208 212 216 220 224 229 230 234 235 244 245];
% 
% % Initialise tallies
% sd_tally = 0;
% d_tally = 0;
% s_tally = 0;
% iso_tally = 0;
% qz_tally = 0;
% 
% for SPK = unique(SPKdata(:,1))' % For each spike
%     
%     SPKidx = find(SPKdata(:,1) == SPK);
%     elecs = SPKdata(SPKidx,2)';
%     
%     
%     if (sum(ismember(elecs,p_elecs2)) > 0) % If it touches quiescent zone
%         qz_tally = qz_tally + 1;
%         
%         if (sum(ismember(elecs,g_elecs))>0) && (sum(ismember(elecs,i_elecs))>0) % If touches both
%             sd_tally = sd_tally + 1;
%             
%         elseif (sum(ismember(elecs,i_elecs))>0) % If touches duodenum
%             d_tally = d_tally + 1;
%             
%         elseif (sum(ismember(elecs,g_elecs))>0) % If touches stomach
%             s_tally = s_tally + 1;
%             
%         else % Else isolated to pylorus
%             iso_tally = iso_tally + 1;
%         end
%         
%     end
%     
% end
% 
% % Convert to percentages
% sd_perc = sd_tally / qz_tally * 100;
% d_perc = d_tally / qz_tally * 100;
% s_perc = s_tally / qz_tally * 100;
% iso_perc = iso_tally / qz_tally * 100;
% 
% sd_tally
% d_tally
% s_tally
% iso_tally
% 
% disp('DONE!');

%% Code archive

% % Reconfigure electrode structures
% pyl_elec = 13; % Manually determined using GEMS
% [elecConfig,flipConfig] = RecentreElecConfigs(oldElecConfig,pyl_elec,pyl_flip);
% 
% % Extract valid contractions & metrics
% [CONdata, validCON] = ValidateContractions(flipConfig,excelCONdata,exp4_shft,t_exp4);
% 
% CONmetrics = [];
% for CONclst = validCON
%     metrics = ContractionMetrics(CONclst,flipConfig,CONdata,exp4_shft,t_exp4);
%     CONmetrics = [CONmetrics; metrics];
% end
% 
% % Correlations between spikes and contractions:
% CON_SPK_ML = [];
% CON_SPK_ML2 = [];
% energy_temp = 0;
% energy_temp2 = 0;
% dur_temp = [];
% dur_temp2 = [];
% Fs_SPK = 512;
% 
% % Shift time values in SPKdata so it aligns with CONdata
% SPKdata_shft = SPKdata;
% SPKdata_shft(:,3:4) = SPKdata(:,3:4) - 280;
% 
% % For each validated contraction:
% for CONclst = validCON
%     
%     % Find temporally related spikes for this contraction (within 2 sec)
%     [CON_SPK_Tmprl] = ContSpkTemporal(CONdata,CONclst,SPKdata_shft);
%     % IDX: : 1=SPKclst, 2=meanAT(s), 3=meanATdiff(s).
%     
%     if ~isempty(CON_SPK_Tmprl) % If we have any valid spikes
%         SPK_clsts = CON_SPK_Tmprl(:,1);
%         for j = 1:length(SPK_clsts) % For each valid spike
%             
%             % Reject pairing if it does not overlap at least 10 mm
%             CON_SPK_Sptl = ContSpkSpatial(CONdata,CONclst,SPKdata_shft,SPK_clsts(j),oldElecConfig,elecConfig,flipConfig);
%             % IDX: 1=SPKclst, 2=propLength(mm), 3=duration(s),
%             % 4=1stATdiff(s), 5=overlap(mm), 6=minStartAT(s),
%             % 7=maxFinalAT(s)
%             
%             if ~isempty(CON_SPK_Sptl)
%                 % Save to master list
%                 CON_SPK_ML(end+1,1:7) = [CONclst,SPK_clsts(j),CON_SPK_Tmprl(j,3),CON_SPK_Sptl(2),CON_SPK_Sptl(3),CON_SPK_Sptl(4),CON_SPK_Sptl(5)];
%                 % IDX: 1=CONclst, 2=SPKclst, 3=meanATdiff,
%                 % 4=propLength(mm), 5=duration(s),6=1stATdiff(s),
%                 % 7=overlap(mm),8=aveAmp, 9=energy, 10=totalEnergy,
%                 % 11=totalDur
%                 
%                 % save to temp array 1st Start and Final End times
%                 dur_temp(end+1,:) = [CON_SPK_Sptl(6),CON_SPK_Sptl(7)];
%                 
%                 % Get spike metrics
%                 SPKmetrics = SpikeMetrics(SPKdata_shft,toappSPKdata,SPK_clsts(j),Fs_SPK);
%                 
%                 % Save aveAmp and energy to master list
%                 CON_SPK_ML(end,8:9) = SPKmetrics(2:3);
%                 
%                 % Save temporary energy so it can be summed
%                 energy_temp = energy_temp + SPKmetrics(3);
%                 
%                 if (CON_SPK_Sptl(4) < 0) % If spike occurs before contraction
%                     % Temp arrays for total duration and energy
%                     dur_temp2(end+1,:) = [CON_SPK_Sptl(6),CON_SPK_Sptl(7)];
%                     energy_temp2 = energy_temp2 + SPKmetrics(3);
%                     
%                     % Also save to a more exclusive master list
%                     CON_SPK_ML2(end+1,1:9) = CON_SPK_ML(end,1:9);
%                     
%                 end
%             end
%             
%         end
%         
%     end
%     if energy_temp ~= 0
%         CON_SPK_ML(end,10) = energy_temp;
%     end
%     if energy_temp2 ~= 0
%         CON_SPK_ML2(end,10) = energy_temp2;
%     end
%     energy_temp = 0; % Reset for next SPKclst
%     energy_temp2 = 0;
%     
%     % Save 1st Start and Final End times to master list
%     if ~isempty(dur_temp)
%         CON_SPK_ML(end,11) = max(dur_temp(:,2)) - min(dur_temp(:,1));
%     end
%     if ~isempty(dur_temp2)
%         CON_SPK_ML2(end,11) = max(dur_temp2(:,2)) - min(dur_temp2(:,1));
%     end
%     dur_temp = []; % Reset for next SPKclst
%     dur_temp2 = [];
%     
% end
% 
% 
% % Correlations between SW and contractions :::::::::::::::::::::::::::::::
% 
% SWclsts = unique(SW_SPK_ML(:,1));
% valid_CONlist = [];
% valid_CONlist2 = [];
% for i = 1:length(SWclsts) % For each SW
%     
%     % Find temporally associated contractions:
%     valid_CON = ContSwTemporal(SWclsts(i),SWdata,CONmetrics,SW_SPK_ML,CON_SPK_ML);
%     if ~isempty(valid_CON)
%         valid_CONlist = [valid_CONlist; valid_CON];
%     end
%     
%     % Find spatially associated contractions:
%     valid_CON2 = ContSwSpatial(CONdata,SWdata,valid_CON,oldElecConfig,elecConfig,flipConfig);
%     if ~isempty(valid_CON)
%         valid_CONlist2 = [valid_CONlist2; valid_CON2];
%     end
%     
% end
% 
% % Now to pull metrics for each of these validated SWs:
% SWclsts = unique(valid_CONlist2(:,1));
% SW_CON_metrics = [];
% for SWclst = SWclsts' % For each slow wave
%     [metrics] = SlowWaveMetrics(SWclst,SWdata,oldElecConfig,elecConfig);
%     SW_CON_metrics = [SW_CON_metrics; metrics];
% end
% 
% % Repeat for contractions:
% CONclsts = unique(valid_CONlist2(:,5));
% CONmetrics2 = [];
% for CONclst = CONclsts' % For each contraction
%     [metrics] = ContractionMetrics(CONclst,flipConfig,CONdata,exp4_shft,t_exp4);
%     CONmetrics2 = [CONmetrics2; metrics];
% end
% 
% % > What do we want to pull?
% % Get SW amp, dur, energy - toapp data
% % Get SW vertical propagation (same as SPKs)
% 
% % Then correlate:
% % > Relate SW amp to spike amp / dur / energy / number of burst correlated with contraction
% % > Relate SW amp to contraction amp / dur / energy
% % > Relate SW propagation to contractile propagation
% 
% 
% disp('DONE!');