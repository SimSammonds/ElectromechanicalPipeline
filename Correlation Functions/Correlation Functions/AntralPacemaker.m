function [ap_list,ap_AVE,dp_list,dp_AVE] = Pacemaker(toappData,elecConfig,oldElecConfig,g_clts,i_clts,QZWidth)
%AntralPacemaker is a function which will identify the location of any
%antral pacemakers for a marked SW file.
% Inputs: toappData - SW data, elecConfig - obv, g_clts = gastric wave
%         clsts
% Outputs: AT - activation times, R - elec location ROW, C - "" COLUMN

% Author: Sam Simmonds
% Date: 21th July 2023

% Find the pyloric channels:
[~,~,p_list] = TrueChans(toappData,g_clts);
[p_list] = ConfirmPlist(p_list,elecConfig,oldElecConfig,QZWidth);

% Define outputs (to prevent errors when not assigned)
OUT = [];
AVE = [];
SD = [];

% For each wave, find earliest AT:
for waveNum = g_clts
    TEMP = []; % Clear the temp outputs
    
    % Pull and sort markers for this wave:
    wave_data = toappData.toapp.TimeAmplCluster{1,waveNum};
    elecs = wave_data(:,4);
    ATs = wave_data(:,3);
    [ATs,idx] = sort(ATs,'ascend');  % Sort ATs by ascending order
    elecs = elecs(idx);
    
    % Check if the first 4 elecs are pacemakers for this wave
    for i = 1:4
        n_list = FindNeighbours(elecs(i),oldElecConfig);
        NumValidNbrs = 0;
        if ~isempty(n_list)
            %NumValidNbrs = sum(n_list<=ATs(i)); % This doesnt make sense
            % Compare adjacent AT(s)
            
            for j = 1:length(n_list) % Each neighbour
                % Check this elec has an associated AT for this wave
                n_idx = find(elecs == n_list(j));
                n_AT = ATs(n_idx);
                if n_AT >= ATs(i)
                    NumValidNbrs = NumValidNbrs + 1;
                end
            end %for
            if NumValidNbrs >= 4
                % We have found our pacemaker
                TEMP = [waveNum,elecs(i),ATs(i)];
            end %if
        end %if
        
        if ismember(waveNum,TEMP)
            break % End early if we find our pacemaker
        end
    end %for
    
    % Find location relative to QZ electrodes:
    if ~isempty(TEMP)
        [~,~,dist] = nearestElecs(TEMP(2), p_list, oldElecConfig);
        OUT = [OUT; waveNum, dist, elecs(i), ATs(i)];
    end
end

if ~isempty(OUT)    % Quantify mean +/- SD
    AVE = mean(OUT(:,2));
    SD = std(OUT(:,2));
    ap_AVE = [AVE,SD];
end

end

function [g_list,i_list,p_list] = TrueChans(T_data,g_clts)

numWaves = length(T_data.toapp.TimeAmplCluster) - 1;
TruthData = T_data.toapp.TimeAmplCluster;

g_list = [];
i_list = [];
m_elec = [];

for wave = 1:numWaves
    [R,~] = size(TruthData{1,wave});
    
    for row = 1:R   % for each electrode that this wave passes
        elec = TruthData{1,wave}(row,4);
        amp = TruthData{1,wave}(row,7);
        m_elec = [m_elec, elec]; % All marked elecs
        
        if ismember(wave,g_clts)
            g_list = [g_list, elec]; % Save elec to g_elecs
        else
            i_list = [i_list, elec];
        end
    end
end

g_list = unique(g_list);
i_list = unique(i_list);
m_elec = unique(m_elec);

all_elec = 1:256;
p_list = all_elec(~ismember(all_elec,m_elec));

end

function [p_list2] = ConfirmPlist(p_list,elecConfig,oldElecConfig,QZWidth)
% A quick function to confirm that electrodes within our p_list are not
% simply unmarked electrodes within the gastric or intestinal regions. This
% function will rely on the QZ width input to vary tolerance.
p_list2 = [];

for i = 1:length(p_list) % For every elec in p-List
    % Pull elec R and C
    [R,C] = find(oldElecConfig == p_list(i));
    % Check that this elec is within the QZ region
    dist = elecConfig(R,C);
    if abs(dist) <= (round(QZWidth/10))*5
        % Save to new list
        p_list2 = [p_list2; p_list(i)];
    end
    
end
end

function [n_list] = FindNeighbours(elec,oldElecConfig)
% A function to find 8 neighbouring electrodes
% Try not to error when on a boundary!

n_list = [];
out = [];

% Find indices of current elec
[r,c] = find(oldElecConfig == elec);

% We know that a neighbour can be rejected if:
    % Neighbour is = 257
    % r < 1
    % c < 1 or c > 29
    
% Find 8 neighbours:
if ((r - 1) >= 1) % North
    n_list = [n_list; oldElecConfig(r-1,c)];
    if (c - 1) >= 1 % Northwest
        n_list = [n_list; oldElecConfig(r-1,c-1)];
    end
    if (c + 1) <= 29 % Northeast
        n_list = [n_list; oldElecConfig(r-1,c+1)];
    end
end

if ((r + 1) >= 1) % South
    n_list = [n_list; oldElecConfig(r+1,c)];
    if (c - 1) >= 1 % Southwest
        n_list = [n_list; oldElecConfig(r+1,c-1)];
    end
    if (c + 1) <= 29 % Southeast
        n_list = [n_list; oldElecConfig(r+1,c+1)];
    end
end
        
if ((c - 1) >= 1) % West
    n_list = [n_list; oldElecConfig(r,c-1)];
end
if ((c + 1) >= 1) % East
    n_list = [n_list; oldElecConfig(r,c+1)];
end

% Now go ahead and reject all neighbours which are == 257
if ~isempty(n_list)
    for i = 1:length(n_list)
        if n_list(i) ~= 257
            out = [out; n_list(i)];
        end
    end
end

clear n_list
n_list = out;

end

function [gElec, iElec, dist] = nearestElecs(gElecs, iElecs, elecConfig)
% nearestElecs is a function to find the closest pair of electrodes from a
% list. Uses the magnitude method, rather than the vertical method.

% Find gastric indices
gMaster = []; % Establish gastric masterlist
for m = 1:length(gElecs)
    % find x position
    [r,c] = find(elecConfig == gElecs(m));
    gMaster = [gMaster; gElecs(m), r, c];
end

% Find intestinal indices
iMaster = [];
for n = 1:length(iElecs)
    % find x position
    [r,c] = find(elecConfig == iElecs(n));
    iMaster = [iMaster; iElecs(n), r, c];
end

% Find shortest distance
best = 999.9;
for k = 1:length(gElecs)
    g_x = gMaster(k,2);
    g_y = gMaster(k,3);
    
    for l = 1:length(iElecs)
        % Calculate magnitude distance
        x = g_x - iMaster(l,2);
        y = g_y - iMaster(l,3);
        mag = sqrt((x*x) + (y*y));
        if mag < best
            best = mag;
            gElec = gMaster(k,1);
            iElec = iMaster(l,1);
        end
    end
end

dist = best*5; % (mm)

end










