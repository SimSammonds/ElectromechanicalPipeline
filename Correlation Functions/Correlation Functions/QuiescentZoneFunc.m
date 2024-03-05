function [AVE,SD,pyl_elec] = QuiescentZoneFunc(SWdata,elecConfig,gas_list,int_list)
%QuiescentZoneFunc is a function to extract all valuable quiescent zone data from a SW GEMS file
% Inputs: SWdata (DataOrganised file), elecConfig (matrix), gas_list
% (matrix), int_list (matrix)
% gas_list and int_list have the structure: each row == unique gastric SW
% each col; gas == singular gas SW, int == 2-3 associated int SWs
% Outputs: AVE (average QZ width), SD (standard deviation of QZ width),
% pyl_elec (row position associated with middle of quiescent zone ~pylorus)

width = [];
g_elecs = [];
i_elecs = [];

% For each gastric SW
for i = 1:length(gas_list)
    % Find all gas elec
    idx = find(SWdata(:,1) == gas_list(i)); % Idx for current gastric SW
    gElecs = SWdata(idx,2);
    
    % For each paired int SW
    for j = 1:length(int_list(i,:))
        % Find all int elecs
        idx = find(SWdata(:,1) == int_list(i,j)); % Idx for current int SW
        iElecs = SWdata(idx,2);
        
        % Save the smallest distance from all these combinations
        [gElec, iElec, dist] = nearestElecs(gElecs,iElecs,elecConfig);
        width = [width, dist];
        g_elecs = [g_elecs, gElec];
        i_elecs = [i_elecs, iElec];

    end
end

% Mean +/- SD width
AVE = mean(width);
SD = std(width);
MIN = min(width);
g_elec = mode(g_elecs);
i_elec = mode(i_elecs);

% Find idx
[gr,gc] = find(elecConfig == g_elec);
[ir,ic] = find(elecConfig == i_elec);

% Find middle
ideal_r = round((ir-gr)/2);
pyl_elec = ideal_r + gr; % Row position associated with middle of pylorus

%-------------------------------------------------------------------------
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

end


% I know the boundaries, and can easily find the row / col positions. I can
% find the average middle row, and go from there? I can find the average
% middle row PER column, which may allow me to create a non-horizontal QZ.
% I like this idea (although it complicates everything else that follows!)
% ^^^ Let's try it later


































