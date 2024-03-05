function [CONdata, valid_conts2] = ValidateContractions(flipConfig,CONdata,EXPdata,time)
% ContractionMetrics aims to validate contractions by checking amplitude and propagation distance

% Author: Sam Simmonds
% Date: 8th Novemember 2022

% Checking the number of occluding contractions:
numMin = [];

% Amplitude metrics and threshold checking:
valid_conts1 = [];
tol = 0.01;
ampl_thresh = 1; % Threshold for valid contraction (mm)
for i = 1:length(CONdata(:,1))
    
    % Find start and end times
    startTime = find(abs(time-CONdata(i,3))<tol);
    endTime = find(abs(time-CONdata(i,5))<tol);
    
    % Save duration (to be averaged later)
    duration = endTime - startTime;
    
    % Find relevant flip electrode
    elec = CONdata(i,2);
    
    % Determine if this elec is gastric, pyl, or duodenal
    if flipConfig(elec) == 0 % pyloric
        location = 0;
    elseif flipConfig(elec) > 1 % gastric
        location = 1;
    else % duodenal
        location = -1;
    end
    
    % Extract range of relevant data for this electrode
    cont = EXPdata(startTime:endTime,elec);
    ampl = max(cont) - min(cont); % To be averaged later
    
    % Save all additional data, even if current elec/cont is rejected
    CONdata(i,7:9) = [ampl, duration, location];
    
    if ampl >= ampl_thresh  % save clst and elec labels
        valid_conts1(end+1,:) = [CONdata(i,1),CONdata(i,2)];

        % Save a count of the number of occluding contractions:
        if min(cont) <= 4.7; numMin = [numMin,CONdata(i,1)]; end
    end
end

% Locational checking and propagation thresholding
prop_thresh = 2; % 2 electrodes
valid_conts2 = [];
for i = 1:length(valid_conts1) % Check if surpasses propagation threshold
    cur_clst = valid_conts1(i,1);
    if numel(find(valid_conts1(:,1) == cur_clst)) >= prop_thresh
        valid_conts2(end+1) = [cur_clst]; % Valid
    end
end

valid_conts2 = unique(valid_conts2); % Make unique

% Disp the number of occluding contractions:
disp('Num occluding contractions: ')
disp(numel(unique(numMin)))
disp('Total CONs: ')
disp(numel(unique(CONdata(:,1))))
end
