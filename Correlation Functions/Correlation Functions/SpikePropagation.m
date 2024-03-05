function [SPKprop,totals] = SpikePropagation(SPKdata,oldElecConfig,elecConfig)
% SpikePropagation is a function to determine qualitative propagation
% properties of a group of spike bursts.
% Author: Sam Simmonds
% Date Created: 3rd April 2023

% Inputs: SPK data set, elecConfig, elecConfigDistance
% Ouputs: Matrix of spike origin and propagation direction + organs
% traversed

SPKprop = [];

% For each SPK burst, identify origin
for spk = unique(SPKdata(:,1))'
    
    destination = ["","",""];
    
    % Find start time and electrode of origin 
    idx = find(SPKdata(:,1) == spk); % All current spk idx
    firstAT = min(SPKdata(idx,3));
    origin_idx = find(SPKdata(idx,3) == firstAT); % Origin spk idx
    elec = SPKdata(idx(origin_idx),2);
    
    % Classify electrode - organ
    % Find the indeces in normal elecConfig
    if length(elec) > 1
        elec = elec(1);
    end
    [r,c] = find(oldElecConfig == elec);
    
    if abs(elecConfig(r,c)) <= 5
        % Pyloric
        origin = "pyloric";
    elseif r == 1               % gastric
        origin = "gastric";
    elseif isnan(elecConfig(r-1,c))
        origin = "gastric";
    elseif (elecConfig(r-1,c) > 0)
        origin = "gastric";
    else
        % Duodenal
        origin = "duodenal";
    end
    
    % Classify propagation - organ (regardless of origin)
    for i = 1:length(idx) % For each time this spk was measured
        elec2 = SPKdata(idx(i),2); % Find electrode
        [r,c] = find(oldElecConfig == elec2);
        
        if abs(elecConfig(r,c)) <= 5 % Pyloric
            destination(1) = "pyloric";
        elseif r == 1               % gastric
            destination(2) = "gastric";
        elseif isnan(elecConfig(r-1,c))
            destination(2) = "gastric";
        elseif (elecConfig(r-1,c) > 0)
            destination(2) = "gastric"; % Gastric
        else    % Duodenal
            destination(3) = "duodenal";
            
        end
    end
    
    SPKprop = [SPKprop; num2str(spk), origin, destination(1), destination(2), destination(3)];
    
    
end

    % Now find number that are:
    total = 0;
    allthree = 0;
    pylandstom = 0;
    pylandduo = 0;
    isolated = 0;
    
    for i = 1:length(SPKprop(:,1))
        
        if (SPKprop(i,3) == "pyloric") % If present in pylorus
            total = total + 1;

            if (SPKprop(i,4) == "gastric") && (SPKprop(i,5) == "duodenal")
                % In all three regions
                allthree = allthree + 1;
            elseif (SPKprop(i,4) == "gastric") % In pylorus and stomach
                pylandstom = pylandstom + 1;
            elseif (SPKprop(i,5) == "duodenal") % In pylorus and duodenum
                pylandduo = pylandduo + 1;
            else % Isolated to pyloruss
                isolated = isolated + 1;
            end
        end
    end
    
    % Total in each category
    totals = [total, allthree,pylandstom,pylandduo,isolated];
    

end















