function [validCON2] = ContSwSpatial(CONdata,SWdata,validCON,oldElecConfig,elecConfig,flipConfig)
% ContSwSpatial is a function to determine which slow waves are spatially
% correlated with contractions.
% A list of potential candidates is fed in as validCON.
% validCON idx: 1=SWclst, 2=SW_AT, 3=SPKclst, 4=SPK_AT, 5=CONclst, 6=CON_AT,
% 7=startATdiff(CON and SW), 8=overlap.
% Author: Sam Simmonds
% Date: 9th Novemember 2022

validCON2 = [];
[r,~] = size(validCON);
thresh = 10; % (mm)

for i = 1:r % For each SW-CON pair
    
    SWclst = validCON(i,1); % Clusters
    SWidx = find(SWdata(:,1) == SWclst); % Indices
    SWelecs = SWdata(SWidx,2); % Electrodes
    
    SPKclst = validCON(i,3);
    
    for j = 1:length(SWelecs) % For each electrode
        [elec_r, elec_c] = find(oldElecConfig == SWelecs(j)); % Find co-ords
        SWdist(j) = elecConfig(elec_r,elec_c); % Get distance
    end
    
    CONclst = validCON(i,5); % Repeat for CON
    CONidx = find(CONdata(:,1) == CONclst);
    CONelecs = CONdata(CONidx,2);
    CONdist = flipConfig(CONelecs);
    
    % Times for export metrics
    CON_AT = min(CONdata(CONidx,3));
    SW_AT = min(SWdata(SWidx,3)) - 280; % 0-120s
    
    % Find overlap (mm)
    min_list = [min(SWdist), min(CONdist)];
    max_list = [max(SWdist), max(CONdist)];
    overlap = (min(max_list) - max(min_list));
    
    % Does it surpass threshold?
    if overlap >= thresh
        % Calculate additional metrics
        
        startATdiff = CON_AT - SW_AT;
        % Output
        validCON2(end+1,:) = [SWclst,validCON(i,2),SPKclst,validCON(i,4),CONclst,validCON(i,6),startATdiff,overlap];
    end
    
end

end

