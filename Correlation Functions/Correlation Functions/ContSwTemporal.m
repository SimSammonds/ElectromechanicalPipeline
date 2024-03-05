function [valid_CON] = ContSwTemporal(SWclst,SWdata,CONmetrics,SW_SPK_ML,CON_SPK_ML)
% RelateSWsContractions is a function to find the relationships between
% slow waves, spikes, and contractions.
% valid_CON idx:
% 1=SWclst, 2=SPKclst, 3=CONclst, 4=SW_AT, 5=SPK_AT, 6=CON_AT
% Pass in CON_SPK_ML2, so all CON_SPK correlations have SWs occuring
% prior to contractions.

% Author: Sam Simmonds
% Date: 9th Novemember 2022

valid_CON = [];

% Pull each associated spike
SPKidx = find(SW_SPK_ML(:,1) == SWclst);
SPKclsts = SW_SPK_ML(SPKidx,2);

% Pull SW AT
SWidx = find(SWdata(:,1) == SWclst);
SW_AT = min(SWdata(SWidx,3)) - 280; % 0-120s

for i = 1:length(SPKclsts) % For each spike correlation
    % Check each associated contraction
    CONidx = find(CON_SPK_ML(:,2) == SPKclsts(i));
    CONclsts = CON_SPK_ML(CONidx,1);
    
    SPK_AT = SW_SPK_ML(SPKidx(i),3) + SW_AT;
    
    if ~isempty(CONclsts)
        for j = 1:length(CONclsts)
            % Pull CON AT
            CONidx2 = find(CONmetrics(:,1) == CONclsts(j));
            CON_AT = CONmetrics(CONidx2,6);
            
            % Check if SWs occur first, and if they are within threshold
            ATdiff = CON_AT - SW_AT; % Pos = SW activates first
            
            if (ATdiff >= 0) && (ATdiff <= 4)
                % Valid - export
                valid_CON(end+1,:) = [SWclst,SW_AT,SPKclsts(i),SPK_AT,CONclsts(j),CON_AT];
            end
        end
    end
end

end

