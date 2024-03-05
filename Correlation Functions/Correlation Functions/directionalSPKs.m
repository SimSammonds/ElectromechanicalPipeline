function [isotSPKs,longSPKs,circSPKs] = directionalSPKs(SPKdata,oldElecConfig,tol)
% directionalSPKs is a function to extract the clusters, and label them as
% LONGITUDINAL, CIRCULAR, or ~ISOTOPIC
% Use this function to test the hypothesis that SPKs occur in predominantly
% longitudinal or circumferential patches.

% INPUTS:
% SPK data; 1=SPKclst, 2=elec, 3=StartAT, 4=EndAT

% OUTPUTS:
% longSPKs; 1=SPKclst, 2=elec, 3=StartAT, 4=EndAT
% circSPKs; 1=SPKclst, 2=elec, 3=StartAT, 4=EndAT
% isotSPKs; 1=SPKclst, 2=elec, 3=StartAT, 4=EndAT

isotSPKs = [];
longSPKs = [];
circSPKs = [];

clsts = unique(SPKdata(:,1));
for clst = clsts'
    r_list = [];
    c_list = [];

    % Find all instances of clst
    idx = SPKdata(:,1) == clst;
    elecs = SPKdata(idx,2);
    for i = 1:length(elecs)
        [r,c] = find(oldElecConfig == elecs(i));
        r_list = [r_list, r];
        c_list = [c_list, c];
    end
    
    % Find if max long or circ prop is higher
    longProp = max(r_list)-min(r_list);
    circProp = max(c_list)-min(c_list);

    if ((longProp-tol) <= circProp) && ((circProp-tol) <= longProp)
        % More or less isotopic
        isotSPKs = [isotSPKs; SPKdata(idx,:)];

    elseif longProp > circProp  % Longitudinal
        longSPKs = [longSPKs; SPKdata(idx,:)];

    else % Circular
        circSPKs = [circSPKs; SPKdata(idx,:)];
    end
end

end