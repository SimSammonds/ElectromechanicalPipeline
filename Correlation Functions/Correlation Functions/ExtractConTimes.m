function [OUT] = ExtractConTimes(c1,c2)
% This function will take the saved cursor_info and pull out the start or
% end contraction times

times = zeros(length(c1),2);

for i = 1:length(c1)
    temp = c1(i).Position;
    times(i,1) = temp(1);
    
    elecLabel(i) = str2num(c1(i).Target.Parent.YLabel.String);
end

% Match end times
for i = 1:length(c2)
    
    temp = c2(i).Position;
    % Find appropriate electrode
    idx = find(times(:,1) == str2num(c2(i).Target.Parent.YLabel.String));
    times(i,2) = temp(1); % Save end time
    
end

% Combine
OUT = [elecLabel',times(:,1),zeros(length(c1),1),times(:,2)];
OUT = sortrows(OUT);

end

