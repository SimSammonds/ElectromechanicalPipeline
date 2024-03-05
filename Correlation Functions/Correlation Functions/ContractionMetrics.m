function [metrics] = ContractionMetrics(CONclst,flipConfig,CONdata,EXPdata,time)
% ContractionMetrics is a function to pull useful data from GEMS' toapp
% data for a specific contraction.
% metrics idx:
% 1=clst, 2=propLength, 3=aveAmpl, 4=totalPower(?), 5=totalDuration,
% 6=1stStartTime, 7=direction

% Author: Sam Simmonds
% Date: 14th Novemember 2022

% Define tolerance
tol = 0.01;
% Find relevant indices for relation to CONdata
idx = find(CONdata(:,1) == CONclst);

% Find propagation distance using elec numbers
elecs = CONdata(idx,2);
propLength = flipConfig(min(elecs)) - flipConfig(max(elecs));

aveAmpl = mean(CONdata(idx,7));
energy = 0;
totalDuration = max(CONdata(idx,5)) - min(CONdata(idx,3));
initialStartTime = min(CONdata(idx,3));
finalEndTime = max(CONdata(idx,5));

% Direction
if initialStartTime == CONdata(idx(1),3)
    direction = -1; % Then propagates downwards
elseif initialStartTime == CONdata(idx(end),3)
    direction = 1; % Then propagates upwards
else
    direction = 0; % Propagates in both directions
end

% Energy calculation
for j = 1:length(elecs)
    % Find start and end times
    startTime = find(abs(time-CONdata(idx(j),3))<tol);
    endTime = find(abs(time-CONdata(idx(j),5))<tol);
    
    % Create sig and time vectors for integration
    sig = EXPdata(startTime:endTime,elecs(j));
    sig2 = sig.^2;
    T = time(startTime:endTime);
    
    % Perform energy integration
    energy = trapz(T,sig2') + energy;
end

metrics = [CONclst,propLength,aveAmpl,energy,totalDuration,...
    initialStartTime,finalEndTime];


end

