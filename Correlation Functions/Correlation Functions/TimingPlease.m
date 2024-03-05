function [OUT,BDFzeroes] = TimingPlease(tZero,time,TempCor,BDFend)
% I need a function which can automatically find the appropriate FLIP and
% BDF times for prokinetics and vagotomy, relative to the infusion start
% time or surgical time.

% Inputs:
% tZero = scalar; infusion/surgery time in FLIP entries.
% time = vector; times to output; relative to infusion/surgery time; minutes,
% TempCor = Nx2 matrix; BDF times linked to FLIP entry times -> do this manually before calling function.
% BDFend = Nx2 matrix; [exp#, BDF end time (BDf timing = s)] ;

% Outputs: [exp#, BDFstart, BDFend, EntryStart, EntryEnd, T-0 (min)]

% Author: Sam Simmonds
% Date: 23rd June 2023

Fs_flip = 10;   % (Hz)
OUT = [];

% First, find BDF zeros and exp range in flip entries:
BDFzeroes = zeros(height(TempCor),1);
exp_entries = zeros(height(TempCor),2);

for i = 1:length(TempCor)   % for each entry in TempCor
    BDFzeroes(i) = TempCor(i,2) - TempCor(i,1)*Fs_flip; % start times (entry)
    
    exp_entries(i,1) = BDFzeroes(i);    % Start for exp i (entry)
    exp_entries(i,2) = BDFzeroes(i)+ (BDFend(i,2)*Fs_flip);       % End for exp i (entry)
end

% Find associated BDF file for each time stamp:
for j = 1:length(time)
    
    % Find time in entries
    time_entries(j) = tZero + (time(j)*60*Fs_flip); % (entry)
    
    % Find which domain this is covered by
    for i = 1:height(exp_entries)
        if (time_entries(j) >= exp_entries(i,1)) && (time_entries(j) <= exp_entries(i,2))
            OUT(j,1) = BDFend(i,1);   % exp#
            OUT(j,2) = (time_entries(j) - BDFzeroes(i))/Fs_flip; % BDF start (s)
            OUT(j,3) = OUT(j,2) + 80; % BDF end (s)
            OUT(j,4) = time_entries(j); % CON start (entry)
            OUT(j,5) = OUT(j,4) + 80*Fs_flip; % CON end (entry)
            OUT(j,6) = time(j); % time stamp (minutes)
        end
    end
    if isempty(time_entries(j))
        disp('WARNING! NO APPROPRIATE BDF FOR THIS TIME STAMP')
    end
end % Repeat

end

