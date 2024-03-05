function [OUT] = CONmetricsForHisto(CONdata)
% This function reads in CONdata and outputs averages for each location, as
% well as a total average

% for all gastric; col 9 == 1;
gasIDX = CONdata(:,9) == 1;
% for all pyloric; col 9 == 0;
pylIDX = CONdata(:,9) == 0;
% for all duodenal; col 9 == -1;
duoIDX = CONdata(:,9) == -1;

gasAVE = mean(CONdata(gasIDX,8));
pylAVE = mean(CONdata(pylIDX,8));
duoAVE = mean(CONdata(duoIDX,8));
totalAVE = mean(CONdata(:,8));

OUT = [gasAVE, pylAVE, duoAVE, totalAVE];

end