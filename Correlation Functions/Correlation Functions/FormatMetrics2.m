function [DataOut] = FormatMetrics2(Type1,Type2,metrics1,metrics2,ML)
% FormatMetrics is a function to take correlations and metrics, and to
% organise them in such a way that they can be easily exported and
% immediately plotted in excel - save time reorganising data later.
% Type1 should be either 'SW' or 'CON', and Type2 should be 'SPK' or 'CON'

% Inputs: metrics1/2: [clst, DUR, AMP, NRG]

% Outputs: DataOut: [Dur1, Amp1, Nrg1, Dur2, Amp2, Nrg2]


% Author: Sam Simmonds
% Date: 15th November 2022
[r,~] = size(ML);
DataOut = NaN(r,6);

for i = 1:length(ML(:,1)) % For every entry in the ML
    DUR_1 = metrics1(metrics1(:,1) == ML(i,1),2);
    AMP_1 = metrics1(metrics1(:,1) == ML(i,1),3);
    NRG_1 = metrics1(metrics1(:,1) == ML(i,1),4);
    DUR_1 = DUR_1(1); % Remove duplicates
    AMP_1 = AMP_1(1);
    NRG_1 = NRG_1(1);

    DUR_2 = metrics2(metrics2(:,1) == ML(i,2),2);
    AMP_2 = metrics2(metrics2(:,1) == ML(i,2),3);
    NRG_2 = metrics2(metrics2(:,1) == ML(i,2),4);
    DUR_2 = DUR_2(1); % Remove duplicates
    AMP_2 = AMP_2(1);
    NRG_2 = NRG_2(1);

    DataOut(i,:) = [DUR_1, AMP_1, NRG_1, DUR_2, AMP_2, NRG_2];
end

end

