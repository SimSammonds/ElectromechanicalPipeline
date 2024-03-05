function [DataOut] = mergeMetrics2(M1,M2,M3,M4,M5,M6,M7,M8,M9)
% mergeMetrics is a function to format only the useful (correlatable)
% metrics from each script.
% Date: 23rd March 2023
% Inputs: from formatMetrics

[r(1),~] = size(M1);
[r(2),~] = size(M2);
[r(3),~] = size(M3);
[r(4),~] = size(M4);
[r(5),~] = size(M5);
[r(6),~] = size(M6);
[r(7),~] = size(M7);
[r(8),~] = size(M8);
[r(9),~] = size(M9);

DataOut = NaN(max(r),54); % Establish matrix

if ~isempty(M1) % Isotopic / Gastric
    DataOut(1:r(1),1:6) = M1(:,1:6);     % _DUR/_AMP, _DUR/_AMP
end
if ~isempty(M2) % Longitudinal / Gastric
    DataOut(1:r(2),7:12) = M2(:,1:6);     % _DUR/_AMP, _DUR/_AMP
end
if ~isempty(M3) % Circular / Gastric
    DataOut(1:r(3),13:18) = M3(:,1:6);     % _DUR/_AMP, _DUR/_AMP
end
if ~isempty(M4) % Isotopic / Intestinal
    DataOut(1:r(4),19:24) = M4(:,1:6);     % S_DUR/_AMP, _DUR/_AMP
end
if ~isempty(M5) % Longitudinal / Intestinal
    DataOut(1:r(5),25:30) = M5(:,1:6);     % _DUR/_AMP, _DUR/_AMP
end
if ~isempty(M6) % Circular / Intestinal
    DataOut(1:r(6),31:36) = M6(:,1:6);     % _DUR/_AMP, _DUR/_AMP
end
if ~isempty(M7) % Long / CON
    DataOut(1:r(7),37:42) = M7(:,1:6);     % _DUR/_AMP, _DUR/_AMP
end
if ~isempty(M8) % Circ / CON
    DataOut(1:r(8),43:48) = M8(:,1:6);     % _DUR/_AMP, _DUR/_AMP
end
if ~isempty(M9) % Isot / CON
    DataOut(1:r(9),49:54) = M9(:,1:6);     % _DUR/_AMP, _DUR/_AMP
end

end







