function [Mean, SD, Upper, Lower] = Mean_SD_3D(Variable)
%Calculate mean, SD, and upper and lower bounds of the input variable

Mean = nanmean(Variable, 2);
SD = nanstd(Variable, 1, 2);
Upper = Mean+SD;
Lower = Mean-SD;
end

