function [Mean, SD, Upper, Lower] = Mean_SD(Variable)
%Calculate mean, SD, and upper and lower bounds of the input variable

Mean = nanmean(Variable)';
SD = nanstd(Variable)';
Upper = Mean+SD;
Lower = Mean-SD;
end

