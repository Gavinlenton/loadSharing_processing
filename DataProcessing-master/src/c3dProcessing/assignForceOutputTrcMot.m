function [dataFinal, force_data2] = assignForceOutputTrcMot(data, fileName)
%Assign forces to a each foot and use this data to create .trc and .mot
%files for OpenSim
%   % Assign force to the feet and generate the .trc and .mot files. Data
%   must be a file generated from the function btk_loadc3d

% Choose force plate filter frequency
data.FilterFreq = 26;

% Combine force measures from both plates
%           dataForcesAssigned2 = assign_forces_LS(data,{'RCAL';'LCAL'},{'calcn_r';'calcn_l'},[20, 0.25],data.FilterFreq);

dataForcesAssigned3 = assign_forces_Gerber_method(data,[20, 0.25]);

% Create the .trc and .mot files
%           [dataFinal, force_data2] = btk_c3d2trc_treadmill_LS(dataForcesAssigned,'off');

%           [dataFinal, force_data2] = btk_c3d2trc_treadmill_LS_new(dataForcesAssigned2,'off');

if dataForcesAssigned3.GRF.FP(1).F(200,3) ~= 0
	[dataFinal, force_data2] = btk_c3d2trc_treadmill_LS_new(dataForcesAssigned3,'off');
	
else
	fprintf('Trial %s does has dodgy FP data, not writing it for further processing', fileName);
	dataFinal = dataForcesAssigned3;
	force_data2 = [];
end
% Alternative if things aren't working
%           [dataFinal] = btk_c3d2trc(dataForcesAssigned,'off');

end

