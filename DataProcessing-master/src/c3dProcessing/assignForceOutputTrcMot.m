function [dataFinal, force_data2] = assignForceOutputTrcMot(data, croppedTrialNum)
%Assign forces to a each foot and use this data to create .trc and .mot
%files for OpenSim
%   % Assign force to the feet and generate the .trc and .mot files. Data
%   must be a file generated from the function btk_loadc3d

          %Choose force plate filter frequency
          data.FilterFreq = 30;
          
          %Assign forces to a foot
%           dataForcesAssigned = assign_forces(data,{'RCAL';'LCAL'},{'calcn_r';'calcn_l'},[30, 0.25],data.FilterFreq);
          
          dataForcesAssigned2 = assign_forces_LS(data,{'RCAL';'LCAL'},{'calcn_r';'calcn_l'},[30, 0.25],data.FilterFreq);
          
          %Create the .trc and .mot files
%           [dataFinal, force_data2] = btk_c3d2trc_treadmill_LS(dataForcesAssigned,'off');
          
          [dataFinal, force_data2] = btk_c3d2trc_treadmill_LS_new(dataForcesAssigned2,'off', croppedTrialNum);
          
          % Alternative if things aren't working
%           [dataFinal] = btk_c3d2trc(dataForcesAssigned,'off');

end

