function [dataFinal] = assignForceOutputTrcMot(data)
%Assign forces to a each foot and use this data to create .trc and .mot
%files for OpenSim
%   % Assign force to the feet and generate the .trc and .mot files. Data
%   must be a file generated from the function btk_loadc3d

          %Choose force plate filter frequency
          data.FilterFreq = 30;
          %Assign forces to a foot
          dataForcesAssigned = assign_forces(data,{'RCAL';'LCAL'},{'calcn_r';'calcn_l'},[30, 0.25],data.FilterFreq);
          %Create the .trc and .mot files
          dataFinal = btk_c3d2trc_treadmill_LS(dataForcesAssigned,'off');

end

