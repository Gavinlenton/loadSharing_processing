function [dataFinal] = assignForceOutputTrcMot(fileName, pname)
%Assign forces to a each foot and use this data to create .trc and .mot
%files for OpenSim
%   Detailed explanation goes here

     % Assign force to the feet and generate the .trc and .mot files
          %Load the new acquisition
          data1 = btk_loadc3d([pname, fileName], 50);
          %Choose force plate filter frequency
          data1.FilterFreq = 30;
          %Assign forces to a foot
          data1 = assign_forces(data1,{'RCAL';'LCAL'},{'calcn_r';'calcn_l'},[30, 0.25],data1.FilterFreq);
          %Create the .trc and .mot files
          dataFinal = btk_c3d2trc_treadmill(data1,'off');

end

