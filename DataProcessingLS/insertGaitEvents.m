function insertGaitEvents(acquisition, RightHS, LeftHS, RightTO, LeftTO)
%Insert the heel-strike events into the c3d file
%   Input the column vectors corresponding to the heel-strike events for
%   each foot and use these time points to add heel-strikes into the c3d
%   file

% Loop through all heel strikes and append event into the file
for i = 1:length(RightHS)
     
     btkAppendEvent(acquisition,'Foot Strike', RightHS(i,1), 'Right', '', 'Right Foot Strike', 1);
     btkAppendEvent(acquisition,'Foot Off', RightTO(i,1), 'Right', '', 'Right Foot Off', 2);
     btkAppendEvent(acquisition,'Foot Strike', LeftHS(i,1), 'Left', '', 'Left Foot Strike', 1);
     btkAppendEvent(acquisition,'Foot Off', LeftTO(i,1), 'Left', '', 'Left Foot Off', 2);
     
end


end

