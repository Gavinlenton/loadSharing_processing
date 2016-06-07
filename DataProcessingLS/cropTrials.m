function [rightHS, rightTO] = cropTrials(acqLS, c3dFile_name, data)
%Take c3d file as input and crop into multiple gait cycles
%Load c3d files and insert heel strike event. Then crop the
% trial based on the heel strike

% First find time index for all events
[rightHS, leftHS, rightTO, leftTO] = findHeelStrike(data);

% I want to create a cloned acquisition, then extract each gait cycle
% and create a new c3d file of that cropped acquisition

for ii = 1:length(rightHS)-1
     
     % Clone acquisition
     acq_newLS = btkCloneAcquisition(acqLS);
     
     % Insert new events into clone
     insertGaitEvents(acq_newLS, rightHS, leftHS, rightTO, leftTO)
     
     % Check if events were actually appended
     [times, labels, descriptions, ids] = btkGetEventsValues(acq_newLS);
     
     if isempty(times)
          
          uiwait(msgbox('Warning: Events do not exist'));
          
     else
          
          % Crop the new acquisition based on time between heel strikes
          numFrames = rightHS(ii+1,:) - rightHS(ii,:);
          btkCropAcquisition(acq_newLS, rightHS(ii), numFrames);

          %Write the new acquisition
          filename = [c3dFile_name(1:end-4), num2str(ii), '.c3d'];
          btkWriteAcquisition(acq_newLS, filename);
         
     end
     
end

end

