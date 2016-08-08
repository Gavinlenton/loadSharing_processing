function [rightHS, rightTO] = cropTrials(acqLS, c3dFile_name, data)
%Take c3d file as input and crop into multiple gait cycles
%Load c3d files and insert heel strike event. Then crop the
% trial based on the heel strike

% First find time index for all events
[rightHS, rightTO] = findHeelStrike(data);

% Create new c3d files of gait cycles

for ii = 1:length(rightHS)-1
     
	% Clone acquisition
	acq_newLS = btkCloneAcquisition(acqLS);
	
	% Crop the new acquisition based on time between heel strikes
	numFrames = rightHS(ii+1,:) - rightHS(ii,:);
	btkCropAcquisition(acq_newLS, rightHS(ii), numFrames);
	
	%Write the new acquisition
	filename = [c3dFile_name(1:end-4), num2str(ii), '.c3d'];
	btkWriteAcquisition(acq_newLS, filename);
         
end

end

