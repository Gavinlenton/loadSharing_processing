function [HSRight, HSLeft, TORight, TOLeft] = findHeelStrike(data)
%Identify frames in which the heel marker is at the maximum distance from
%the sacrum marker
%   Run through the marker data and identify the peaks in which the heel
%   marker is the furthest distance from the sacrum marker. 
%   NOTE - need to perform this function for both the left and right leg
%   separately. 

% Find heel-strike
maxDistanceR = (data.marker_data.Markers.RCAL - data.marker_data.Markers.SAC1);
maxDistanceL = (data.marker_data.Markers.LCAL - data.marker_data.Markers.SAC2);

% Find toe-off
toeOffR = data.marker_data.Markers.SAC1 - data.marker_data.Markers.RMT1;
toeOffL = data.marker_data.Markers.SAC2 - data.marker_data.Markers.LMT1;

% Use the findpeaks function to determine when the heel marker is furthest
% from the sacrum marker. Specify the minimum distance between peaks so the
% function does not compute multiple peaks close to eachother. Can also
% specify a minimum peak height. 

[pks1, HSRight] = findpeaks(maxDistanceR(:,2), 'MinPeakDistance', 80);
[pks2, HSLeft] = findpeaks(maxDistanceL(:,2), 'MinPeakDistance', 80);

[pks1, TORight] = findpeaks(toeOffR(:,2), 'MinPeakDistance', 80);
[pks2, TOLeft] = findpeaks(toeOffL(:,2), 'MinPeakDistance', 80);


% Make sure events correspond with frame correctly. 
firstFrame = data.marker_data.First_Frame;

HSRight = HSRight + firstFrame;
TORight = TORight + firstFrame;
HSLeft = HSLeft + firstFrame;
TOLeft = TOLeft + firstFrame;


end

