function [emgCaptured] = checkSessionEMG(subject, c3dinSession)
%checkEMG Summary of this function goes here
%   Input the subject and session name to determine if EMG exists for that
%   session. Function will return a 1 for yes and 0 for no. i indicates the
%   loop iteration, assuming the session are processed in chronological
%   order

% Dates change so be careful of that, check that the sessions are loaded in
% chronological order

load('Z:\s2921887\Google Drive\Load Sharing Main Data Collection\emgExist4Subject.mat')

% Find the subject number
Key   = 'Subject';
Index = strfind(subject, Key);
Value = sscanf(subject(Index(1) + length(Key):end), '%g', 1);

% Add onto the subject
subjectName = ['Subject', num2str(Value)];

emgUsedInSession = emgSubjects.(subjectName).(c3dinSession);

if emgUsedInSession == 1
     emgCaptured = 1;
     
elseif emgUsedInSession == 0
     emgCaptured = 0;
     
else
     disp('Not able to detect if EMG is in this session')
end

end
