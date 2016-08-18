function [trialsList] = trialsListGeneration_LS(inputDir)
% Function to get the trials list from the input folder 

% This file is part of Batch OpenSim Processing Scripts (BOPS).
% Copyright (C) 2015 Alice Mantoan, Monica Reggiani
% 
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License.
%
% Alice Mantoan, Monica Reggiani
% <ali.mantoan@gmail.com>, <monica.reggiani@gmail.com>


%%

trials=dir(inputDir);
j=1;

for k=3:length(trials)
     
     if (trials(k).isdir==1 && strcmp(trials(k).name,'Figures')==0 && strcmp(trials(k).name,'maxemg')==0 && strcmp(trials(k).name,'EMGs')==0)
          if exist([inputDir, filesep, trials(k).name, filesep, trials(k).name, '_grf.mot'], 'file')
               trialsList{j}=trials(k).name;
               j=j+1;
          end
     end
end