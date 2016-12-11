function  ensembleOut = combineEnsembles(ensembles)
%Combine two ensemble averages into a mean and SD
%   INPUT - structure containing multiple ensembles (e.g., data from a gait
%           cycle).
%   OUTPUT - structure containing the ensemble mean and SD

% Length of the two ensembles

loops = fieldnames(ensembles);

% Loop through the ensembles

for i = 1:length(loops)-1
     n1 = length(ensembles.(i));
     n2 = length(ensembles.(i+1));
     
     for p = 1:1:n1
          var1 = p - 1;
          var2 = p - 1;
          
          ensembleAverage(p, 1) = ((n1 * ensembles.(i)(p)) + (n2 * ensembles.(i+1)(p))) / (n1 + n2);
          ensembleSD(p,1) = (n1 + n2) (var1 * (n1 - 1) + var2 * (n2-1) + n1 * 
     end

end

