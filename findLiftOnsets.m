function [iStart1,iStartEnd] = findLiftOnsets(conditions,minDuration)
% [iStart1,iStartEnd] = findConditionsAreTrue(conditions,minDuration)
%
% Find start and end indices of intervals during which certain conditions are true.
% Example: find when x and y positions of hand is in a target. 
% Outputs are NaN if interval is not found.

% Anouk de Brouwer
% June 2016

% make sure the input is in columns
if size(conditions,1)<size(conditions,2)
    conditions = conditions';
end

% find start and end indices of when all conditions are fulfilled
allConditionsTrue = prod(conditions,2); % logical AND in case there are multiple columns
% get column of indx for all rising edges, when condition is on
startTrue = find(diff(allConditionsTrue)>0); 
% get column of indx for all falling edges, when condition turns off
endTrue = find(diff(allConditionsTrue)<0);

% make sure the number of start and end indices is the same
% if condition is true from start of search window
if allConditionsTrue(1) && (isempty(startTrue) || (startTrue(1)>endTrue(1)))
    startTrue = [0; startTrue]; 
end
% if condition is true until end of search window
if allConditionsTrue(end) && (isempty(endTrue) || (startTrue(end)>endTrue(end)))
    endTrue = [endTrue; length(allConditionsTrue)];
end

% check duration of intervals during which all conditions are true
trueStartEnd = [startTrue endTrue];
if ~isempty(trueStartEnd)
    durationTrue = trueStartEnd(:,2) - trueStartEnd(:,1);
    trueStartEnd = trueStartEnd(durationTrue>=minDuration,:);
    if ~isempty(trueStartEnd)
        % increment start idx by 1 to reflect actual location
        iStart1 = trueStartEnd(1)+1;
        iStartEnd = [trueStartEnd(:,1)+1 trueStartEnd(:,2)];
    else
        iStart1 = NaN;
        iStartEnd = NaN;
    end
else
    iStart1 = NaN;
    iStartEnd = NaN;
end
