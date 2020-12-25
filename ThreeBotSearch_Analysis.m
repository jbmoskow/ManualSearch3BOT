% Analyzes the data for the 3BOT Search Experiment
% Programmed by J. Moskowitz in July 2019

%% Returns the contents of the analysis in a struct
function [AnalysisStruct] = ThreeBotSearch_Analysis(subjectID,rawDir)

% load raw data structures
cd(rawDir)
fileinfo = dir(['*' subjectID '.mat']); % find correct filename
data = load(fileinfo.name);

% init 
dataLength = length(data.ObjectInfo);
numTrials = size(data.Frames,1);
liftWeights = cell(numTrials,1); % weight of objects lifted
liftLocations = cell(numTrials,1); % Nx2 matrix of the x and y-pos of objects lifted
liftTimes = cell(numTrials,1); % Nx2 matrix of lift onset and offset (ms)
robotMeanForce = cell(numTrials,1); % mean load force (FZ) on the hand during lifts (N);
robotPeakDeltaForce = cell(numTrials,1); % Peak rate of change of force during lift (N/s);
robotPeakLiftHeight = cell(numTrials,1); % Peak height of the robot during lift (cm);
trialNum = zeros(numTrials,1); % trial number from raw data
blockNum = zeros(numTrials,1); % block number
searchTime = zeros(numTrials,1); % search time (s)

%% Grab non-trial specific data
subjsplit = split(fileinfo.name,'_');
expDate = [subjsplit{2} '_' subjsplit{3}];

AnalysisStruct.Subject = subjectID;
AnalysisStruct.Date = expDate;
   
% object weights
uniqueWeights = unique(data.ObjectInfo(:,1,:));
% excluding zero, and rounded
AnalysisStruct.ObjectWeights = round(uniqueWeights(2:end),3,'significant'); 

count = 1;
%% Loop through each trial
for t = 1:numTrials
    
    disp(t);
    
    % Trial Types are 0 = rest, 1 = done, 2 = search
    currTrialType = data.TrialData.TrialType(t);
    
    if currTrialType == 2 % non-practice, non-bias
        
        liftData = reshape(data.ObjectInfo(t,1,:),[1 dataLength]);
        objectX = reshape(data.ObjectInfo(t,2,:),[1 dataLength]);
        objectY = reshape(data.ObjectInfo(t,3,:),[1 dataLength]);
        targetData = reshape(data.ObjectInfo(t,4,:),[1 dataLength]);
        
        [~,idxLifts] = findLiftOnsets(liftData,1);
        [~,idxTarget] = findLiftOnsets(targetData,1);
        
        if ~isnan(idxLifts) % there was some actual lifts
            % find where the actual target was among the objects lifted
            % in case an additional object was lifted post-target or target was
            % not found
            lastTargSearched = find(idxLifts(:,1) == idxTarget(1,1));
        else
            lastTargSearched = NaN;
        end
        
        % the target was actually found
        if ~isnan(lastTargSearched) 
            %% Store the trial type and block num
            
            trialNum(count) = t;
            blockNum(count) = data.TrialData.block_count(t);
            
            %% Store object weights lifted for this trial
            
            for i = 1:size(idxLifts,1) % for each lift
                liftWeights{count}(i) = liftData(idxLifts(i,1));
            end
            
            % update these matrices to reflect correct values
            liftWeights{count} = liftWeights{count}(1:lastTargSearched);
            
            % round weight values
            liftWeights{count} = round(liftWeights{count},3,'significant');
            
            %% Store location of object lifted

            % values obtained from line 238-244 in cfg_SearchObjects.m
            % removes offset so all trials have common lift locations
            baseTargLocs = [-8,-8;-8,-2.66666666666667;-8,2.66666666666667;-8,8;-2.66666666666667,-8;-2.66666666666667,-2.66666666666667;-2.66666666666667,2.66666666666667;-2.66666666666667,8;2.66666666666667,-8;2.66666666666667,-2.66666666666667;2.66666666666667,2.66666666666667;2.66666666666667,8;8,-8;8,-2.66666666666667;8,2.66666666666667;8,8];
            
            % use idx from unique lifts to find position upon lift init
            trialTargLocsX = objectX(idxLifts(1:lastTargSearched,1));
            trialTargLocsY = objectY(idxLifts(1:lastTargSearched,1));
            
            targAt = zeros(length(trialTargLocsX),1);
            for i = 1:length(trialTargLocsX) % for each lift
                % calculate distance of this position to each target
                distToTargets = sqrt((baseTargLocs(:,1)-trialTargLocsX(i)).^2 + ...
                    (baseTargLocs(:,2)-trialTargLocsY(i)).^2);
                [~,targAt(i)] = min(distToTargets);
            end
            
            liftLoc = [baseTargLocs(targAt,1) baseTargLocs(targAt,2)];
            
            % ASK MARTIN WHERE THE REST POSITION IS LOCATED
            liftLocations{count} =  liftLoc;
            
            %% Store time point of when the lift was initiated and released
            
            liftTimes{count} = idxLifts(1:lastTargSearched,:);
            
            %% Store mean Load Force (Robot FZ) during lift
            
            peakLoads = zeros(lastTargSearched,1);
            peakChangeLoads = zeros(lastTargSearched,1);
            peakHeights = zeros(lastTargSearched,1);
            for i = 1:lastTargSearched
                % grab robot force Z (load) during lift
                RFZ = data.RobotForces(t,3,idxLifts(i,1):idxLifts(i,2));
                
                % grab robot position z during lift
                RPZ = data.RobotPosition(t,3,idxLifts(i,1):idxLifts(i,2));
                
                % calculate peak of this load
                peakLoads(i) = abs(min(RFZ));
                        
                % filter it using a 4th order Butterworth filter, with a
                % 14 Hz low-pass filter
                fs = 1000;
                [num,den] = butter(4,2*14/fs,'low');
                butter_filtered = filtfilt(num,den,RPZ(:));
                
                % calculate peak rate of change of this Force 
                % (converted to N/s)
                peakChangeLoads(i) = abs(min(diff(butter_filtered)));
                
                % Calculate peak lift height
                peakHeights(i) = max(RPZ) - min(RPZ);
                
%                 figure(99)
%                 plot(diff(butter_filtered));
%                 hold on;
                
            end
            
            robotMeanForce{count} = peakLoads;
            robotPeakDeltaForce{count} = peakChangeLoads;
            robotPeakLiftHeight{count} = peakHeights;
            
            %% Store time at which target was found
            
            timeToTarget = idxLifts(end,1) / 1000;
            
            searchTime(count) = timeToTarget;
            
            
            %% increment
            count = count + 1;
            
            %% exemplar
            RPX = data.RobotPosition(t,1,:);
            RPY = data.RobotPosition(t,2,:);
            RPZ = data.RobotPosition(t,3,:);
            ThreeBotSearch_Exemplar(subjectID,t,RPX(:),RPY(:),RPZ(:),idxLifts,data.TrialData(t,:))
            
            keyboard

        end
    end
end

AnalysisStruct.TrialTable = table(blockNum,trialNum,searchTime,...
    liftWeights,liftLocations,liftTimes,robotMeanForce,...
    robotPeakDeltaForce,robotPeakLiftHeight);

% cut off empty trials
rows = AnalysisStruct.TrialTable.trialNum ~= 0;
AnalysisStruct.TrialTable = AnalysisStruct.TrialTable(rows,:);

end