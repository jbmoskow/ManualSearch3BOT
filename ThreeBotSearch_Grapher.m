% Graphs the data for the 3BOT Search experiment.
% Programmed by J. Moskowitz in July 2019

%% Graph the results based on analyzed data

function [GroupData] = ThreeBotSearch_Grapher(SUBS,FeatureSUBS,SpatialSUBS,GroupStruct,graphs)

% init
res = [1000,600];
GroupData = [];
fontSize = 14;

% coloring for individual participants 
c_colors = [228,26,28
55,126,184
77,175,74
152,78,163
255,127,0
25,255,150
166,86,40
247,129,191
153,153,153
30,50,240
1,255,255
40,100,100];  
c_colors = c_colors ./ 255;

RGB_LIGHTGREEN = [0/255 255/255 0/255]; 
RGB_DARKGREEN = [0/255 100/255 0/255];

% get weights
binWeights = GroupStruct{1}.ObjectWeights;

% get condition type
cond = zeros(1,length(GroupStruct));
for s = 1:length(SUBS)
    if ismember(SUBS{s},{'MC1','CVB1','JS1','LZ1','MS1','TG1','ET1','KC1','IH1',...
     'TSJ1','HA1','MB1'})
        cond(s) = 1; % feature condition
    elseif ismember(SUBS{s},{'SP1','YZ1','ET2','AK1','HB1','YZ2','AL1','KF1','AN1','NL1',...
    'CO1','BI1'})
        cond(s) = 2; % spatial condition
    end
end

%% Plot the proportion of lifts for each object weight across studies
if ismember(graphs,1)
    
    GroupData = NaN(length(FeatureSUBS),2);
    
    count = 1;
    for s = FeatureSUBS % for each feature condition subject
        
        T = GroupStruct{s}.TrialTable;
        
        expTrials = T.trialNum > 3;
        
        vars = {'liftWeights'}; % for these columns only
        T2 = T{expTrials,vars}; % select just this part of the table
        
        % concatenate all cells and convert to double
        allLifts = cat(2,T2{:});
        
        % bin counts for each lift weight
        binCounts = histcounts(allLifts,length(binWeights));
        binProp = binCounts / sum(binCounts); % convert to prop
        
        GroupData(count,1) = binProp(1); % assign to group data
        
        count = count + 1;
    end
    
    count = 1;
    for s = SpatialSUBS % for each spatial condition subject
        
        T = GroupStruct{s}.TrialTable;
        
        expTrials = T.trialNum > 3;
        
        vars = {'liftWeights'}; % for these columns only
        T2 = T{expTrials,vars}; % select just this part of the table
        
        % concatenate all cells and convert to double
        allLifts = cat(2,T2{:});
        
        % bin counts for each lift weight
        binCounts = histcounts(allLifts,length(binWeights));
        binProp = binCounts / sum(binCounts); % convert to prop
        
        GroupData(count,2) = binProp(1); % assign to group data
        
        count = count + 1;
    end
    
    % plotting
    figure(1);
    f1 = gcf;
    f1.Position = [0 0 res(1) res(2)];
    hold on;
    
    boxData = iosr.statistics.boxPlot(GroupData,'showScatter',true,...
        'symbolMarker','*','sampleSize',false,...
        'scatterMarker','.','scatterColor',rgb('orange'),...
        'scatterSize',800,'outlierSize',60);
    
    % labels and axes
    title('Proportion of Lifts for Lighter Object in 3BOT Search')
    xlabel('Experimental Group')
    xticklabels({'Feature Condition','Spatial Condition'})
    ylabel('Proportion of Lighter Object Lifts');
    axis([0 3 0 1]);
    
    % if object was randomly selected
    equalPropValue = 1 / 2;
    line([0 3],[equalPropValue equalPropValue],'Color','k',...
        'LineStyle','--','LineWidth',2)
    %     legend([SUBS 'Group']);
    
    count = 1;
    for c = 1:2 % for each condition
            p = 1;
        for s = 1:length(GroupData) % for each participant
            hold on
            if boxData.statistics.outliers_IX(s,c) == 1 % if there's an outlier and it's this subject
                xloc = boxData.x(c);
                yloc = GroupData(s,c);
            else
                xloc = boxData.handles.scatters(c).XData(min(p,length(boxData.handles.scatters(c).XData)));
                yloc = boxData.handles.scatters(c).YData(min(p,length(boxData.handles.scatters(c).YData)));
                p = p + 1;           
            end
            
            % FEATURE CONDITION
            if ismember(SUBS{count},{'MC1','LZ1','MS1','JS1','CVB1','MB1'})
                s1 = scatter(xloc,yloc,...
                    100,'MarkerFaceColor',RGB_DARKGREEN);
            elseif ismember(SUBS{count},{'TG1','ET1','KC1','IH1','TSJ1','HA1'})
                s2 = scatter(xloc,yloc,...
                    100,'MarkerFaceColor',RGB_LIGHTGREEN);
                % SPATIAL CONDITION
            elseif ismember(SUBS{count},{'SP1','ET2','HB1','AL1','KF1','CO1'})
                s1 = scatter(xloc,yloc,...
                    100,'MarkerFaceColor',RGB_DARKGREEN);
            elseif ismember(SUBS{count},{'YZ1','AK1','YZ2','NL1','AN1','BI1'})
                s2 = scatter(xloc,yloc,...
                    100,'MarkerFaceColor',RGB_LIGHTGREEN);
            end
            
            count = count + 1;
        end
    end
    
    legend([s1 s2],'Dark / Right Objects Lighter','Bright / Left Objects Lighter');
    
    set(gca,'FontSize',fontSize,'FontName','Myriad Pro');
    
end

%% Plot the prop. of lighter objects lifted for first N lifts
if ismember(graphs,2)

    % init
    N = 10;
    GroupData = NaN(length(FeatureSUBS),N);
    GroupProp = NaN(length(FeatureSUBS),N);
    GroupPer = NaN(length(FeatureSUBS),N);
    GroupBinCounts = [];
    
    for s = FeatureSUBS % for each subject
        
        T = GroupStruct{s}.TrialTable;
        vars = {'liftWeights'}; % for these columns only
        T2 = T{:,vars}; % select just this part of the table
        
        weightForNLifts = NaN(length(T2),N); % init double
        propLightLifts = NaN(length(T2),N); 
        
        for i = 1:length(T2) % for each trial
            totalRemain = 16;
            lightRemain = 8;
            for j = 1:N % for each lift
                try
                    weightForNLifts(i,j) = T2{i}(j);
                    propLightLifts(i,j) = lightRemain / totalRemain;
                    totalRemain = totalRemain - 1;
                    % if this was a light lift deduct 1 from the remaining 
                    if weightForNLifts(i,j) == binWeights(1) 
                        lightRemain = lightRemain - 1;
                    end
                catch % if lifts don't go that far
                    break
                end
            end
        end
        
        % Calculate the amount of data for each lift number (excluding NaNs)
        weightForNLiftsLen = zeros(1,N);
        PerNotNaN = zeros(1,N);
        for i = 1:N
            notNaN = ~isnan(weightForNLifts(:,i));
            weightForNLiftsLen(i) = length(find(notNaN == 1));
            PerNotNaN(i) = length((find(notNaN == 1))) / ...
                length(weightForNLifts(:,i)) * 100;
        end
        
        % bin counts for each lift weight
        for i = 1:N
            binCounts(i,:) = histcounts(weightForNLifts(:,i),length(binWeights));
            binProp(i,:) = binCounts(i,:) / sum(binCounts(i,:)); % convert to prop
        end
        
        % Mean Proportion of Light Lifts Remaining
        meanPropLight = nanmean(propLightLifts);
        
        % Plotting
        figure(2);
        f1 = gcf;
        set(f1,'Position',[0,0,res(1),res(2)])
        hold on;
        
        plot(1:N,binProp(:,1),'-','Marker','.','Color',rgb('lightgray'),...
            'LineWidth',0.5,'MarkerSize',10);
        
        GroupData(s,:) = binProp(:,1);
        GroupProp(s,:) = meanPropLight;
        GroupPer(s,:) = PerNotNaN;
        
        GroupBinCounts = cat(1,GroupBinCounts,binCounts);
        
     
    end
    
    % calculate and plot group means and std. error
    meanProp = nanmean(GroupData); % mean
    meanPer = nanmean(GroupPer); % mean
    errProp = nanstd(GroupData) ./ sqrt(length(GroupStruct)); % std error
    
    meanPropLight = nanmean(GroupProp);
    errPropLight = nanstd(GroupProp) ./ sqrt(length(GroupStruct));
    
    
    errorbar(1:N,meanProp,errProp,'LineWidth',3,'Color','k');
    
    errorbar(1:N,meanPropLight,errPropLight,'Color','k','LineStyle','--','LineWidth',2)
    
    % Add mean percentages of trials that included each lift number
    for i = 1:N
        text(i-0.1,0.7,sprintf('%.0f%%',meanPer(i)),...
            'FontSize',15,'FontName','Myriad Pro');
    end
    
    % labels and axes
    title(sprintf('Prop. of Lighter Objects Lifted for first %d Lifts - Feature Cond.',N));
    xlabel('Lift Number')
    ylabel('Prop. of Light Object Lifts');
    xticks(1:N);
    axis([0 N+1 0 1]);
   
%     legend([SUBS 'Group']);
    
    set(gca,'FontSize',fontSize,'FontName','Myriad Pro');

end

%% Plot the proportion of lifts for each object location
if ismember(graphs,3)
   
    figNum = 3;
    for group = 1:2
        
        if group == 1
            SUBIdx = FeatureSUBS;
        else
            SUBIdx = SpatialSUBS;
        end
        
        for locationType = 1:2    % 1 = xlocation, 2 = ylocation
            
            GroupData = NaN(length(FeatureSUBS),4);
            GroupBinCount = [];
            
            for s = SUBIdx % for each subject
                
                T = GroupStruct{s}.TrialTable;
                vars = {'liftLocations'}; % for these columns only
                T2 = T{:,vars}; % select just this part of the table
                
                % concatenate all cells and convert to double
                allLifts = cat(1,T2{:});
                allLifts = allLifts(:,locationType); % just the x-locs
                
                % bin counts for each lift weight
                binCounts = histcounts(allLifts,4); % there should only be 4 locs
                binProp = binCounts / sum(binCounts); % convert to prop
                
                % plotting
                figure(figNum);
                f1 = gcf;
                f1.Position = [0 0 res(1) res(2)];
                hold on;
                
                xval = unique(allLifts);
                
                if ismember(GroupStruct{s}.Subject,{'SP1','ET2','HB1','AL1','KF1','CO1'}) && locationType == 1
                    plot(xval,binProp,'-','Marker','.','Color',rgb('dimgray'),...
                        'LineWidth',0.5,'MarkerSize',10);
                elseif ismember(GroupStruct{s}.Subject,{'YZ1','AK1','YZ2','NL1','AN1','BI1'}) && locationType == 1
                    plot(xval,binProp,'-','Marker','.','Color',rgb('lightgray'),...
                        'LineWidth',0.5,'MarkerSize',10);
                else
                    plot(xval,binProp,'-','Marker','.','Color',rgb('lightgray'),...
                        'LineWidth',0.5,'MarkerSize',10);
                end
                
                fitobject = fit(xval,binProp','poly1')
                
                slopes(s) = fitobject.p1;
                
                GroupData(s,:) = binProp; % assign to group data
                GroupBinCount = cat(1,GroupBinCount,binCounts);
                
            end
            
            % calculate and plot group mean and std. error
            if group == 2 % this is looking at spatial condition participants
                for s = 1:length(GroupStruct)
                    if ismember(GroupStruct{s}.Subject,{'SP1','ET2','HB1','AL1','KF1','CO1'})
                        spatialCond(s) = 1; % objects are light on the right
                    elseif ismember(GroupStruct{s}.Subject,{'YZ1','AK1','YZ2','NL1','AN1','BI1'})
                        spatialCond(s) = 2; % objects are light on the left
                    end
                end
            end
            
            if group == 2 && locationType == 1
                meanPropRightLight = nanmean(GroupData(spatialCond == 1,:));
                errPropRightLight = nanstd(GroupData(spatialCond == 1,:)) ./ sqrt(6);
                
                meanPropLeftLight = nanmean(GroupData(spatialCond == 2,:));
                errPropLeftLight = nanstd(GroupData(spatialCond == 2,:)) ./ sqrt(6);
                
                e1 = errorbar(xval,meanPropRightLight,errPropRightLight,'LineWidth',3,'Color',RGB_DARKGREEN);
                e2 = errorbar(xval,meanPropLeftLight,errPropLeftLight,'LineWidth',3,'Color',RGB_LIGHTGREEN);
                
                legend([e1 e2],'Right Objects Lighter','Left Objects Lighter');
                
                [h,p] = ttest(slopes(spatialCond == 1),0)
                [h,p] = ttest(slopes(spatialCond == 2),0)
            else
                
                meanProp = nanmean(GroupData); % mean
                errProp = nanstd(GroupData) ./ sqrt(length(GroupStruct)); % std error
                errorbar(xval,meanProp,errProp,'LineWidth',3,'Color','k');
                
                [h,p] = ttest(slopes,0)
            end
            
            % labels and axes
            if locationType == 1
                title('Proportion of Lifts for Each Object X-location')
                xlabel('Object X-Location (cm)')
            else
                title('Proportion of Lifts for Each Object Y-location')
                xlabel('Object Y-Location (cm)')
            end
            ylabel('Proportion of Total Lifts');
            axis([-9 9 0 0.6]);
            
            % if object was randomly selected
            equalPropValue = 1 / 4;
            line([-8 8],[equalPropValue equalPropValue],'Color','k',...
                'LineStyle','--','LineWidth',2)
            %     legend([SUBS 'Group']);
            
            set(gca,'FontSize',fontSize,'FontName','Myriad Pro');
            
            figNum = figNum + 1;
        end
    end

end

%% Show Lift Height as a function of trial for Both Feature and Spatial condition
if ismember(graphs,4)
   
    FeatureDataLight = NaN(length(FeatureSUBS),90);
    FeatureDataHeavy = NaN(length(FeatureSUBS),90);
    SpatialDataLight = NaN(length(FeatureSUBS),90);
    SpatialDataHeavy = NaN(length(FeatureSUBS),90);

    for s = FeatureSUBS % for each feature condition subject
        
        T = GroupStruct{s}.TrialTable;
        vars = {'liftWeights','robotPeakLiftHeight'}; % for these columns only
        T2 = T(:,vars); % select just this part of the table
        
        liftHeightHeavy = NaN(height(T2),16); % init double with max lift length
        liftHeightLight = NaN(height(T2),16); % init double with max lift length
        
        for i = 1:height(T2) % for each trial
            for j = 1:size(T2.liftWeights{i},2)
                if T2.liftWeights{i}(j) == binWeights(1) % lighter weight
                    liftHeightLight(i,j) = T2.robotPeakLiftHeight{i}(j);
                else
                    liftHeightHeavy(i,j) = T2.robotPeakLiftHeight{i}(j);
                end
            end
        end
        
        % Mean and std. error for each lift
        meanLiftHeightLight = nanmean(liftHeightLight,2); % mean 
        meanLiftHeightHeavy = nanmean(liftHeightHeavy,2); % mean        
              
        % Plotting
        figure(7);
        f1 = gcf;
        set(f1,'Position',[0,0,res(1),res(2)])
        hold on;
        
        FeatureDataLight(s,1:length(meanLiftHeightLight)) = meanLiftHeightLight';
        FeatureDataHeavy(s,1:length(meanLiftHeightHeavy)) = meanLiftHeightHeavy';
    
    end
    
    for s = SpatialSUBS % for each spatial condition subject
        
        T = GroupStruct{s}.TrialTable;
        vars = {'liftWeights','robotPeakLiftHeight'}; % for these columns only
        T2 = T(:,vars); % select just this part of the table
        
        liftHeightHeavy = NaN(height(T2),16); % init double with max lift length
        liftHeightLight = NaN(height(T2),16); % init double with max lift length
        
        for i = 1:height(T2) % for each trial
            for j = 1:size(T2.liftWeights{i},2)
                if T2.liftWeights{i}(j) == binWeights(1) % lighter weight
                    liftHeightLight(i,j) = T2.robotPeakLiftHeight{i}(j);
                else
                    liftHeightHeavy(i,j) = T2.robotPeakLiftHeight{i}(j);
                end
            end
        end
        
        % Mean and std. error for each lift
        meanLiftHeightLight = nanmean(liftHeightLight,2); % mean 
        meanLiftHeightHeavy = nanmean(liftHeightHeavy,2); % mean        
        
        SpatialDataLight(s-12,1:length(meanLiftHeightLight)) = meanLiftHeightLight';
        SpatialDataHeavy(s-12,1:length(meanLiftHeightHeavy)) = meanLiftHeightHeavy';
    
    end
    
    % calculate and plot mean and std. error for blocks of trials
    n = 10; % Number of elements in each block
    % average across block size n
    meanLightFeatureG = NaN(length(FeatureSUBS),n-1);
    meanHeavyFeatureG = NaN(length(FeatureSUBS),n-1);
    meanLightSpatialG = NaN(length(FeatureSUBS),n-1);
    meanHeavySpatialG = NaN(length(FeatureSUBS),n-1);
    count = 1;
    for i = 1:n:length(FeatureDataLight)
        meanLightFeatureG(:,count) = nanmean(FeatureDataLight(:,i:i+n-1),2);
        meanHeavyFeatureG(:,count) = nanmean(FeatureDataHeavy(:,i:i+n-1),2);
        meanLightSpatialG(:,count) = nanmean(SpatialDataLight(:,i:i+n-1),2);
        meanHeavySpatialG(:,count) = nanmean(SpatialDataHeavy(:,i:i+n-1),2);
        count = count + 1;
    end
   
    % ALL STANDARD ERRORS SHOULD BE ACROSS SUBJECTS, NOT BLOCKS/TRIALS
    
    meanLightFeature = nanmean(meanLightFeatureG);
    errLightFeature = nanstd(meanLightFeatureG) ./ sqrt(length(FeatureSUBS)); % std error
    meanHeavyFeature = nanmean(meanHeavyFeatureG); % mean
    errHeavyFeature = nanstd(meanHeavyFeatureG) ./ sqrt(length(FeatureSUBS)); % std error
    
    meanLightSpatial = nanmean(meanLightSpatialG); % mean
    errLightSpatial = nanstd(meanLightSpatialG) ./ sqrt(length(FeatureSUBS)); % std error
    meanHeavySpatial = nanmean(meanHeavySpatialG); % mean
    errHeavySpatial = nanstd(meanHeavySpatialG) ./ sqrt(length(FeatureSUBS)); % std error
    
   % calculate across trial means
    CmeanLightFeature = mean(nanmean(meanLightFeatureG,2)); % mean
    CerrLightFeature = nanstd(nanmean(meanLightFeatureG,2)) ./ sqrt(length(FeatureSUBS)); % std error
    CmeanHeavyFeature = mean(nanmean(meanHeavyFeatureG,2)); % mean
    CerrHeavyFeature = nanstd(nanmean(meanHeavyFeatureG,2)) ./ sqrt(length(FeatureSUBS)); % std error
    
    CmeanLightSpatial = mean(nanmean(meanLightSpatialG,2)); % mean
    CerrLightSpatial = nanstd(nanmean(meanLightSpatialG,2)) ./ sqrt(length(FeatureSUBS)); % std error
    CmeanHeavySpatial = mean(nanmean(meanHeavySpatialG,2)); % mean
    CerrHeavySpatial = nanstd(nanmean(meanHeavySpatialG,2)) ./ sqrt(length(FeatureSUBS)); % std error
    
    subplot(1,2,1)
    hold on;
    e1 = errorbar(1:9,meanLightFeature,errLightFeature,'LineWidth',2,'Color',rgb('darkgray'));
    e2 = errorbar(1:9,meanHeavyFeature,errHeavyFeature,'LineWidth',2,'Color',rgb('dimgray'));
    e3 = errorbar(1:9,meanLightSpatial,errLightSpatial,'LineWidth',2,'Color',rgb('lightskyblue'));
    e4 = errorbar(1:9,meanHeavySpatial,errHeavySpatial,'LineWidth',2,'Color',rgb('dodgerblue'));

    % labels and axes
    title('Avg. Lift Height by Condition and Weight Across Trials')
    xlabel('Block Number')
    ylabel('Lift Height (cm)');
    axis([0 10 3 12]);

    legend([e1 e2 e3 e4],'Light Lifts Feature Condition','Heavy Lifts Feature Condition',...
        'Light Lifts Spatial Condition','Heavy Lifts Spatial Condition')

    set(gca,'FontSize',fontSize,'FontName','Myriad Pro');

    % mean of time series
    subplot(1,2,2); 
    hold on;
    e5 = errorbar(1:2,[CmeanLightFeature CmeanHeavyFeature],[CerrLightFeature CerrHeavyFeature],'LineWidth',2,'Color',rgb('black'));
    e6 = errorbar(1:2,[CmeanLightSpatial CmeanHeavySpatial],[CerrLightSpatial CerrHeavySpatial],'LineWidth',2,'Color',rgb('blue'));
    axis([0.5 2.5 3 12]);
    title('Avg. Lift Height by Condition and Weight')
    xlabel('Lift Type')
    xticks([1 2]);
    xticklabels({'Light Lift','Heavy Lift'});
    ylabel('Lift Height (cm)');
    
    legend([e5 e6],'Feature Condition','Spatial Condition');
    
    set(gca,'FontSize',fontSize,'FontName','Myriad Pro');
    
end

%% mean search times
if ismember(graphs,5)
    
    GroupData = zeros(length(GroupStruct),1);
    
    for s = 1:length(GroupStruct)
        
        T = GroupStruct{s}.TrialTable;
        
        expTrials = T.trialNum > 3;
        
        vars = {'searchTime'}; % for these columns only
        T2 = T{expTrials,vars}; % select just this part of the table
        
        GroupData(s,:) = mean(T2);
    
    end
    
    
end

%% mean lift count
if ismember(graphs,6)
    
    GroupData = zeros(length(GroupStruct),1);
    
    for s = 1:length(GroupStruct)
        
        T = GroupStruct{s}.TrialTable;
        
        expTrials = T.trialNum > 3;
        
        vars = {'liftWeights'}; % for these columns only
        T2 = T{expTrials,vars}; % select just this part of the table
        
        numLifts = zeros(length(T2),1);
        
        for i = 1:length(T2)
            numLifts(i) = length(T2{i});
        end
        
        GroupData(s,:) = mean(numLifts);
    
    end
    
    
end


