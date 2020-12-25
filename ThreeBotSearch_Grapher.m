% Graphs the data for the 3BOT Search experiment.
% Programmed by J. Moskowitz in July 2019

%% Graph the results based on analyzed data

function [GroupData] = ThreeBotSearch_Grapher(SUBS,GroupStruct,graphs,NEWSUBS,OLDSUBS)

% init
res = [400,400];
GroupData = [];
fontSize = 14;

% if ~ismember(graphs,11)
%     GroupStruct = GroupStruct(1:length(NEWSUBS));
%     SUBS = SUBS(1:length(NEWSUBS));
% end

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
    

%% Plot the proportion of lifts for each object weight
if ismember(graphs,1)
    
    GroupData = NaN(length(GroupStruct),length(binWeights));
    
    for s = 1:length(GroupStruct) % for each subject
        
        T = GroupStruct{s}.TrialTable;
        
        expTrials = T.trialNum > 3;
        
        vars = {'liftWeights'}; % for these columns only
        T2 = T{expTrials,vars}; % select just this part of the table
        
        % concatenate all cells and convert to double
        allLifts = cat(2,T2{:});
        
        % bin counts for each lift weight
        binCounts(s,:) = histcounts(allLifts,length(binWeights));
        binProp = binCounts(s,:) / sum(binCounts(s,:)); % convert to prop
        
        % plotting
%         figure(1);
%         f1 = gcf;
%         f1.Position = [0 0 res(1) res(2)];
%         hold on;
        
%         plot(binWeights,binProp,'-','Marker','.','Color',c_colors(s,:),...
%             'LineWidth',0.5,'MarkerSize',10);
        
        GroupData(s,:) = binProp; % assign to group data
        
    end
    
    % plotting
    figure(1);
    f1 = gcf;
    f1.Position = [0 0 res(1)/2 res(2)];
    ax = gca;
    hold on;
    
    % calculate and plot group mean and std. error
    meanProp = nanmean(GroupData); % mean
    errProp = nanstd(GroupData,1) ./ sqrt(length(GroupStruct)); % std error
    
    boxData = iosr.statistics.boxPlot(binWeights(1),GroupData(:,1),'showScatter',true,...
    'symbolMarker','*','sampleSize',true,...
        'scatterMarker','.','scatterColor',rgb('orange'),...
        'scatterSize',800,'outlierSize',60);  
    
    % labels and axes
    title('Proportion of Lifts for Lighter Object')
    ylabel('Proportion of Light Object Lifts');
    axis([-.8 1 0.4 1]);
        
    % if object was randomly selected
    equalPropValue = 1 / length(binWeights);
    line([-2 2],[equalPropValue equalPropValue],'Color','k',...
        'LineStyle','--','LineWidth',2)
%     legend([SUBS 'Group']);
    
    counter = 1;
    for s = 1:length(GroupData)
        hold on
        if boxData.statistics.outliers_IX(s) == 1 % if there's an outlier and it's this subject
            xloc = boxData.x;
            yloc = GroupData(s,1);
        else
            xloc = boxData.handles.scatters.XData(counter);
            yloc = boxData.handles.scatters.YData(counter);
            counter = counter + 1;
        end
            
%         text(xloc+0.01,yloc+0.01,SUBS{s},'FontSize',12) % label P name
        % FEATURE CONDITION
        if ismember(SUBS{s},{'MC1','LZ1','MS1','JS1','CVB1','MB1'})
            s1 = scatter(xloc,yloc,...
                100,'MarkerFaceColor',RGB_DARKGREEN);
        elseif ismember(SUBS{s},{'TG1','ET1','KC1','IH1','TSJ1','HA1'})
            s2 = scatter(xloc,yloc,...
                100,'MarkerFaceColor',RGB_LIGHTGREEN);
        % SPATIAL CONDITION
        elseif ismember(SUBS{s},{'SP1','ET2','HB1','AL1','KF1','CO1'})
            s1 = scatter(xloc,yloc,...
                100,'MarkerFaceColor',RGB_DARKGREEN);
        elseif ismember(SUBS{s},{'YZ1','AK1','YZ2','NL1','AN1','BI1'})
            s2 = scatter(xloc,yloc,...
                100,'MarkerFaceColor',RGB_LIGHTGREEN);
        end
    end
    
    % this is a feature condition subject
    if ismember(SUBS{1},{'MC1','CVB1','JS1','LZ1','MS1','TG1','ET1','KC1','IH1',...
            'TSJ1','HA1','MB1'})
        legend([s1 s2],'Dark Objects Lighter','Bright Objects Lighter');
    % this is a spatial condition subject
    elseif ismember(SUBS{1},{'SP1','YZ1','ET2','AK1','HB1','YZ2','AL1','KF1','AN1','NL1',...
            'CO1','BI1'})
        legend([s1 s2],'Right Objects Lighter','Left Objects Lighter'); 
    end
     
    set(gca,'FontSize',fontSize,'FontName','Myriad Pro');

end

%% Plot the proportion of lifts for each object weight across studies
if ismember(graphs,11)
    
    GroupData = NaN(length(NEWSUBS),length(binWeights));
    
    for s = 1:length(NEWSUBS) % for each new subject
        
        T = GroupStruct{s}.TrialTable;
        
        expTrials = T.trialNum > 3;
        
        vars = {'liftWeights'}; % for these columns only
        T2 = T{expTrials,vars}; % select just this part of the table
        
        % concatenate all cells and convert to double
        allLifts = cat(2,T2{:});
        
        % bin counts for each lift weight
        binCounts = histcounts(allLifts,length(binWeights));
        binProp = binCounts / sum(binCounts); % convert to prop
        
        % plotting
%         figure(1);
%         f1 = gcf;
%         f1.Position = [0 0 res(1) res(2)];
%         hold on;
        
%         plot(binWeights,binProp,'-','Marker','.','Color',c_colors(s,:),...
%             'LineWidth',0.5,'MarkerSize',10);
        
        GroupData(s,1) = binProp(1); % assign to group data
        
    end
    
    % get new weights
    binWeights = GroupStruct{length(NEWSUBS)+1}.ObjectWeights;
    
    for s = 1:length(OLDSUBS) % for each old subject
      
        i = length(NEWSUBS) + s;

        T = GroupStruct{i}.TrialTable;
        
        vars = {'liftWeights'}; % for these columns only
        T2 = T{:,vars}; % select just this part of the table
        
        % concatenate all cells and convert to double
        allLifts = cat(2,T2{:});
        
        % bin counts for each lift weight
        binCounts = histcounts(allLifts,length(binWeights));
        binProp = binCounts / sum(binCounts); % convert to prop
        
        % average across the first half lightest objects
        binProp = sum(binProp(1:(length(binProp)/2)));
        
        % plotting
        %         figure(1);
        %         f1 = gcf;
        %         f1.Position = [0 0 res(1) res(2)];
        %         hold on;
        
        %         plot(binWeights,binProp,'-','Marker','.','Color',c_colors(s,:),...
        %             'LineWidth',0.5,'MarkerSize',10);
        
        GroupData(s,2) = binProp; % assign to group data
        
    end
    
    % plotting
    figure(1);
    f1 = gcf;
    f1.Position = [0 0 res(1) res(2)];
    hold on;
    
    % calculate and plot group mean and std. error
    meanProp = nanmean(GroupData); % mean
    errProp = nanstd(GroupData,1) ./ sqrt(length(GroupStruct)); % std error
    
    boxData = iosr.statistics.boxPlot(GroupData,'showScatter',true,...
    'symbolMarker','*','sampleSize',true,...
        'scatterMarker','.','scatterColor',rgb('orange'),...
        'scatterSize',800,'outlierSize',60);  
    
    % labels and axes
    title('Proportion of Lifts for Lighter Object in 3BOT Search')
    xlabel('Experimental Group')
    xticklabels({'Jennifer \newlineFeature Condition','Summer 2019'})
    ylabel('Proportion of Total Lifts');
    axis([0 3 0 1]);
        
    % if object was randomly selected
    equalPropValue = 1 / 2;
    line([0 3],[equalPropValue equalPropValue],'Color','k',...
        'LineStyle','--','LineWidth',2)
%     legend([SUBS 'Group']);
    
    for s = 1:length(NEWSUBS)
        text(boxData.handles.scatters(1).XData(s)+0.01,...
            boxData.handles.scatters(1).YData(s)+0.01,SUBS{s},'FontSize',12)
        if ismember(SUBS{s},{'MC1','LZ1','MS1','JS1','CVB1','MB1'})
            hold on
            s1 = scatter(boxData.handles.scatters(1).XData(s),...
                boxData.handles.scatters(1).YData(s),...
                100,'MarkerFaceColor',RGB_DARKGREEN);
        elseif ismember(SUBS{s},{'TG1','ET1','KC1','IH1','TSJ1','HA1'})
            s2 = scatter(boxData.handles.scatters(1).XData(s),...
                boxData.handles.scatters(1).YData(s),...
                100,'MarkerFaceColor',RGB_LIGHTGREEN);
        end
    end
    
    legend([s1 s2],'Dark Objects Lighter','Bright Objects Lighter');
 
    set(gca,'FontSize',fontSize,'FontName','Myriad Pro');

end

%% Plot the prop. of lighter objects lifted for first N lifts
if ismember(graphs,2)

    % init
    N = 10;
    GroupData = NaN(length(GroupStruct),N);
    GroupProp = NaN(length(GroupStruct),N);
    GroupPer = NaN(length(GroupStruct),N);
    GroupBinCounts = [];
    
    for s = 1:length(GroupStruct) % for each subject
        
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
        
%         % Mean and std. error for each lift
%         meanNLifts = nanmean(weightForNLifts); % mean
%         errNLifts = nanstd(weightForNLifts) ./ sqrt(weightForNLiftsLen); % std error
        
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
    title(sprintf('Prop. of Lighter Objects Lifted for first %d Lifts',N));
    xlabel('Lift Number')
    ylabel('Prop. of Light Object Lifts');
    xticks(1:N);
    axis([0 N+1 0 1]);
%     legend(SUBS);
   
     % if object was randomly selected
%     equalPropValue = 1 / length(binWeights);
%     line([0 10],[equalPropValue equalPropValue],'Color','k',...
%         'LineStyle','--','LineWidth',2)

    
    
%     legend([SUBS 'Group']);
    
    set(gca,'FontSize',fontSize,'FontName','Myriad Pro');

end

%% Plot the proportion of lifts for each object distance
if ismember(graphs,3)
    
    GroupData = NaN(length(GroupStruct),4);
    GroupBinCount = [];
    
    for s = 1:length(GroupStruct) % for each subject
        
        T = GroupStruct{s}.TrialTable;
        vars = {'liftLocations'}; % for these columns only
        T2 = T{:,vars}; % select just this part of the table
        
        % concatenate all cells and convert to double
        allLifts = cat(1,T2{:});
        allLifts = allLifts(:,2); % just the y-locs
        
        % bin counts for each lift weight
        binCounts = histcounts(allLifts,4); % there should only be 4 locs
        binProp = binCounts / sum(binCounts); % convert to prop
        
        % plotting
        figure(3);
        f1 = gcf;
        f1.Position = [0 0 res(1) res(2)];
%         subplot(1,2,1);
        hold on;
        
        xval = unique(allLifts);
        
        plot(xval,binProp,'-','Marker','.','Color',rgb('lightgray'),...
            'LineWidth',0.5,'MarkerSize',10);
        
        fitobject = fit(xval,binProp','poly1')
        
        slopes(s) = fitobject.p1;
        
        GroupData(s,:) = binProp; % assign to group data
        GroupBinCount = cat(1,GroupBinCount,binCounts);
       
    end
    
    % calculate and plot group mean and std. error
    meanProp = nanmean(GroupData); % mean
    errProp = nanstd(GroupData) ./ sqrt(length(GroupStruct)); % std error
    
    errorbar(unique(allLifts),meanProp,errProp,'LineWidth',3,'Color','k');
    
    [h,p] = ttest(slopes,0)
    
    % labels and axes
    title('Proportion of Lifts for Each Object Y-location')
    xlabel('Object Y-Location (cm)')
    ylabel('Proportion of Total Lifts');
    axis([-9 9 0 0.6]);
    
        % if object was randomly selected
    equalPropValue = 1 / 4;
    line([-8 8],[equalPropValue equalPropValue],'Color','k',...
        'LineStyle','--','LineWidth',2)
    %     legend([SUBS 'Group']);
    
%     subplot(1,2,2);
%     hold on;
%     
%     boxData = iosr.statistics.boxPlot(slopes','showScatter',true,...
%     'symbolMarker','*','sampleSize',true,...
%         'scatterMarker','.','scatterColor',rgb('orange'),...
%         'scatterSize',800,'outlierSize',60);  
%     line([0.25 1.75], [0 0],'Color','k',...
%          'LineStyle','--','LineWidth',2);
%     axis([0.25 1.75 -0.05 .05])
   
    set(gca,'FontSize',fontSize,'FontName','Myriad Pro');

end

%% Plot the proportion of lifts for each object location
if ismember(graphs,3.1)
    
    locationType = 1; % 1 = xlocation, 2 = ylocation
    
    GroupData = NaN(length(GroupStruct),4);
    GroupBinCount = [];
    
    for s = 1:length(GroupStruct) % for each subject
        
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
        figure(3);
        f1 = gcf;
        f1.Position = [0 0 res(1) res(2)];
%         subplot(1,2,1)
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
    if cond(1) == 2 % this is looking at spatial condition participants
        for s = 1:length(GroupStruct)
            if ismember(GroupStruct{s}.Subject,{'SP1','ET2','HB1','AL1','KF1','CO1'})
                spatialCond(s) = 1; % objects are light on the right
            elseif ismember(GroupStruct{s}.Subject,{'YZ1','AK1','YZ2','NL1','AN1','BI1'})
                spatialCond(s) = 2; % objects are light on the left
            end
        end
    end
        
    if cond(1) == 2 && locationType == 1
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
    
%     subplot(1,2,2)
%     hold on;
    
%     if cond(1) == 2 && locationType == 1 % if this is the spatial condition
%         
%         boxData = iosr.statistics.boxPlot(1:2,[slopes(spatialCond == 1)' slopes(spatialCond == 2)'],'showScatter',true,...
%         'symbolMarker','*','sampleSize',true,...
%             'scatterMarker','.','scatterColor',rgb('orange'),...
%             'scatterSize',800,'outlierSize',60);  
%         line([0.25 2.75],[0 0],'Color','k',...
%             'LineStyle','--','LineWidth',2);
%         axis([0.25 2.75 -.05 .05])
%         
%     else
%         
%         boxData = iosr.statistics.boxPlot(slopes','showScatter',true,...
%             'symbolMarker','*','sampleSize',true,...
%             'scatterMarker','.','scatterColor',rgb('orange'),...
%             'scatterSize',800,'outlierSize',60);
%         line([0.25 1.75],[0 0],'Color','k',...
%             'LineStyle','--','LineWidth',2);
%         axis([0.25 1.75 -.05 .05])
%     end

   
    set(gca,'FontSize',fontSize,'FontName','Myriad Pro');

end

%% Plot the avg. X and Y location of each object lifted for first N lifts
if ismember(graphs,3.2)

    % init
    N = 10;
    GroupXLoc = NaN(length(GroupStruct),N);
    GroupYLoc = NaN(length(GroupStruct),N);
    GroupPer = NaN(length(GroupStruct),N);
    spatialCond = NaN(length(GroupStruct),1);
    for s = 1:length(GroupStruct) % for each subject
        
        T = GroupStruct{s}.TrialTable;
        vars = {'liftLocations'}; % for these columns only
        T2 = T{:,vars}; % select just this part of the table
        
        xLocForNLifts = NaN(length(T2),N); % init double
        yLocForNLifts = NaN(length(T2),N);
        
        for i = 1:length(T2) % for each trial
            for j = 1:N % for each lift
                try
                    xLocForNLifts(i,j) = T2{i}(j,1);
                    yLocForNLifts(i,j) = T2{i}(j,2);
                catch % if lifts don't go that far
                    break
                end
            end
        end
        
        % Calculate the amount of data for each lift number (excluding NaNs)
        xLocForNLiftsLen = zeros(1,N);
        for i = 1:N
            notNaN = ~isnan(xLocForNLifts(:,i));
            xLocForNLiftsLen(i) = length(find(notNaN == 1));
        end
        
        % Mean and std. error for each lift
        meanxLoc = nanmean(xLocForNLifts); % mean
        errxLoc = nanstd(xLocForNLifts) ./ sqrt(xLocForNLiftsLen); % std error
        
        meanyLoc = nanmean(yLocForNLifts); % mean
        erryLoc = nanstd(yLocForNLifts) ./ sqrt(xLocForNLiftsLen); % std error
        
        % Plotting
        figure(4);
        f1 = gcf;
        set(f1,'Position',[0,0,res(1),res(2)])
        hold on;
        
        yyaxis left
        plot(1:N,meanyLoc,'-','Marker','.','Color',rgb('lightgray'),...
            'LineWidth',0.5,'MarkerSize',10);
        
        yyaxis right
        if ismember(GroupStruct{s}.Subject,{'SP1','ET2','HB1','AL1','KF1','CO1'})
            p2 = plot(1:N,meanxLoc,'-','Marker','.','Color',rgb('lightgray'),...
                'LineWidth',0.5,'MarkerSize',10);
            spatialCond(s) = 1;
        elseif ismember(GroupStruct{s}.Subject,{'YZ1','AK1','YZ2','NL1','AN1','BI1'})
            p3 = plot(1:N,meanxLoc,'-','Marker','.','Color',rgb('lightgray'),...
            'LineWidth',0.5,'MarkerSize',10);
            spatialCond(s) = 2;
        else
            p2 = plot(1:N,meanxLoc,'--','Marker','.','Color',rgb('lightgray'),...
                'LineWidth',0.5,'MarkerSize',10);
        end
        
        GroupXLoc(s,:) = meanxLoc;
        GroupYLoc(s,:) = meanyLoc;
     
    end
    
    % calculate and plot group means and std. error
    if cond(1) == 2 % spatial condition
        meanX1 = nanmean(GroupXLoc(spatialCond == 1,:));
        errX1 = nanstd(GroupXLoc(spatialCond == 1,:)) ./ sqrt(6); % std error
        meanX2 = nanmean(GroupXLoc(spatialCond == 2,:));
        errX2 = nanstd(GroupXLoc(spatialCond == 2,:)) ./ sqrt(6); % std error
    else
        meanX = nanmean(GroupXLoc);
        errX = nanstd(GroupXLoc) ./ sqrt(length(GroupStruct)); % std error
    end
   
    meanY = nanmean(GroupYLoc); % mean
    errY = nanstd(GroupYLoc) ./ sqrt(length(GroupStruct)); % std error
    
    yyaxis left
    e1 = errorbar(1:N,meanY,errY,'LineWidth',3,'Color',rgb('black'));
    ylabel('Avg. Y Location (cm)');
    
    yyaxis right  
    if cond(1) == 2
        e2 = errorbar(1:N,meanX1,errX1,'LineWidth',3,'Color',RGB_DARKGREEN);
        e3 = errorbar(1:N,meanX2,errX2,'LineWidth',3,'LineStyle','-','Color',RGB_LIGHTGREEN);
    else
        e2 = errorbar(1:N,meanX2,errX2,'LineWidth',3,'Color',rgb('cyan'));
    end
    ylabel('Avg. X Location (cm)');
    
    % labels and axes
    title(sprintf('Avg. X and Y loc for first %d Lifts in 3BOT Search',N));
    xlabel('Lift Number')
    xticks(1:N);
    axis([0 N+1 -10 10]);
   
     % if object was randomly selecte
    line([0 10],[0 0],'Color','k',...
        'LineStyle','--','LineWidth',2)
    
    % legend([SUBS 'Group']);
    if ismember(GroupStruct{1}.Subject,{'SP1','YZ1','ET2','AK1','HB1','YZ2','AL1','KF1','AN1','NL1',...
    'CO1','BI1'})
        legend([e1 e2 e3],'Y-Loc','Right Objects Lighter X-Loc','Left Objects Lighter X-Loc');
    else
        legend([e1 e2],'Y-Loc','X-Loc');
    end
    
    set(gca,'FontSize',fontSize,'FontName','Myriad Pro');

end

%% Export data for ANOVA
if ismember(graphs,4)
    
    GroupData = NaN(4*16*11,4);
    k = 1;
    
    for s = 1:length(GroupStruct) % for each subject
        
        T = GroupStruct{s}.TrialTable;
        vars = {'liftLocations'}; % for these columns only
        T2 = T{:,vars}; % select just this part of the table
        
        % concatenate all cells and convert to double
        allLifts = cat(1,T2{:});
        allDists = allLifts(:,2); % just the y-locs
        
        vars = {'liftWeights'}; % for these columns only
        T2 = T{:,vars}; % select just this part of the table
        allWeights = cat(2,T2{:});
        
        allData = cat(2,allDists,allWeights');
        
        binDists = unique(allDists);
        
        count = 1;
        binCounts = zeros(1,16*4);
        weightsOrd = zeros(1,16*4);
        distsOrd = zeros(1,16*4);
        for i = 1:length(binWeights) % for each weight
            for j = 1:length(binDists) % for each distance
                
                % find indices for these unique weights and distance
                idx = find(allData(:,2) == binWeights(i) & ...
                    allData(:,1) == binDists(j));
                
                % bin counts for each lift weight and distance
                binCounts(count) = length(idx);
                weightsOrd(count) = binWeights(i);
                distsOrd(count) = j;
                
                count = count + 1;
                   
            end
        end 
          
        % assign to group data
        GroupData(k:k+63,1) = s; % subject num 
        GroupData(k:k+63,2) = weightsOrd'; % weight
        GroupData(k:k+63,3) = distsOrd'; % distance
        GroupData(k:k+63,4) = binCounts'; % count
        
        k = k + 64;
        
    end
    
    % Convert to table
    tableForANOVA = array2table(GroupData,...
        'VariableNames',{'sub','weight','dist','count'});
    
    rm = fitrm(tableForANOVA,'count ~ weight*dist');
    
    
    
end

%% Plot the avg. lift weight for each participant
if ismember(graphs,5)
    
    GroupData = NaN(length(GroupStruct),1);
    
    for s = 1:length(GroupStruct) % for each subject
        
        T = GroupStruct{s}.TrialTable;
        vars = {'liftWeights'}; % for these columns only
        T2 = T{:,vars}; % select just this part of the table
        
        % concatenate all cells and convert to double
        allLifts = cat(2,T2{:});
        
        % bin counts for each lift weight
        GroupData(s) = mean(allLifts);
        
    end
    
    % calculate and plot group mean and std. error
    GraphData = mean(GroupData); % mean
    GraphErr = std(GroupData) / sqrt(length(GroupData)); % std error
      
    %%% PLOTTING %%%
    
    figure(4);
    clf
    f1 = gcf;
    f1.Position = [0 0 res(1) res(2)];
    hold on;
    
    b1 = barwitherr(GraphErr',GraphData');
    set(b1,'FaceColor',rgb('lightgray'),'EdgeColor','none')
    set(gca,'FontSize',18,'FontName','Myriad Pro');
    set(b1,'LineWidth',2);
    xticks(1);
    b1.BarWidth = 0.1;
    set(gca,'XTickLabel',{})
    title('Average Lift Weight');
    ylabel('Avg. Weight of Objects Lifted (kg)')

    % Draw individual subject traces
    hold on;
    for i = 1:size(GroupData,1) % for each subject
        % draw a line indicating the participant data
        p1 = scatter(1,GroupData(i,1),36,rgb('darkslategray'));
    end
    
    % if object was randomly selected
    randomWeight = sum(binWeights) / length(binWeights);
    line([0 2],[randomWeight randomWeight],'Color','k',...
        'LineStyle','--','LineWidth',2)
    
    % Statistics
    [h,p,~,~] = ttest(GroupData,randomWeight)
    
end

% Plot the correlation between peak rate of change in load force and peak
% height to object weight
if ismember(graphs,6)
    
    GroupData = NaN(length(GroupStruct),length(binWeights));
    
    for s = 1:length(GroupStruct) % for each subject
         
        T = GroupStruct{s}.TrialTable;
        vars = {'liftWeights','robotPeakDeltaForce',...
            'robotPeakLiftHeight','robotMeanForce'}; % for these columns only
        T2 = T{:,vars}; % select just this part of the table
        
        % concatenate all cells and convert to double
        allLifts1 = cat(2,T2{:,1});
        allLifts2 = cat(1,T2{:,2})';
        allLifts3 = cat(1,T2{:,3})';
        allLifts4 = cat(1,T2{:,4})';
        
        allLifts = [allLifts1;allLifts2;allLifts3;allLifts4];
        
        % plotting
        figure(5);
        f1 = gcf;
        f1.Position = [0 0 res(1) res(2)];
        subplot(3,4,s)
        hold on;
        
        plot(allLifts(1,:),allLifts(3,:),'LineStyle','none',...
            'Marker','.','Color',rgb('crimson'),...
            'LineWidth',0.5,'MarkerSize',10);
        
        xlimit = xlim();
        ylimit = ylim();
        line([0 xlimit(2)],[0 ylimit(2)],'LineStyle','--','LineWidth',0.5);
        
%         axis([0 1.4 0 10]);
        
        xlabel('Object Weight (kg)');
%         ylabel('Peak Force Change (N/s)');
        ylabel('Peak Lift Height (cm)');
%         ylabel('Peak Force (N)');
        
    end

end

if ismember(graphs,7)
  
    % init
    N = 14;
    numTrials = 60;
    GroupData = NaN(length(GroupStruct),numTrials);
    
    for s = 1:length(GroupStruct) % for each subject
        
        T = GroupStruct{s}.TrialTable;
        vars = {'liftWeights'}; % for these columns only
        T2 = T{:,vars}; % select just this part of the table
        
        weightForNLifts = NaN(length(T2),N); % init double
        
        for i = 1:length(T2) % for each trial
            for j = 1:N % for each lift
                try
                    weightForNLifts(i,j) = T2{i}(j);
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
        end
        
        % Mean and std. error for each lift
        meanNLifts = nanmean(weightForNLifts,2); % mean
        
        % Plotting
        figure(2);
        f1 = gcf;
        set(f1,'Position',[0,0,res(1),res(2)])
        hold on;
        
        plot(1:length(meanNLifts),meanNLifts,'-','Marker','none','Color',rgb('darkslategray'),...
            'LineWidth',0.5);
        
        GroupData(s,1:length(meanNLifts)) = meanNLifts';
     
    end
    
    % calculate and plot group means and std. error
    meanWeight = nanmean(GroupData); % mean
    errWeight = nanstd(GroupData) ./ sqrt(length(GroupStruct)); % std error
    
    errorbar(1:numTrials,meanWeight,errWeight,'LineWidth',3,'Color','k');
    
    % labels and axes
    title(sprintf('Avg. Weight for First %d Lifts',N));
    xlabel('Trial Number')
    ylabel('Avg. Weight of Objects (kg)');
    axis([0 numTrials+1 0 1]);
   
    % if object was randomly selected
        randomWeight = sum(binWeights) / length(binWeights);
    line([0 numTrials+1],[randomWeight randomWeight],'Color','k',...
        'LineStyle','--','LineWidth',2)
    
    set(gca,'FontSize',fontSize,'FontName','Myriad Pro');
    
end

%% Show Lift Height as a function of trial for Both Feature and Spatial condition
if ismember(graphs,8)
   
    FeatureDataLight = NaN(length(NEWSUBS),90);
    FeatureDataHeavy = NaN(length(NEWSUBS),90);
    SpatialDataLight = NaN(length(NEWSUBS),90);
    SpatialDataHeavy = NaN(length(NEWSUBS),90);

    for s = 1:length(NEWSUBS) % for each feature condition subject
        
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
        figure(2);
        f1 = gcf;
        set(f1,'Position',[0,0,res(1),res(2)])
        hold on;
        
%         plot(1:length(meanLiftHeightLight),meanLiftHeightLight,'-','Marker','.','Color',rgb('lightgray'),...
%             'LineWidth',0.5,'MarkerSize',10);
%         
%         plot(1:length(meanLiftHeightHeavy),meanLiftHeightHeavy,'-','Marker','.','Color',rgb('gray'),...
%             'LineWidth',0.5,'MarkerSize',10);
        
        FeatureDataLight(s,1:length(meanLiftHeightLight)) = meanLiftHeightLight';
        FeatureDataHeavy(s,1:length(meanLiftHeightHeavy)) = meanLiftHeightHeavy';
    
    end
    
    for s = 1:length(OLDSUBS) % for each spatial condition subject
      
        z = length(NEWSUBS) + s;
        
        T = GroupStruct{z}.TrialTable;
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

%         plot(1:length(meanLiftHeightLight),meanLiftHeightLight,'-','Marker','.','Color',rgb('lightgray'),...
%             'LineWidth',0.5,'MarkerSize',10);
%         
%         plot(1:length(meanLiftHeightHeavy),meanLiftHeightHeavy,'-','Marker','.','Color',rgb('gray'),...
%             'LineWidth',0.5,'MarkerSize',10);
        
        SpatialDataLight(s,1:length(meanLiftHeightLight)) = meanLiftHeightLight';
        SpatialDataHeavy(s,1:length(meanLiftHeightHeavy)) = meanLiftHeightHeavy';
    
    end
    
    % calculate and plot mean and std. error for blocks of trials
    n = 10; % Number of elements in each block
    % average across block size n
    meanLightFeatureG = NaN(length(NEWSUBS),n-1);
    meanHeavyFeatureG = NaN(length(NEWSUBS),n-1);
    meanLightSpatialG = NaN(length(NEWSUBS),n-1);
    meanHeavySpatialG = NaN(length(NEWSUBS),n-1);
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
    errLightFeature = nanstd(meanLightFeatureG) ./ sqrt(length(NEWSUBS)); % std error
    meanHeavyFeature = nanmean(meanHeavyFeatureG); % mean
    errHeavyFeature = nanstd(meanHeavyFeatureG) ./ sqrt(length(NEWSUBS)); % std error
    
    meanLightSpatial = nanmean(meanLightSpatialG); % mean
    errLightSpatial = nanstd(meanLightSpatialG) ./ sqrt(length(NEWSUBS)); % std error
    meanHeavySpatial = nanmean(meanHeavySpatialG); % mean
    errHeavySpatial = nanstd(meanHeavySpatialG) ./ sqrt(length(NEWSUBS)); % std error
    
   % calculate across trial means
    CmeanLightFeature = mean(nanmean(meanLightFeatureG,2)); % mean
    CerrLightFeature = nanstd(nanmean(meanLightFeatureG,2)) ./ sqrt(length(NEWSUBS)); % std error
    CmeanHeavyFeature = mean(nanmean(meanHeavyFeatureG,2)); % mean
    CerrHeavyFeature = nanstd(nanmean(meanHeavyFeatureG,2)) ./ sqrt(length(NEWSUBS)); % std error
    
    CmeanLightSpatial = mean(nanmean(meanLightSpatialG,2)); % mean
    CerrLightSpatial = nanstd(nanmean(meanLightSpatialG,2)) ./ sqrt(length(NEWSUBS)); % std error
    CmeanHeavySpatial = mean(nanmean(meanHeavySpatialG,2)); % mean
    CerrHeavySpatial = nanstd(nanmean(meanHeavySpatialG,2)) ./ sqrt(length(NEWSUBS)); % std error
    
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
if ismember(graphs,9)
    
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
if ismember(graphs,10)
    
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


