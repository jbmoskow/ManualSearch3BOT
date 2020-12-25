% Plots the individual trial by trial data showing target locations
% and hand position. Also generates a subplot below containing the X, Y,
% and Z hand position of the robot

%% Creatures plot to show exemplar trial for 3BOTSearch
function [] = ThreeBotSearch_Exemplar(SUB,trialNum,RPX,RPY,RPZ,idxLifts,TrialData)

% SUB - subject ID
% trialNum - current trial number
% RPX - robot X position
% RPY - robot Y position
% RPZ - robot Z position
% idxLifts - Nx2 array containing onset and offset time for each lift

% TrialData contains
% ObjectHomePosX - X locations of the targets
% ObjectHomePosY - Y locations of the targets
% ObjectRadius - radius of the locations of the targets
% ObjectWeight - weight of the targets
% ObjectTarget - 0 for non-targets, and 1 otherwise

figure(102)
clf
f1 = gcf;
set(f1,'Position',[0,50,800/2,600])

% init
RGB_LIGHTGREEN = [0/255 255/255 0/255]; 
RGB_DARKGREEN = [0/255 100/255 0/255];
centers = [TrialData.ObjectHomePosX' TrialData.ObjectHomePosY'];

%% get condition type
if ismember(SUB,{'MC1','CVB1','JS1','LZ1','MS1','TG1','ET1','KC1','IH1',...
        'TSJ1','HA1','MB1'})
    cond = 1; % this is a feature condition subject
elseif ismember(SUB,{'SP1','YZ1','ET2','AK1','HB1','YZ2','AL1','KF1','AN1','NL1',...
        'CO1','BI1'})
    cond = 2; % this is a spatial condition subject
end

if cond == 1
    if ismember(SUB,{'MC1','LZ1','MS1','JS1','CVB1','MB1'})
        lightColor = RGB_DARKGREEN;
        heavyColor = RGB_LIGHTGREEN;
    elseif ismember(SUB,{'TG1','ET1','KC1','IH1','TSJ1','HA1'})
        lightColor = RGB_LIGHTGREEN;
        heavyColor = RGB_DARKGREEN; 
    end
end
    
% assign colors to target objects
colors = zeros(length(TrialData.ObjectWeight),3);
for i = 1:length(TrialData.ObjectWeight)
   if TrialData.ObjectWeight(i) == 0.1 
       colors(i,:) = lightColor;
   else
       colors(i,:) = heavyColor;
   end
end

% get index of light objects to use for plotting
idxLight = TrialData.ObjectWeight == 0.1;

% get index of last lift
lastLift = idxLifts(end,2);

%% Draw stimuli
sh1 = subplot(2,1,1);
hold on;
viscircles(centers(idxLight,:),TrialData.ObjectRadius(idxLight),'Color',lightColor); % light objects
viscircles(centers(~idxLight,:),TrialData.ObjectRadius(~idxLight),'Color',heavyColor); % heavy objects


%% Plot X and Y hand position
plot(RPX(1:lastLift),RPY(1:lastLift),'k--');
axis([-15 15 -15 15])
axis square

title(sprintf('Subject %s, Trial Number: %d',SUB,trialNum));

%% Plot X Y,and Z positions
subplot(2,1,2)
hold on;
plot(RPX(1:lastLift),'k--');
plot(RPY(1:lastLift),'b--');
plot(RPZ(1:lastLift),'r--');
ylimit = ylim;
for i = 1:size(idxLifts,1) % for each lift plot the onset time
    line([idxLifts(i,1) idxLifts(i,1)],[ylimit(1) ylimit(2)],...
        'LineStyle','--','LineWidth',2);
end
% x1 = xticklabels;
% for i = 1:length(x1)
%     newlabels{i} = num2str(str2double(x1{i}) / 1000);
% end
% xticklabels(newlabels); % convert x-axis to seconds
xlabel('Time (ms)')
ylabel('Robot Handle Location (cm)')

legend('X-pos','Y-pos','Z-pos','Lift Onsets')

sh1.Position = [0 0.48 1 0.45];


end