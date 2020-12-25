% MAIN Script to import, analyze, and graph data from the
% 3BOT Object Search experiment

% Coded in July 2019
% Modified in January 2020 to accomodate new data from Jennifer thesis project
% J. Moskowitz

%% select subs

% 3BOT Lifting (Jennifer) Participants - Feature Condition
FeatureCondition = {'MC1','CVB1','JS1','LZ1','MS1','TG1','ET1','KC1','IH1',...
    'TSJ1','HA1','MB1'};

% 3BOT Lifting (Jennifer) Participants - Spatial Condition
SpatialCondition = {'SP1','YZ1','ET2','AK1','HB1','YZ2','AL1','KF1','AN1','NL1',...
    'CO1','BI1'};

SUBS = FeatureCondition;

% SUBS = {'MC1','CVB1','JS1','LZ1','MS1','TG1','ET1','KC1','IH1',...
%     'TSJ1','HA1','MB1','SP1','YZ1','ET2','AK1','HB1','YZ2','AL1','KF1','AN1','NL1',...
%     'CO1','BI1'};

% exemplar
% SUBS = {'CVB1'}; % trial 28

% Summer 2019 SUBS
NEWSUBS = {'MC1','CVB1','JS1','LZ1','MS1','TG1','ET1','KC1','IH1',...
     'TSJ1','HA1','MB1'};
OLDSUBS = {'SP1','YZ1','ET2','AK1','HB1','YZ2','AL1','KF1','AN1','NL1',...
     'CO1','BI1'};

%% select plots  
graphs = []; % plots
% 1 = proportion of lifts for lightest object from Jennifer thesis
% 11 = proportion of lifts for lightest object from both studies (Jennifer
% and summer 2019)
% 111 - proportion of lifts for lightest object from both spatial and
% feature conditions
% 2 = prop. of lifts for lightest object for first N lifts
% 3 = proportion of lifts by y-loc
% 3.1 = proportion of lifts by x-loc
% 3.2 = avg. x-loc and y-loc for first N lifts
% 4 = ANOVA
% 5 = avg. lift weight
% 6 = correlation between object weight and load force delta/lift height
% 7 = avg. weight of lift objects for first 3 lifts across trials
% 8 = avg. lift height by weight for each trial
% 9 = mean search times
% 10 = mean lift number
saveGraphs = 0; % 1 = save graphs to .pdf, 0 = do not save
exemplar = 0; % run ThreeBotSearch_Analysis
rawDir = 'C:/Users/Josh/Desktop/3Bot/Data/Jennifer';
procDir = 'C:/Users/Josh/Desktop/Dropbox/Projects/3BOT/Data/';
workDir = 'C:/Users/Josh/Desktop/Dropbox/Projects/3BOT/MATLAB code/Jennifer';
%% main

for s = 1:length(SUBS)
    
    if (exist([procDir SUBS{s} 'proc.mat'],'file') ~= 2) || exemplar % check if proc. data exists
        % perform analysis and output results to structure
        AnalysisStruct2.S{s} = ThreeBotSearch_Analysis(SUBS{s},rawDir);
        AnalysisStruct = AnalysisStruct2.S{s};
        save([procDir SUBS{s} 'proc.mat'],'AnalysisStruct');
    else
        cd(procDir);
        load([SUBS{s} 'proc.mat']);
    end
    
    disp(['Current S# loaded:' num2str(s)]);
    
    GroupStruct{s} = AnalysisStruct;
    
end

for x = 1:length(graphs)
    
    cd(workDir);
    [GroupData] = ThreeBotSearch_Grapher(SUBS,GroupStruct,graphs(x),NEWSUBS,OLDSUBS);
    
    if saveGraphs == 1
        % export to pdf
        
        export_fig('3BOTSearch.pdf','-append','-painters','-p15');
        if mod(x,5) == 0 % close all currently open figures to refresh memory
            close all;
        end
    end
    
end

