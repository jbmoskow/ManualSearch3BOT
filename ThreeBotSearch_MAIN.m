% MAIN Script to import, analyze, and graph data from the
% 3BOT Object Search experiment

% Coded in July 2019
% Modified in January 2020 to accomodate new data
% J. Moskowitz

%% select subs

% 3BOT Lifting Participants - Feature Condition
FeatureCondition = {'MC1','CVB1','JS1','LZ1','MS1','TG1','ET1','KC1','IH1',...
    'TSJ1','HA1','MB1'};

% 3BOT Lifting Participants - Spatial Condition
SpatialCondition = {'SP1','YZ1','ET2','AK1','HB1','YZ2','AL1','KF1','AN1','NL1',...
    'CO1','BI1'};

FeatureIdx = 1:12;
SpatialIdx = 13:24;
SUBS = [FeatureCondition SpatialCondition];

%% select plots  
graphs = [1 2 3 4]; % plots
% 1 = proportion of lifts for lightest object
% 2 = prop. of lifts for lightest object for first N lifts
% 3 = proportion of lifts by y-loc and x-loc
% 4 = avg. lift height by weight for each trial
% 5 = mean search times
% 6 = mean lift number
saveGraphs = 0; % 1 = save graphs to .pdf, 0 = do not save
exemplar = 0; % run ThreeBotSearch_Analysis
procDir = 'C:/Users/Joshua/Desktop/Data/';
workDir = 'C:/Users/Joshua/Desktop';
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
    [GroupData] = ThreeBotSearch_Grapher(SUBS,FeatureIdx,SpatialIdx,GroupStruct,graphs(x));
    
    if saveGraphs == 1
        % export to pdf
        
        export_fig('3BOTSearch.pdf','-append','-painters','-p15');
        if mod(x,5) == 0 % close all currently open figures to refresh memory
            close all;
        end
    end
    
end

