%% Summary
% Matlab code uploaded to github in support of the paper submitted to 
% geophyiscal research letters titled:
% Locating the top of the saturated zone with airborne electromagnetic data
% the code and the accompanying data in the following github repo can be
% used to recreate all of the figures and results from the aforementioned
% paper.

% Author: Noah Dewar, Stanford University, Department of Geophysics
% March 19th, 2019

% License 
% This work is licensed under the Creative Commons 
% Attribution-NonCommercial-ShareAlike 4.0 International License. To view 
% a copy of this license, visit http://creativecommons.org/licenses/by-nc-sa/4.0/ 
% or send a letter to Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

%% Setup
% set up some required variables, import the data
clear

% load the data
load('MatlabData.mat');

% collect the layer depths
layerDepths=[0;cumsum(THICKNESS(2:end))];
layerTops=layerDepths(1:end);
layerBottoms=[layerDepths(2:end);700];

% depth window to look within for the water table, in meters below the
% surface
depthWindow=[30 100];

% set the windowing distance, window distance increases by 500 instead of
% 50 between 1000 and 6000 to reduce data volume for github repo
windowDistance=[60 100:50:1000 1500:500:6000];

% set some more stuff
methods_to_try=12;

% clean up
clear THICKNESS

%% WINDOWING
% do the windowing, estimate the WTE

% preclean
clear estWatertable statsArray

% set some arrays up
statsArray=cell(methods_to_try,numel(estimationLocations{:,1}));
estWatertable=zeros(numel(windowDistance),numel(estimationLocations{:,1}),methods_to_try/2);

for k=1:numel(windowDistance)    
for i=1:numel(estimationLocations{:,1})
    
    % load from stored arrays the res and subsample depths
    tempRes=resArray{k,i};
    subsampleDepths=subsampleDepthsArray{k,i};
    
    % find some stats like the spread between the min and the max, and
    % the distance between the 25th and 75th percentile, and their
    % gradients  
    statsArray{1,i}=std(tempRes);
    statsArray{2,i}=diff(statsArray{1,i});
    statsArray{3,i}=mean(tempRes);
    statsArray{4,i}=diff(statsArray{3,i});
    statsArray{5,i}=min(tempRes);
    statsArray{6,i}=diff(statsArray{5,i});
    statsArray{7,i}=max(tempRes);
    statsArray{8,i}=diff(statsArray{7,i});
    statsArray{9,i}=max(tempRes)-min(tempRes);
    statsArray{10,i}=diff(statsArray{9,i});
    statsArray{11,i}=prctile(tempRes,75)-prctile(tempRes,25);
    statsArray{12,i}=diff(statsArray{11,i});
    
    % make an index for within what depths to look for the biggest thing
    depth_index=double(subsampleDepths<=subsampleDepths(1)-depthWindow(1))+double(subsampleDepths>=subsampleDepths(1)-depthWindow(2));
    
    % remove the last entry for both the subsample depths and the depth
    % index as the gradient removed a row
    depth_index(end)=[];
    subsampleDepths(end)=[];
    
    % put the estimated watertable depth as the biggest change in
    % percentile spread after the first couple big changes
    for m=1:methods_to_try/2
        est_depth=subsampleDepths(abs(statsArray{m*2,i})==max(abs(statsArray{m*2,i}(depth_index==2))));
        estWatertable(k,i,m)=est_depth(1); 
    end
    
end
end

% clean up
clear current_*  i j k m ii tempRes ...
   subsample_depths est_depth

%% Error Calc
% see how well it did

% clear some variables just cause
clear bestOne

% set how far from the wells the estimations points are allowed to be and
% set up some arrays
validationBuffer=2000;
validationDate=datetime(2015,9,10);
index=1:methods_to_try/2;

% set what windows to use and the validation dataset
windows_to_use=1:numel(windowDistance);
estimationWells_used=estimationWells(estimationLocations{1:numel(estimationWells{:,1}),end}<validationBuffer,:);
estimationWellsTID_used=estimationWellsTID(estimationLocations{numel(estimationWells{:,1})+1:end,end}<validationBuffer,:);

% filter to only the wells that were measured in october
estimationWells_used=estimationWells_used(estimationWells_used{:,8}>validationDate,:);
estimationWellsTID_used=estimationWellsTID_used(estimationWellsTID_used{:,5}>validationDate,:);

% remove wells with Q for questionable quality
estimationWellsTID_used(estimationWellsTID_used{:,6}==8,:)=[];
validationData=[estimationWells_used.WSEL*0.3048;estimationWellsTID_used.WSEL]';

% get the DEM for the wells
estimationWells_usedDem=[estimationWells_used{:,13}*0.3048;estimationWellsTID_used.GSE];

% restrict the estimation locations to be within the validation buffer
index_to_use=estimationLocations{:,end}<validationBuffer;
dateIndex=[estimationWells{:,8}>validationDate;estimationWellsTID{:,5}>validationDate];
index_to_use=logical(max(index_to_use+dateIndex-1,0));

% Remove the same wells that have a questionable quality rating given by
% TID as above, but from the estimation locations
index_to_use(11)=0; index_to_use(17)=0; index_to_use(30)=0;

% bootstrap to see how many calibration wells are needed
number_to_use=1:10;
bootstrap=1;
bootstrapIters=10000;

% set up some arrays
bestOne=zeros(5,numel(windows_to_use));
bestArray=cell(bootstrapIters,numel(number_to_use));
meanError=zeros(length(estWatertable(1,1,:)),numel(windows_to_use));
madError=zeros(length(estWatertable(1,1,:)),numel(windows_to_use));
rmsError=zeros(length(estWatertable(1,1,:)),numel(windows_to_use));
rSquared=zeros(length(estWatertable(1,1,:)),numel(windows_to_use));
if bootstrap
    for k=1:numel(number_to_use)
        for l=1:bootstrapIters
        curr_index=randi(26,number_to_use(k),1);
        validationDataCurrent=validationData(curr_index);
            for m=1:numel(windows_to_use)
                currWatertable=estWatertable(windows_to_use(m),index_to_use,:);
                currWatertable=currWatertable(:,curr_index,:);
                for i=1:length(estWatertable(1,1,:))
                    meanError(i,m)=mean(abs(currWatertable(:,:,i)-validationDataCurrent));
                    madError(i,m)=median(abs(currWatertable(:,:,i)-validationDataCurrent));
                    rmsError(i,m)=rms(currWatertable(:,:,i)-validationDataCurrent);
                end
            bestOne(1,m)=min(rmsError(:,m));
            bestOne(2,m)=meanError(find(min(rmsError(:,m))==rmsError(:,m),1),m);
            bestOne(3,m)=madError(find(min(rmsError(:,m))==rmsError(:,m),1),m);
            bestOne(4,m)=rSquared(find(min(rmsError(:,m))==rmsError(:,m),1),m);
            bestOne(5,m)=index(find(min(rmsError(:,m))==rmsError(:,m),1));
            end
        bestArray{l,k}=bestOne; 
        end
    end
else
    for m=1:numel(windows_to_use)
        currWatertable=estWatertable(windows_to_use(m),index_to_use,:);
        for i=1:length(estWatertable(1,1,:))
            meanError(i,m)=mean(abs(estWatertable(windows_to_use(m),index_to_use,i)-validationData));
            madError(i,m)=median(abs(estWatertable(windows_to_use(m),index_to_use,i)-validationData));
            rmsError(i,m)=rms(estWatertable(windows_to_use(m),index_to_use,i)-validationData);
            mdl=fitlm(validationData,estWatertable(windows_to_use(m),index_to_use,i));
            rSquared(i,m)=mdl.Rsquared.ordinary;
        end
    bestOne(1,m)=min(rmsError(:,m));
    bestOne(2,m)=meanError(find(min(rmsError(:,m))==rmsError(:,m),1),m);
    bestOne(3,m)=madError(find(min(rmsError(:,m))==rmsError(:,m),1),m);
    bestOne(4,m)=rSquared(find(min(rmsError(:,m))==rmsError(:,m),1),m);
    bestOne(5,m)=index(find(min(rmsError(:,m))==rmsError(:,m),1));
    end
end

% best result is from the min gradient (3) and 75th to 25h spread
bestMethod=bestOne(5,min(bestOne(1,:))==bestOne(1,:));
finalWindow_to_use=sum((1:numel(windowDistance)).*double(min(bestOne(1,:))==bestOne(1,:)));
clear final_estimate
finalEstimate=estWatertable(finalWindow_to_use,index_to_use,bestMethod);

% make variables for depth below the surface
finalEstimateDepth=estimationLocations.dem(index_to_use)'-finalEstimate;
validation_data_depth=estimationWells_usedDem'-validationData;

% find out the size of the skytem pickle for each estimation depth/location
layerThickness=layerBottoms-layerTops;
finalEstimateResolution=zeros(1,numel(finalEstimate));
for i=1:numel(finalEstimate)
    finalEstimateResolution(i)=layerThickness((abs(layerBottoms-finalEstimateDepth(i)))==min(abs(layerBottoms-finalEstimateDepth(i))));
end

% clean up
clear current_* nan_index i j k m avg distances* index mdl *_index curr* ...
    

%% Bootstrap Analysis
% loop through the bootstrap iters and see what was the best
index=1:numel(windowDistance)+1;
minArrayRMS=zeros(bootstrapIters,numel(number_to_use));
minArrayMethod=zeros(bootstrapIters,numel(number_to_use));
minArrayWindow=zeros(bootstrapIters,numel(number_to_use));
method_counts=zeros(8,numel(number_to_use));
window_counts=zeros(8,numel(number_to_use));
for j=1:numel(number_to_use)
for i=1:bootstrapIters
    minArrayRMS(i,j)=min(bestArray{i,j}(1,:));
    temp=bestArray{i,j}(5,min(bestArray{i,j}(1,:))==bestArray{i,j}(1,:));
    minArrayMethod(i,j)=temp(1);
    temp=index(min(bestArray{i,j}(1,1:end-1))==bestArray{i,j}(1,1:end-1));
    minArrayWindow(i,j)=temp(1);
end
for i=1:8
    method_counts(i,j)=sum(minArrayMethod(minArrayMethod(:,j)==i,j)./i);
end
for i=1:numel(windowDistance)+1
    window_counts(i,j)=sum(minArrayWindow(minArrayWindow(:,j)==i,j)./i);
end
end

% clean up
clear i j minArrayWindow minArrayMethod

%% Figures
% some details or additional graphics may have been added to figures in
% Microsoft Powerpoint. The figures made by this code may not exactly match
% the figures presented in the submitted version of the manuscript

% code for figure 2
lineColors=lines(6);
figure(2)
set(gcf,'pos',[100 100 1000 900]);
plot(windowDistance,rmsError(2,1:end),'LineWidth',2,'Color',lineColors(2,:)); hold on
plot(windowDistance,rmsError(3,1:end),'LineWidth',2,'Color',lineColors(3,:));
plot(windowDistance,rmsError(4,1:end),'LineWidth',2,'Color',lineColors(4,:));
plot(windowDistance,rmsError(6,1:end),'LineWidth',2,'Color',lineColors(6,:)); 
plot(windowDistance,rmsError(5,1:end),'LineWidth',2,'Color',lineColors(5,:)); 
plot(windowDistance,rmsError(1,1:end),'LineWidth',2,'Color',lineColors(1,:));
plot([min(windowDistance)-1000 max(windowDistance)],[min(min(rmsError)) min(min(rmsError))],'--k')
plot(windowDistance,rmsError(6,1:end),'LineWidth',2,'Color',lineColors(6,:)); 
legend({'Minimum','Mean','Maximum','75th to 25th percentile difference','Maximum minimum difference','Standard deviation',['Lowest RMS error, ',num2str(min(round(100*bestOne(1,:))/100)),' m']},'Location','SouthEast');
ylabel('RMS error, meters'); xlabel('Search radius, meters')
grid on; box on; set(gca,'FontSize',16);
xlim([min(windowDistance) max(windowDistance)]);

% code for figure 3
figure(3)
set(gcf,'pos',[100 100 1000 900]);
scatter(validationData,finalEstimate,60,estimationLocations.distance(index_to_use),'filled'); hold on;
plot(0:0.1:100,0:0.1:100,'-r')
errorbar(validationData,finalEstimate,finalEstimateResolution,'b','LineStyle','none');
h=colorbar; ylabel(h, 'Separation distance, meters');
legend({'Data','One to one line','SkyTEM pixel size'},'location','northwest')
ylabel('Estimated TSZ elevation MSL, meters'); xlabel('Measured water table elevation MSL, meters')
grid on; box on; set(gca,'FontSize',16); colormap(jet)
xlim([0 60])
ylim([0 60])
