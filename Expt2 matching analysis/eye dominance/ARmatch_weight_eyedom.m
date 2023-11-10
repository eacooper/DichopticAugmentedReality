%Experiment 2 matching result
%find and plot the weight of the higher contrast eye, separated by
%different stimulus patterns

clear all;
close all;
addpath('../../data/Expt2_AR');
addpath('../');
%set(groot,'defaultfigureposition',[275,243,1158,420]);
filenames = dir(['../../data/Expt2_AR/BAR_*','mat']);
filenames = {filenames.name};
N = size(filenames,2); %number of subjects

%initiate grid search space 
stepsize = 0.01;
weights  = 0:stepsize:1;

ARweights = [];

%load file with whether LE is dominant
load('ARexpt_eyedom.mat'); %1 = LE dominant

%% process the data
for s = 1:N

    %load data
    filename = filenames{s};
    load(filename);
    data = [dat2.stim dat2.resp(:,1)];
    
    %find weight for high contrast stim, separated by LE or RE seeing high
    %contrast stim
    for eyeinfotype = 1:5
        eyeinfo = eyedom(s,eyeinfotype); %get eye dominance info based on each stim (1-4) or all stim (5)

        if eyeinfo ==1 %LE dominant
            DEdata = data(data(:,2)>data(:,3),:);
            NDEdata = data(data(:,2)<data(:,3),:);
        else
            DEdata = data(data(:,2)<data(:,3),:); %RE higher
            NDEdata = data(data(:,2)>data(:,3),:);
        end
        
        
        for eye = 1:2 %separate data based on dominant eye saw high contrast or not
            if eye ==1
                eyedata = DEdata;
            else
                eyedata = NDEdata;
            end
            
            for stim_ind = 1:5 %4 stimulus patterns, 5 as taking all 4 types
                
                %get the subset of data for this condition
                if stim_ind~=5
                    dataToFit= eyedata(eyedata(:,1)==stim_ind,:);
                else
                    dataToFit= eyedata;
                end
                
                dataToFit = dataToFit(dataToFit(:,2)~=dataToFit(:,3),:); % removing catch trials (nondichoptic)
                
                dataToFitlow = min(dataToFit(:,2:3),[],2);
                dataToFithigh = max(dataToFit(:,2:3),[],2);
                
                %do the fitting
                ModelPred=genBino(dataToFitlow,dataToFithigh,weights);
                
                humandata = dataToFit(:,4);
                subjdata = ones(size(ModelPred)).*humandata; %repeat human data for multiple cols
                
                %compare to the human data by minimizing RMSE
                diffsq = (subjdata - ModelPred).^2;
                rmse_matrix = sqrt(mean(diffsq,1));
                minRMSE = min(rmse_matrix);
                bestW_ind = find(rmse_matrix == minRMSE);
                w = weights(bestW_ind);
                
                %store the result
                ARweights = [ARweights; s, eyeinfotype, eye, stim_ind, w, minRMSE];
            end
        end
    end
end



%% make figure
figure(2);
eyeinfotype = 5; %take eye dominance info based on one of the stim (1-4) or all the stim (5)
DE = [];
NDE = [];
for stim = 1:4 %for each stimulus type
    
    %dominant eye plot
    domeye = ARweights(ARweights(:,2)==eyeinfotype & ARweights(:,3) == 1 & ARweights(:,4)==stim,:);
    subplot(1,2,1); hold on;
    plot(ones(1,size(domeye,1))*(stim)-0.25+rand(1,size(domeye,1))*0.4,domeye(:,5),'k.');
    DE = [DE,domeye(:,5)];
    title('dominant eye high contrast');
    axis square;
    
    %non-dominant eye plot
    nondomeye = ARweights(ARweights(:,2)==eyeinfotype & ARweights(:,3) == 2 & ARweights(:,4)==stim,:);
    subplot(1,2,2); hold on;
    plot(ones(1,size(nondomeye,1))*(stim)-0.25+rand(1,size(nondomeye,1))*0.4,nondomeye(:,5),'k.');
    NDE = [NDE,nondomeye(:,5)];
    title('non-dominant eye high contrast');
    axis square;

    
end

%add in box and whisker plot on top of individual data
subplot(1,2,1);
boxplot(DE,'Colors',[0.99 0.64 0.32; 0.88 0 0.5;0.52 0.67 0.36;0 0.4 0.99],'Symbol','');
xlim([0.5 4.5]);
ylim([-0.1 1.1]);
set(gca,'XTickLabel',{"grating","noise","simple","complex"});
subplot(1,2,2);
boxplot(NDE,'Colors',[0.99 0.64 0.32; 0.88 0 0.5;0.52 0.67 0.36;0 0.4 0.99],'Symbol','');
xlim([0.5 4.5]);
ylim([-0.1 1.1]);
set(gca,'XTickLabel',{"grating","noise","simple","complex"});





