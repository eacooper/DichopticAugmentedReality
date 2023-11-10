%Experiment 2 matching result
%find and plot the weight of the higher contrast eye, separated by
%different stimulus patterns

clear all;
close all;
addpath('../data/Expt2_AR');
set(groot,'defaultfigureposition',[275,243,1158,420]);
filenames = dir(['../data/Expt2_AR/BAR_*','mat']);
filenames = {filenames.name};
N = size(filenames,2); %number of subjects

%initiate grid search space 
stepsize = 0.01;
weights  = 0:stepsize:1;

ARweights = [];
subjList = {};

for s = 1:N

    filename = filenames{s};
    load(filename);

    data = [dat2.stim dat2.resp(:,1)];
    
    for stim_ind = 1:5 %4 stimulus patterns, 5 as taking all 4 types combined
        
        %get the subset of data for this condition
        if stim_ind~=5
            dataToFit= data(data(:,1)==stim_ind,:);
        else
            dataToFit= data;
        end
        
        dataToFit = dataToFit(dataToFit(:,2)~=dataToFit(:,3),:); % removing catch trials (nondichoptic)
        
        dataToFitlow = min(dataToFit(:,2:3),[],2); %low contrast stim
        dataToFithigh = max(dataToFit(:,2:3),[],2); %high contrast stim
        
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
        ARweights = [ARweights; s, stim_ind, w, minRMSE];

        if stim_ind~=5 %create a list of subjects for R analysis
            subjList = [subjList; ['S',num2str(s)]];
        end
        
    end
end

%% save weight results for R
T_temp = ARweights(ARweights(:,2)~=5,:); %remove all stim type combined weight
T = array2table(T_temp);
% Assign the specific headings
T.Properties.VariableNames(1:4) = {'Subj','Stimulus','Weight','RMSE'};
Stimtxt = {'grating','noise','simple','complex'};
ss = Stimtxt(T.Stimulus)';
T.Stimulus = categorical(ss);
T.Subj = subjList;

writetable(T,'./R stats/AR_weightBypattern.csv');


%% plot the weights
for stim = 1:4
    stimdata = ARweights(ARweights(:,2)==stim,:);
    subplot(1,2,1); hold on;
    
    %plot the mean and 95% CI
    datamean = mean(stimdata(:,3));
    dataerror = std(stimdata(:,3))*1.96/sqrt(size(stimdata,1));
    errorbar(stim,datamean,dataerror,'ko','MarkerFaceColor',[0 0 0],'LineWidth',2,'MarkerSize',6);    
    plot(ones(1,size(stimdata,1))*(stim)-0.1+rand(1,size(stimdata,1))*0.2,stimdata(:,3),'.','MarkerEdgeColor',[0.5 0.5 0.5]);
    disp(datamean)
end

xlim([0.5 4.5]);
ylim([-0.1 1.1]);
yticks([0 0.5 1])
xticks([1 2 3 4]);
plot([0.5 5.5],[1 1],'k--');
plot([0.5 5.5],[0 0],'k--');
plot([0.5 5.5],[0.5 0.5],'k--');

set(gca,'xticklabel',{'grating','noise','simple','complex'});
xtickangle(45);
ylabel('Weight of High Contrast');
axis square;
box on;



