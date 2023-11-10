%Experiment 2 matching result
%find and plot the weight of the higher contrast eye, separated by
%interocualr contrast ratio (ICR)

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

ICRlist = [1 2 4 5]; %ICR values (5 indicates monocular)

for s = 1:N

    filename = filenames{s};
    load(filename);
    % calculate ICR for each trial
    ICRval = round(max(dat2.stim(:,2:3),[],2)./min(dat2.stim(:,2:3),[],2));
    ICRval(isinf(ICRval)) = 5; %recode inf (monocular) to 5
    
    data = [dat2.stim dat2.resp(:,1) ICRval];


    for ICRind = 2:4 %do dichoptic trials only

        %get the subset of data for this condition
        dataToFit= data(data(:,5)==ICRlist(ICRind),:);

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
        ARweights = [ARweights; s, ICRlist(ICRind), w, minRMSE];
        subjList = [subjList; ['S',num2str(s)]];


    end
end


%% save weight results for R
T = array2table(ARweights);
% Assign the specific headings
T.Properties.VariableNames(1:4) = {'Subj','ICRcat','Weight','RMSE'};
ICRtxt = {'one','two','three','four','mono'};
ss = ICRtxt(T.ICRcat)';
T.Stimulus = categorical(ss);
T.Subj = subjList;

writetable(T,'./R stats/AR_weightByICR.csv');


%% plot the weights

figure(1); subplot(1,3,1); hold on;

for ICRind = 1:4

    stimdata = ARweights(ARweights(:,2)==ICRlist(ICRind),:);

    %plot the mean and 95% CI
    datamean = mean(stimdata(:,3));
    dataerror = std(stimdata(:,3))*1.96/sqrt(size(stimdata,1));
    errorbar(ICRlist(ICRind),datamean,dataerror,'ko','MarkerFaceColor',[0 0 0],'LineWidth',2,'MarkerSize',6);
    plot(ones(1,size(stimdata,1))*ICRlist(ICRind)-0.1+rand(1,size(stimdata,1))*0.2,stimdata(:,3),'.','MarkerEdgeColor',[0.5 0.5 0.5]);

end


xlim([0.5 5.5]);
ylim([-0.1 1.1]);
yticks([0 0.5 1]);
xticks([1 2 3 4 5]);
set(gca,'xticklabel',{'1','2','','4','mono'});
%xtickangle(45);
ylabel('Weight of High Contrast');
xlabel('Interocular contrast ratio (ICR)');
axis square;

plot([0.5 5.5],[1 1],'k--');
plot([0.5 5.5],[0 0],'k--');
plot([0.5 5.5],[0.5 0.5],'k--');
box on;

