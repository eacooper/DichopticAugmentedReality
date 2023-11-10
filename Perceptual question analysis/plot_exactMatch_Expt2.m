% Experiment 2 exact match plot,the proportion of exact match based on ICR or stimulus type

clear all;
close all;
%figure size
set(groot,'defaultfigureposition',[275,243,1158,420]);
%data paths
addpath('../data/Expt2_AR/'); 
filenames = dir(['../data/Expt2_AR/*','mat']);
filenames = {filenames.name};
N = size(filenames,2); %number of subjects

%% ICR ratio

%create a list of ratios
contrast = [0 0.13 0.25 0.5];
ICR_ratio = [];

c = [];
for ii = 1:size(contrast,2)
    for jj = ii:size(contrast,2)
        c = [c; contrast(ii), contrast(jj)];
    end
end
ICR_ratio = c(:,2)./c(:,1);%higher contrast / lower contrast
ICR_ratio(isnan(ICR_ratio)) = []; %remove 0/0
ICR_ratio = round(ICR_ratio); %rounding the ratios
ICR_ratio = unique(ICR_ratio);  %find the unique ICRs

%store subject results
%first column is the ICRs, columns 2:end will be each subject's data
allsubj_ICR = ICR_ratio;

for s = 1:size(filenames,2)
    
    filename = filenames{s};
    load(filename);

    %add in ICR for each trial
    trial_ICR = round(max(dat2.stim(:,2:3),[],2)./(min(dat2.stim(:,2:3),[],2)));
    sdata = [dat2.stim, trial_ICR, dat2.resp(:,2)];
    
    for ICR = ICR_ratio'
        
        subdata = sdata(sdata(:,4)==ICR,:);
        propMatch = sum(subdata(:,5)==1)/size(subdata,1);
        
        allsubj_ICR(allsubj_ICR(:,1)==ICR,s+1) = propMatch;

    end
    
end

allsubj_ICR(4,1) = 5; %change Inf to 5 for plotting
%make plot
figure(1); 
subplot(1,3,1);
hold on;
ylabel('Proportion exact match','fontweight','bold');
xlabel('ICR','fontweight','bold');
yticks([0 0.25 0.5 0.75 1]);
xticks([1 2 3 4 5]);
xlim([0.5 5.5]); ylim([0 1.01]); axis square; box on;
set(gca,'xticklabels',{'1','2','3','4','monocular'});

%average
stimavg = mean(allsubj_ICR(:,2:end),2);
%95% CI
err1 = 1.96* std(allsubj_ICR(:,2:end),[],2)/sqrt(N);
%plot average
errorbar(allsubj_ICR(1:3,1),stimavg(1:3),err1(1:3),'ko','MarkerFaceColor',[0 0 0],'LineWidth',2,'MarkerSize',6);
errorbar(allsubj_ICR(4,1),stimavg(4),err1(4),'ko','MarkerFaceColor',[0 0 0],'LineWidth',2,'MarkerSize',6);

%plot individual
individualY = allsubj_ICR(:,2:end); individualY = individualY(:);
individualX = repmat(allsubj_ICR(:,1),N,1);
plot(individualX+rand(size(individualX,1),1)*0.1-0.05,individualY,'.','MarkerEdgeColor',[0.5 0.5 0.5]);

%% by stimulus type
allsubj_stim = [1 2 3 4]';

for s = 1:size(filenames,2)

    filename = filenames{s};
    load(filename);

    %add in ICR for each trial
    trial_ICR = round(max(dat2.stim(:,2:3),[],2)./(min(dat2.stim(:,2:3),[],2)));
    sdata = [dat2.stim, trial_ICR, dat2.resp(:,2)];

    for stim = 1:4 %stimulus patterns
        subdata = sdata(sdata(:,1)==stim ,:);
        propMatch = sum(subdata(:,5)==1)/size(subdata,1);
        if isempty(propMatch)
            keyboard
        end
        allsubj_stim(allsubj_stim(:,1) == stim, s+1) = propMatch;
    end
end

%average
stimavg = mean(allsubj_stim(:,2:end),2);
%95% CI
err1 = 1.96* std(allsubj_stim(:,2:end),[],2)/sqrt(N);
    

%make plot
figure(1);
subplot(1,3,2);
hold on;
%plot average
errorbar([1 2 3 4],stimavg,err1,'ko','MarkerFaceColor',[0 0 0],'LineWidth',2,'MarkerSize',6);
%plot individual
individualY = allsubj_stim(:,2:end); individualY = individualY(:);
individualX = repmat([1 2 3 4]',N,1);
plot(individualX+rand(size(individualX,1),1)*0.1-0.05,individualY,'.','MarkerEdgeColor',[0.5 0.5 0.5]);

ylabel('Proportion exact match','fontweight','bold');
xlabel('Stimulus type','fontweight','bold');
yticks([0 0.25 0.5 0.75 1]);
xticks([1 2 3 4]);
set(gca,'xticklabel',{'grating','noise','simple','complex'});
xlim([0.5 4.5]); ylim([0 1.01]); box on; axis square;



