%For the dichoptic trials, what is the average response of not same 
%(i.e. number of effect present for a given trial out of the 5 effects)

clear all;
close all;

%% experiment 1
%data paths
addpath('../data/Expt1'); 
filenames = dir(['../data/Expt1/*','mat']);
filenames = {filenames.name};

NotSame = []; %not same reponses
NotSameHIGH = []; %not same reponses for high IOR trials only

for s = 1:size(filenames,2)
    filename = filenames{s};
    load(filename);
    
    %remove the stimtypes not considered 
    sdata = [dat1.stim dat1.resp];
    %convert to ratio
    sdata(:,4) = max(sdata(:,2:3),[],2)./(min(sdata(:,2:3),[],2));
    %recode found exact match as 'same' for the rest of the prompts
    sdata(isnan(sdata)) = 4;

    %grab the dichoptic trials
    dichoptic = sdata(sdata(:,2)~=sdata(:,3),:);
    %calculate the sum of non 4s numbers (not same) for prompt responses
    %(col 7 -11), per trial
    numNonSamePerTrial = sum(dichoptic(:,7:11)~=4,2);
    
    %grab the IOR = 4 trials
    highIOR = sdata(sdata(:,4) ==4,:);
    numNonSamePerTrialHIGH = sum(highIOR(:,7:11)~=4,2);
    
    %calculate the average number of nonsame response 
    SubjtrialNonSame = mean(numNonSamePerTrial);
    NotSame = [NotSame,SubjtrialNonSame];%store all subjects
    
    %calculate the average number of nonsame response 
    SubjtrialNonSameHIGH = mean(numNonSamePerTrialHIGH);
    NotSameHIGH = [NotSameHIGH,SubjtrialNonSameHIGH]; %store all subjects

end

disp('EXPERIMENT 1: mean, median, p, tstat')
%run simple t test to check if number of not same is greater than 1
[h,p,ci,stats] = ttest(NotSame,ones(size(NotSame)));
[mean(NotSame), median(NotSame),p,stats.tstat]

[h,p,ci,stats] = ttest(NotSameHIGH,ones(size(NotSameHIGH)));
[mean(NotSameHIGH), median(NotSameHIGH),p,stats.tstat]


%% experiment 2
clear all;
%data paths
addpath('../data/Expt2_AR/'); 
filenames = dir(['../data/Expt2_AR/*','mat']);
filenames = {filenames.name};

NonSame2 = [];
NonSame2HIGH = [];

for s = 1:size(filenames,2)
    filename = filenames{s};
    load(filename);
    
    %convert to ratio
    ICR = round(max(dat2.stim(:,2:3),[],2)./min(dat2.stim(:,2:3),[],2));
    sdata2 = [dat2.stim ICR dat2.resp];

    %convert to same if NaN
    sdata2(isnan(sdata2)) = 4;

    %grab the dichoptic trials
    dichoptic = sdata2(sdata2(:,2)~=sdata2(:,3),:);
    
    %col 7 - 11, calculate the sum of non 4s numbers per trial
    numNonSamePerTrial = sum(dichoptic(:,7:11)~=4,2);
    %calculate the average response of nonsame
    SubjtrialNonSame = mean(numNonSamePerTrial);
    NonSame2 = [NonSame2,SubjtrialNonSame];
    
    %IOR = 4
    highIOR = sdata2(sdata2(:,4) ==4,:);
    numNonSamePerTrialHIGH = sum(highIOR(:,7:11)~=4,2);
    %calculate the average response of nonsame
    SubjtrialNonSameHIGH = mean(numNonSamePerTrialHIGH);
    NonSame2HIGH = [NonSame2HIGH,SubjtrialNonSameHIGH];
end

disp('EXPERIMENT 2: mean, median, p, tstat')
%run simple t test to check if number of not same is greater than 1
[h,p,ci,stats] = ttest(NonSame2,ones(size(NonSame2)));
[mean(NonSame2), median(NonSame2),p,stats.tstat]

[h,p,ci,stats] = ttest(NonSame2HIGH,ones(size(NonSame2HIGH)));
[mean(NonSame2HIGH), median(NonSame2HIGH),p,stats.tstat]