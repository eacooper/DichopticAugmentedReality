%check what's the proportion of time that ppl responded with the 'unsure'
%option

clear all;
close all;
%% experiment 1
%data paths
addpath('../data/Expt1'); 
filenames = dir(['../data/Expt1/*','mat']);
filenames = {filenames.name};

%store results
unsure = [];
for s = 1:size(filenames,2)
    filename = filenames{s};
    load(filename);

    %col 7-11 are responses
    sdata = [dat1.stim dat1.resp];
    trialUnsure = sum(sdata(:,7:11)==5,'all'); %total number of unsure response
    unsure = [unsure,trialUnsure];

end
%calculate the percent of responses that were unsure for each subject
unsure = unsure*100/(size(sdata,1)*5);
%find the mean,median and std
disp('*****EXPERIMENT 1*****')
disp('overall mean, median, and standard error of unsure response across subjects')
[mean(unsure), median(unsure),std(unsure)/sqrt(size(sdata,1))]


% check the number of unsure responses for the first question vs the last
% question
unsureQ = [];
for s = 1:size(filenames,2)
    filename = filenames{s};
    load(filename);
    sdata = [dat1.stim dat1.resp];

    trialUnsure = 100*sum(sdata(:,7:11)==5,1)/size(sdata,1);
    
    unsureQ = [unsureQ;trialUnsure];
end
%average across people
disp('average and std for the unsure response for each prompt')
mean(unsureQ,1)
std(unsureQ,[],1)


%% experiment 2
%data paths
addpath('../data/Expt2_AR/'); 
filenames = dir(['../data/Expt2_AR/*','mat']);
filenames = {filenames.name};

%store result
unsure2 = [];
for s = 1:size(filenames,2)
    filename = filenames{s};
    load(filename);
    sdata2 = [dat2.stim dat2.resp];
    trialUnsure2 = sum(sdata2(:,6:10)==5,'all'); %total number of unsure responses
    unsure2 = [unsure2,trialUnsure2];
end
%calculate the percent of responses that were unsure for each subject
unsure2 = unsure2*100/(size(sdata2,1)*5);
disp('*****EXPERIMENT 2*****')
disp('overall mean, median, and standard error of unsure response across subjects')
[mean(unsure2), median(unsure2),std(unsure2)/sqrt(size(sdata2,1))]

% check the number of unsure responses for the first question vs the last
% question
unsureQ2 = [];
for s = 1:size(filenames,2)
    filename = filenames{s};
    load(filename);
    sdata2 = [dat2.stim dat2.resp];
    trialUnsure2 = 100*sum(sdata2(:,6:10)==5,1)/size(sdata2,1);
    unsureQ2 = [unsureQ2;trialUnsure2];
end
%average across people
disp('average and std for the unsure response for each prompt')
mean(unsureQ2,1)
std(unsureQ2,[],1)