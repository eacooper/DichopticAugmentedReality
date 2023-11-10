%glme for expt 2 exact match question

clear all;

%data paths
addpath('../data/Expt2_AR/'); 
filenames = dir(['../data/Expt2_AR/*','mat']);
filenames = {filenames.name};

%initiate data structures
glmdata = {};
subjid = {};
ICR = [];
mono = {};
stimtype = {};
Q1resp = [];

monotxt = {'nonmonocular','monocular'};
stimname = {'grating','noise','simple','complex'};
maxICR = 100;

%% process the data
for s = 1:size(filenames,2)
    
    filename = filenames{s};
    id = extractBetween(filename,'_','.mat'); id = id{:}; %get subj id
    load(filename);

    %add in ICR for each trial
    trial_ICR = round(max(dat2.stim(:,2:3),[],2)./(min(dat2.stim(:,2:3),[],2))); % compute ICR
    sdata = [dat2.stim, trial_ICR, dat2.resp(:,2)];
    ismono = sdata(:,2) == 0 | sdata(:,3) ==0; %whether or not a trial is a monocular trial
    
    for t = 1:size(sdata,1)
        trialdata = sdata(t,:);
        subjid = [subjid; ['S',id]];
        ratio = trialdata(4);
        if isinf(ratio)
            ratio = maxICR; %recode inf to maxICR
        end
        
        ICR = [ICR; round(ratio)];
        mono = [mono;monotxt{ismono(t)+1}];
        stimtype = [stimtype; stimname{trialdata(1)}];
        Q1resp = [Q1resp; trialdata(5)==1]; %exact match        
    end
end


T_full = table(subjid,stimtype,mono,ICR,Q1resp);
T = T_full;

T.mono = categorical(T.mono);
T.mono = reordercats(T.mono,{'nonmonocular','monocular'});
T.stimtype = categorical(T.stimtype);
T.stimtype = reordercats(T.stimtype,stimname);

T_ICR = T(T.ICR~=maxICR,:); %remove the monocular trials

%% mini model for guideline generation, Eq 3 in Discussion, ICR=1,2,4 only
glme_ICR1 = fitglme(T_ICR,'Q1resp ~ 1 + ICR + (1|subjid)',...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace', ...
    'DummyVarCoding','reference','verbose',0);
intercept = glme_ICR1.Coefficients{1,2};
slope = glme_ICR1.Coefficients{2,2};

P = [0.8, 0.9];  %proportion threshold
ICRval = (log(P./(1-P)) - intercept)./slope; %calculate the ICR value for the P threshold

%% main result
glme_main = fitglme(T_ICR,'Q1resp ~ 1 + ICR+stimtype + (1|subjid)',...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace', ...
    'DummyVarCoding','reference','verbose',0);

T_mono = T(T.ICR==maxICR | T.ICR==4,:); %compare the monocular trial
glme_mono = fitglme(T_mono,'Q1resp ~ 1 + mono + (1|subjid)',...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace', ...
    'DummyVarCoding','reference','verbose',0);

main_result = dataset2table(glme_main.Coefficients);
main_result.or = exp(main_result.Estimate); %odds ratio

main_result

mono_result = dataset2table(glme_mono.Coefficients);
mono_result.or = exp(mono_result.Estimate); %odds ratio
mono_result
