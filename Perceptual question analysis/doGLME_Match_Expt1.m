%glme for expt 1 exact match question

clear all;

%data paths
addpath('../data/Expt1'); 
filenames = dir(['../data/Expt1/*','mat']);
filenames = {filenames.name};

%initiate data structures
glmdata = {};
subjid = {};
ICR = [];
mono = {};
stimtype = {};
Q1resp = [];

monotxt = {'nonmonocular','monocular'};
stimname = {'grating','noise','hist match','bandpass','broadband'};
maxICR = 100;

%% process the data
for s = 1:size(filenames,2)
    
    filename = filenames{s};
    id = extractBetween(filename,'_','.mat'); id = id{:}; %get subj id
    load(filename);

    sdata = [dat1.stim dat1.resp(:,2)];
    sdata(:,4) = max(sdata(:,2:3),[],2)./(min(sdata(:,2:3),[],2)); % compute ICR
    ismono = sdata(:,2) == 0 | sdata(:,3) ==0; %whether or not a trial is a monocular trial
    
    for t = 1:size(sdata,1)
        trialdata = sdata(t,:);
        subjid = [subjid; ['S',id]];
        %change ICR to ratio
        ratio = trialdata(4);
        if isinf(ratio)
            ratio = maxICR;
        end
        
        ICR = [ICR; ratio];
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

%% models
%main model, with stim type and ICR (exclude Monocular)
glme_main = fitglme(T_ICR,'Q1resp ~ 1 + stimtype+ICR + (1|subjid)',...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace', ...
    'DummyVarCoding','reference','verbose',0);


%compare monocular with ICR = 4
T_mono = T(T.ICR==maxICR | T.ICR==4,:);

glme_mono = fitglme(T_mono,'Q1resp ~ 1 + mono + (1|subjid)',...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace', ...
    'DummyVarCoding','reference','verbose',0);


main_result = dataset2table(glme_main.Coefficients);
main_result.or = exp(main_result.Estimate); %odds ratio

main_result

mono_result = dataset2table(glme_mono.Coefficients);
mono_result.or = exp(mono_result.Estimate); %odds ratio
mono_result

