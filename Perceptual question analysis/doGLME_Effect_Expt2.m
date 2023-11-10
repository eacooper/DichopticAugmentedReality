% do a GLME for each effect separately
clear all;

%data paths
addpath('../data/Expt2_AR/'); 
filenames = dir(['../data/Expt2_AR/*','mat']);
filenames = {filenames.name};

subjid = {};
Qtype = {};
ICR = [];
mono = {};
stimtype = {};
Qresp = [];

monotxt = {'nonmonocular','monocular'};
stimname = {'grating','noise','simple','complex'};
prompts = {'contrast','brightness','luster','rivalry','depth'};

maxICR = 100;

%% process the data
for s = 1:size(filenames,2)
    
    filename = filenames{s};
    id = extractBetween(filename,'_','.mat'); id = id{:};
    
    load(filename);

    %add in ICR for each trial
    trial_ICR = round(max(dat2.stim(:,2:3),[],2)./(min(dat2.stim(:,2:3),[],2))); % compute ICR
    sdata = [dat2.stim, trial_ICR, dat2.resp(:,2:end)];
    ismono = sdata(:,2) == 0 | sdata(:,3) ==0; %whether or not a trial is a monocular trial
    
    sdata_recode = sdata;
    %recode found exact match as same
    sdata_recode(isnan(sdata_recode)) = 4;
    %recode as same (0, not present) vs other (1, present)
    sdata_recode(:,6:10) = sdata_recode(:,6:10)~=4;
        
    %build lists of variables for table
    for Q = 1:5 %for each effect question
        for t = 1:size(sdata_recode,1)
            trialdata = sdata_recode(t,:);
            subjid = [subjid; ['S',id]];
            Qtype = [Qtype; prompts{Q}];
            ratio = trialdata(4);
            if isinf(ratio)
                ratio = maxICR;
            end
            
            ICR = [ICR; ratio];
            mono = [mono;monotxt{ismono(t)+1}];
            stimtype = [stimtype; stimname{trialdata(1)}];
            Qresp = [Qresp; trialdata(Q+5)];
        end
    end
end


T_full = table(subjid,stimtype,mono,ICR,Qtype,Qresp);


T = T_full;
T.mono = categorical(T.mono);
T.mono = reordercats(T.mono,{'nonmonocular','monocular'});
T.stimtype = categorical(T.stimtype);
T.stimtype = reordercats(T.stimtype,stimname);

%% main ICR model
T2 = T(T.ICR~=maxICR,:); %remove the monocular trials
%separate data for each effect type
T_contrast = T2(strcmp(T2.Qtype,'contrast'),:);
T_brightness = T2(strcmp(T2.Qtype,'brightness'),:);
T_luster = T2(strcmp(T2.Qtype,'luster'),:);
T_rivalry = T2(strcmp(T2.Qtype,'rivalry'),:);
T_depth = T2(strcmp(T2.Qtype,'depth'),:);

glme_contrast = fitglme(T_contrast,'Qresp ~ 1 + ICR + (1|subjid)',...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace', ...
    'DummyVarCoding','reference','verbose',0);
glme_brightness = fitglme(T_brightness,'Qresp ~ 1 + ICR + (1|subjid)',...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace', ...
    'DummyVarCoding','reference','verbose',0);
glme_luster = fitglme(T_luster,'Qresp ~ 1 + ICR + (1|subjid)',...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace', ...
    'DummyVarCoding','reference','verbose',0);
glme_rivalry = fitglme(T_rivalry,'Qresp ~ 1 + ICR + (1|subjid)',...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace', ...
    'DummyVarCoding','reference','verbose',0);
glme_depth = fitglme(T_depth,'Qresp ~ 1 + ICR + (1|subjid)',...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace', ...
    'DummyVarCoding','reference','verbose',0);


disp('MAIN ICR')
%get the odds ratio
coeff = [glme_contrast.Coefficients{2,2},...
    glme_brightness.Coefficients{2,2},...
    glme_luster.Coefficients{2,2},...
    glme_rivalry.Coefficients{2,2},...
    glme_depth.Coefficients{2,2}];

ORs = exp(coeff);

coeff
ORs


%% monocular vs ICR = 4
T3 = T(T.ICR==maxICR | T.ICR==4,:);

%subset of data for each regression
T_contrastM = T3(strcmp(T3.Qtype,'contrast'),:);
T_brightnessM = T3(strcmp(T3.Qtype,'brightness'),:);
T_lusterM = T3(strcmp(T3.Qtype,'luster'),:);
T_rivalryM = T3(strcmp(T3.Qtype,'rivalry'),:);
T_depthM = T3(strcmp(T3.Qtype,'depth'),:);


glme_contrastM = fitglme(T_contrastM,'Qresp ~ 1 + mono + (1|subjid)',...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace', ...
    'DummyVarCoding','reference','verbose',0);
glme_brightnessM = fitglme(T_brightnessM,'Qresp ~ 1 + mono + (1|subjid)',...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace', ...
    'DummyVarCoding','reference','verbose',0);
glme_lusterM = fitglme(T_lusterM,'Qresp ~ 1 + mono + (1|subjid)',...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace', ...
    'DummyVarCoding','reference','verbose',0);
glme_rivalryM = fitglme(T_rivalryM,'Qresp ~ 1 + mono + (1|subjid)',...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace', ...
    'DummyVarCoding','reference','verbose',0);
glme_depthM = fitglme(T_depthM,'Qresp ~ 1 + mono + (1|subjid)',...
    'Distribution','Binomial','Link','logit','FitMethod','Laplace', ...
    'DummyVarCoding','reference','verbose',0);

disp('mono vs nonmono')

%calculate odds ratios
coeff2 = [glme_contrastM.Coefficients{2,2},...
    glme_brightnessM.Coefficients{2,2},...
    glme_lusterM.Coefficients{2,2},...
    glme_rivalryM.Coefficients{2,2},...
    glme_depthM.Coefficients{2,2}];

ORs2 = exp(coeff2);

coeff2
ORs2