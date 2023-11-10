% plot the 5 effects, Figure 10

clear all;
close all;

set(groot,'defaultfigureposition',[275,243,1158,597]);
%data paths
addpath('../data/Expt2_AR/'); 
filenames = dir(['../data/Expt2_AR/*','mat']);
filenames = {filenames.name};

%% effect by ICR
%create a list of ratios
luminance = [0 0.13 0.25 0.5];
ICR_ratio = [];

c = [];
for ii = 1:size(luminance,2)
    for jj = ii:size(luminance,2)
        c = [c; luminance(ii), luminance(jj)];
    end
end
ICR_ratio = c(:,2)./c(:,1);%higher luminance / lower luminance
ICR_ratio(isnan(ICR_ratio)) = []; %remove 0/0

ICR_ratio = round(ICR_ratio); %round it
ICR_ratio = unique(ICR_ratio);  %find the unique ICRs

allsubj_ICR = [];

for s = 1:size(filenames,2)
    filename = filenames{s};
    load(filename);
    
    %add in ICR for each trial
    trial_ICR = round(max(dat2.stim(:,2:3),[],2)./(min(dat2.stim(:,2:3),[],2)));
    sdata = [dat2.stim, trial_ICR, dat2.resp(:,2:end)];
    
    %find the trials with exact match and mark the follow up response as 4s for SAME
    sdata(isnan(sdata)) = 4;
    
    %transform responses to same(1) vs. other responses (0)
    sdata(:,6:10) = sdata(:,6:10)==4;

    for ICR = ICR_ratio'
        ICRdata = sdata(sdata(:,4) == ICR,:);
        
        for Q=1:5
            %find the proportion of trials where the two stimuli look
            %different in some aspect
            propDiff = sum(ICRdata(:,5+Q)==0)/size(ICRdata,1);
            allsubj_ICR = [allsubj_ICR; ICR,s,Q,propDiff];
        end
    end
end

save('expt2_effect.mat',"allsubj_ICR");

%make plots, 1 for monocular and one for nonmoncular
figure(1);
subplot(2,2,1); hold on;
markers = {'s','d','*','^','pentagram'};

colors = [0 0.4470 0.7410; 0.8500 0.3250 0.0980; 0.9290 0.6940 0.1250; 0.4940 0.1840 0.5560;...
    0.4660 0.6740 0.1880; 0.3010 0.7450 0.9330; 0.6350 0.0780 0.1840];

for Q = 1:5
    subdata = allsubj_ICR(allsubj_ICR(:,3)==Q,:);
    subdata = sortrows(subdata);
    
    trialtype = []; %reorg data for each question, each column will be a different ICR
    
    for ICR = ICR_ratio'
        trialtype(:,end+1) = subdata(subdata(:,1)==ICR,4);
    end

    %find the average and 95%CI
    avgdata = mean(trialtype,1); %average across subjects
    err1 = 1.96*std(trialtype,[],1)/sqrt(size(trialtype,1));
    
    errorbar([1 2 4]+Q*0.01-0.02,avgdata(1:3),err1(1:3),[markers{Q}],'Color',colors(Q,:),'MarkerFaceColor',colors(Q,:),'MarkerEdgeColor',colors(Q,:),'MarkerSize',10);
    errorbar(5+Q*0.01-0.02,avgdata(4),err1(4),[markers{Q}],'Color',colors(Q,:),'MarkerFaceColor',colors(Q,:),'MarkerEdgeColor',colors(Q,:),'MarkerSize',10,'HandleVisibility','Off');
end
legend('contrast','brightness','luster','rivalry','depth','Location','Northwest');
xlabel('ICR','fontweight','bold'); ylabel('Prop effect present','fontweight','bold');
xticks([1 2 3 4 5]); 
set(gca,'xticklabels',{'1','2','','4','mono'});
ylim([0 1]); yticks([0 0.25 0.5 0.75 1]);
xlim([0 6]);box on; axis square;


%% by stimulus conditions
clear all;
%data paths
addpath('../data/Expt2_AR/'); 
filenames = dir(['../data/Expt2_AR/*','mat']);
filenames = {filenames.name};

stimname = {'grating','noise','simple','complex','wild card'};
effects = {'contrast','brightness','luster','rivalry','depth'};
c = {'b','b','k','r','r','r','m','k-'};


%average proportion for each stimulus
stimprop= [];

for stim_ind = 1:4

    subjprop = [];  %each stim for all subject

    for f = 1:size(filenames,2)

        % load file
        load(filenames{f});

        %add in ICR for each trial
        trial_ICR = round(max(dat2.stim(:,2:3),[],2)./(min(dat2.stim(:,2:3),[],2)));
        sdata = [dat2.stim, trial_ICR, dat2.resp(:,2:end)];

        % grab the dichoptic trials for this stimulus type
        subdata = sdata(sdata(:,1)==stim_ind & sdata(:,2) ~=sdata(:,3),:);
        subdata(isnan(subdata)) = 4; %recode found exact match as SAME

        %recode the responses as same(0) vs not same(1)
        subdata(:,6:10) = subdata(:,6:10)~=4;
        %calculate proportion effect present out of all trials
        propCalc = sum(subdata(:,6:10),1)./size(subdata,1);

        %store all subject's proportion, col = question, row = subject
        subjprop = [subjprop; propCalc];

    end

    %average across all subjects
    stimprop = [stimprop; mean(subjprop,1)];
end

figure(1)
subplot(2,2,2);
imagesc(stimprop); caxis([0 1]); colormap copper;
a = colorbar; a.Label.String = 'Proportion effect present';
a.Ticks = [0 0.25 0.5 0.75 1];
xticks([1:1:5]);
yticks([1:1:6]); xtickangle(45);
set(gca,'xticklabels',effects,'yticklabels',stimname);
axis square; title('all dichoptic')
disp(stimprop)

