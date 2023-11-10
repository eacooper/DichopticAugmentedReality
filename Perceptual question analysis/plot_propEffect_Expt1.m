% plot the 5 effects, Figure 7 

clear all;
close all;

%set figure size
set(groot,'defaultfigureposition',[275,243,1158,597]);

%data paths
addpath('../data/Expt1'); 
filenames = dir(['../data/Expt1/*','mat']);
filenames = {filenames.name};


%% effect by ICR
%create a list of ratios
contrasts = [0 0.25 0.5 1];
ICR_ratio = [];
c = [];
for ii = 1:size(contrasts,2)
    for jj = ii:size(contrasts,2)
        c = [c; contrasts(ii), contrasts(jj)];
    end
end
ICR_ratio = c(:,2)./c(:,1);%larger contrast / smaller contrast
ICR_ratio(isnan(ICR_ratio)) = []; %remove 0/0
ICR_ratio = unique(ICR_ratio);  %find the unique ICRs

allsubj_ICR = [];

%% process data
for s = 1:size(filenames,2)
    
    filename = filenames{s};
    load(filename);
    
    %replace col 4 with ICR ratio
    sdata = [dat1.stim dat1.resp(:,2:end)];
    sdata(:,4) = max(sdata(:,2:3),[],2)./(min(sdata(:,2:3),[],2));

    %find the trials with exact match and mark the follow up response as 4s for SAME
    sdata(isnan(sdata)) = 4;
    
    %transform responses to same(1) vs. other responses (0)
    sdata(:,6:end) = sdata(:,6:end)==4;
    
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

save('expt1_effect.mat',"allsubj_ICR");

%% make plots
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
    %find the average and 95% CI
    avgdata = mean(trialtype,1); %average across subjects
    err1 = 1.96*std(trialtype,[],1)/sqrt(size(trialtype,1));
    
    %ICR = 1 2 4
    errorbar([1 2 4]+Q*0.01-0.02,avgdata(1:3),err1(1:3),[markers{Q},'-'],'Color',colors(Q,:),'MarkerFaceColor',colors(Q,:),'MarkerEdgeColor',colors(Q,:),'MarkerSize',10);
    %monocular
    errorbar(5+Q*0.01-0.02, avgdata(4),err1(4),[markers{Q},'-'],'Color',colors(Q,:),'MarkerFaceColor',colors(Q,:),'MarkerEdgeColor',colors(Q,:),'MarkerSize',10,'HandleVisibility','off');
    
end

legend('contrast','brightness','luster','rivalry','depth','Location','Northwest');
xlabel('ICR','fontweight','bold'); ylabel('Proportion of time present','fontweight','bold');
xticks([1 2 3 4 5]);
set(gca,'xticklabels',{'1','2','','4','mono'});
yticks([0 0.25 0.5 0.75 1]);
ylim([0 1]);
xlim([0 6]);box on; axis square;


%% for each stim, what is probability of specific effect when there's no match(row=stim, col=effects)
clear all;
%data paths
addpath('../data/Expt1'); 
filenames = dir(['../data/Expt1/*','mat']);
filenames = {filenames.name};

stimname = {'grating','noise','histeq','bandpass','broadband'};
effects = {'contrast','brightness','luster','rivalry','depth'};

%average proportion for each stimulus
stimprop= [];

for stim_ind = 1:5 
    subjprop = []; %each stim for all subject
    
    for f = 1:size(filenames,2)
        % load file
        load(filenames{f});

        % grab the DICHOPTIC trials for this stimulus type
        sdata = [dat1.stim dat1.resp(:,2:end)];
        subdata = sdata(sdata(:,1)==stim_ind & sdata(:,2)~= sdata(:,3),:);
        subdata(isnan(subdata)) = 4; %recode found exact match as SAME

        %recode the responses as same vs nonsame
        subdata(:,6:end) = subdata(:,6:end)~=4;
        %calculate the proportion of no match
        propCalc = sum(subdata(:,6:10),1)./size(subdata,1);
        %store all subject's proportion
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






