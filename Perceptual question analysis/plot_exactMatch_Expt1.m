% Experiment 1 exact match plot,the proportion of exact match based on ICR or stimulus type

clear all;
close all;
%set figure size for 1x3 subplot
set(groot,'defaultfigureposition',[275,243,1158,420]);

%data paths
addpath('../data/Expt1'); 
filenames = dir(['../data/Expt1/*','mat']);
filenames = {filenames.name};
N = size(filenames,2); %number of subjects

%% exact match by ICR
%create a list of ratios
contrasts = [0 0.25 0.5 1];
ICR_ratio = [];
c = []; %store contrast combinations without repeat

for ii = 1:size(contrasts,2)
    for jj = ii:size(contrasts,2)
        c = [c; contrasts(ii), contrasts(jj)];
    end
end

ICR_ratio = c(:,2)./c(:,1);%higher contrast / lower contrast
ICR_ratio(isnan(ICR_ratio)) = []; %remove 0/0
ICR_ratio = unique(ICR_ratio);  %find the unique ICRs

%store subject results
%first column is the ICRs, columns 2:end will be each subject's data
allsubj_ICR = ICR_ratio;

for s = 1:size(filenames,2)
    
    filename = filenames{s};
    load(filename);
    
    %replace col 4 with ICR ratio
    sdata = [dat1.stim dat1.resp(:,2)];
    sdata(:,4) = max(sdata(:,2:3),[],2)./(min(sdata(:,2:3),[],2));
    
    for ICR = ICR_ratio'
        subdata = sdata(sdata(:,4)==ICR,:);
        propMatch = sum(subdata(:,5)==1)/size(subdata,1);

        %each column will be a subject, each row will be a different ICR
        allsubj_ICR(allsubj_ICR(:,1)==ICR,s+1) = propMatch;
    end
    
end


allsubj_ICR(4,1) = 5; %change Inf to 5 for plotting

%save('expt1_nomatch.mat',"allsubj_ICR");


%make plot
figure(1); 
subplot(1,3,1);
hold on;
ylabel('Proportion exact match','fontweight','bold');
xlabel('ICR','fontweight','bold');
yticks([0 0.25 0.5 0.75 1]);
xticks([1 2 3 4 5]);
xlim([0 6]); ylim([0 1.01]); axis square; box on;
set(gca,'xticklabels',{'1','2','3','4','Mono'});

%average
stimavg = mean(allsubj_ICR(:,2:end),2);
%95% CI
err1 = 1.96* std(allsubj_ICR(:,2:end),[],2)/sqrt(N);

%plot average
errorbar(allsubj_ICR(1:3,1),stimavg(1:3),err1(1:3),'k-o','MarkerFaceColor',[0 0 0],'LineWidth',2,'MarkerSize',6);
errorbar(allsubj_ICR(4,1),stimavg(4),err1(4),'ko','MarkerFaceColor',[0 0 0],'LineWidth',2,'MarkerSize',6);

%plot individual
individualY = allsubj_ICR(:,2:end); individualY = individualY(:);
individualX = repmat(allsubj_ICR(:,1),32,1);
plot(individualX+rand(size(individualX,1),1)*0.1-0.05,individualY,'.','MarkerEdgeColor',[0.5 0.5 0.5]);


%% exact match by stimulus type 
%store results
allsubj_stim = [1 2 3 4 5]';
for s = 1:size(filenames,2)

    filename = filenames{s};
    load(filename);

    %replace col 4 with ICR ratio
    sdata = [dat1.stim dat1.resp(:,2)];
    sdata(:,4) = max(sdata(:,2:3),[],2)./(min(sdata(:,2:3),[],2));

    for stimtype = 1:5  %for each stimulus type
        subdata = sdata(sdata(:,1)==stimtype ,:);
        propMatch = sum(subdata(:,5)==1)/size(subdata,1);
        allsubj_stim(allsubj_stim(:,1) ==stimtype,s+1) = propMatch;

    end

end


%average
stimavg = mean(allsubj_stim(:,2:end),2);
%95% CI
err1 = 1.96* std(allsubj_stim(:,2:end),[],2)/sqrt(N);

figure(1);
subplot(1,3,2);
hold on;
%plot average
errorbar([1 2 3 4 5],stimavg,err1,'ko','MarkerFaceColor',[0 0 0],'LineWidth',2,'MarkerSize',6);
%plot individual
individualY = allsubj_stim(:,2:end); individualY = individualY(:);
individualX = repmat([1 2 3 4 5]',32,1);
plot(individualX+rand(size(individualX,1),1)*0.2-0.1,individualY,'.','MarkerEdgeColor',[0.5 0.5 0.5]);

ylabel('Proportion exact match','fontweight','bold');
xlabel('Stimulus type','fontweight','bold');
yticks([0 0.25 0.5 0.75 1]);
xticks([1 2 3 4 5]);
set(gca,'xticklabel',{'grat','noise','histMatch','bandpass','1D noise'});
xlim([0.5 5.5]); ylim([0 1.01]); box on; axis square;

