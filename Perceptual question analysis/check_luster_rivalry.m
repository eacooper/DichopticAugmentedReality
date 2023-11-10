%check whether the dichoptic stimulus was chosen for luster/rivalry?

clear all;
close all;

%% experiment 1

%data paths
addpath('../data/Expt1'); 
filenames = dir(['../data/Expt1/*','mat']);
filenames = {filenames.name};

N = size(filenames,2); %number of subjects

question = {'contrast','brightness','luster','rivalry','depth'};
stimname = {'grating','noise','histeq','bandpass','broadband'};

%to store results for the proportion of time that dichoptic, non dichoptic
%was selected (Nx4x5 matrix, row = subj, col = response type, d3 = question)
Q_checktype = NaN(size(filenames,2),4,5);

%for each subject, calculate the proportion of response
for s = 1:size(filenames,2)

    filename = filenames{s};
    load(filename);

    %Check the saved info for whether top or bottom was the reference
    if strcmp(dat1.reference,'TOP')
        refisTop = 1;
    elseif strcmp(dat1.reference,'BOTTOM')
        refisTop = 0;
    end

    sdata = [dat1.stim dat1.resp];
    %col 7 - 11 are responses 
    %recode top/bottom into dichoptic(1) or nondichoptic(2)
    if refisTop == 0 
        temp = sdata(:,7:11);
        temp(sdata(:,7:11)==1) = 2; %top is nondichoptic        
        temp(sdata(:,7:11)==2) = 1; %bottom is reference
        sdata(:,7:11) = temp;
    else %keep the way it is
    end

    %calculate the proportion of dichoptic non dichoptic response for this subject
    sdata = sdata(sdata(:,6)==2,:);  %out of the no match trials
    for q = 1:5 %for each of the effect question
        dich = sum(sdata(:,6+q) == 1)/size(sdata,1);
        ndich = sum(sdata(:,6+q) == 2)/size(sdata,1);
        same = sum(sdata(:,6+q) == 4)/size(sdata,1);
        unsure = sum(sdata(:,6+q) == 5)/size(sdata,1);

        Q_checktype(s,:,q) = [dich, ndich, same, unsure];

    end


end

%find the average across subjects, and the 95% CI
avg = mean(Q_checktype,1);
errs = 1.96*std(Q_checktype,1)/sqrt(size(Q_checktype,1));

%plot each question
for qind = 1:5 
    effectsplot_line = figure(1);
    subplot(2,5,qind);
    hold on;

    %four response types
    e1 = errorbar([1 2 3 4], avg(:,:,qind),errs(:,:,qind),'o','LineWidth',2,'MarkerSize',4);
    %e1.Color = [0.4940 0.1840 0.5560];
    ylim([0 1]);
    xlim([0 5]);
    ylabel('prop chosen');

    title(question{qind});
    xticks([1 2 3 4]);
    set(gca,'xticklabels',{'Dichoptic','Non-dich','Same','Unsure'});
    xtickangle(40)
    axis square;
    box on;

    %ratio for when top or bottom stimulus was chosen, it was the dichoptic one
    dd = avg(:,:,qind);
    ratio = dd(1)/(dd(1)+dd(2));
    ratio = round(ratio,2);
    disp([question{qind},' dichoptic chosen ', num2str(ratio)])

end

%% experiment AR
clear all;
%data paths
addpath('../data/Expt2_AR/'); 
filenames = dir(['../data/Expt2_AR/*','mat']);
filenames = {filenames.name};

N = size(filenames,2); %number of subjects

question = {'contrast','brightness','luster','rivalry','depth'};

%store results
Q_checktype = NaN(size(filenames,2),4,5);

%for each subject, calculate the proportion of response
for s = 1:size(filenames,2)

    filename = filenames{s};
    load(filename);
        
    %Check the saved info for whether top or bottom was the reference
    if strcmp(dat2.reference,'TOP')
        refisTop = 1;
    elseif strcmp(dat2.reference,'BOTTOM')
        refisTop = 0;
    end

    sdata = [dat2.stim dat2.resp];
    %col 6 - 10 are responses
    %recode top/bottom into dichoptic(1) or nondichoptic(2)
    if refisTop == 0
        temp = sdata(:,6:10);
        temp(sdata(:,6:10)==1) = 2; %top is nondichoptic
        temp(sdata(:,6:10)==2) = 1; %bottom is reference
        sdata(:,6:10) = temp;
    else %keep the way it is
    end

    
    %calculate the proportion of times for each of the response types
    sdata = sdata(sdata(:,5)==2,:); %out of the no match trials 
    for q = 1:5 %for each of the effect question
        dich = sum(sdata(:,5+q) == 1)/size(sdata,1);
        ndich = sum(sdata(:,5+q) == 2)/size(sdata,1);
        same = sum(sdata(:,5+q) == 4)/size(sdata,1);
        unsure = sum(sdata(:,5+q) == 5)/size(sdata,1);

        Q_checktype(s,:,q) = [dich, ndich, same, unsure];

    end


end

%find the average across subjects, and the 95% CI
avg = mean(Q_checktype,1);
errs = 1.96*std(Q_checktype,1)/sqrt(size(Q_checktype,1));

%plot each response
for qind = 1:5
    subplot(2,5,qind+5); 
    hold on;

    %four response types
    e1 = errorbar([1 2 3 4], avg(:,:,qind),errs(:,:,qind),'o','LineWidth',2,'MarkerSize',4);
    e1.Color = [0.4940 0.1840 0.5560]; %use a different color
    ylim([0 1]);
    xlim([0 5]);
    ylabel('prop chosen');

    title(question{qind});
    xticks([1 2 3 4]);
    set(gca,'xticklabels',{'Dichoptic','Non-dich','Same','Unsure'});
    xtickangle(40)
    axis square;
    box on;

    %ratio for when top or bottom stimulus was chosen, it was the dichoptic one
    dd = avg(:,:,qind);
    ratio = dd(1)/(dd(1)+dd(2));
    ratio = round(ratio,2);
    disp([question{qind},' dichoptic chosen ', num2str(ratio)])

end
