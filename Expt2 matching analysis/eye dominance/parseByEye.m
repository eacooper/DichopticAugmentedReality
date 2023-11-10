%separate data by left eye dominant or right eye dominant based on matching
%result, LE dominant means that when left eye was shown higher contrast stim,
%resulted in more winne-take-all than when right eye was shown higher
%contrast stim

clear all;
close all;
addpath('../../data/Expt2_AR');
%set(groot,'defaultfigureposition',[275,243,1158,420]);
filenames = dir(['../../data/Expt2_AR/BAR_*','mat']);
filenames = {filenames.name};
N = size(filenames,2); %number of subjects

%initiate to store result
eyedom = [];

for f = 1:N %for each subject
    
    filename = filenames{f};
    load(filename);
    
    data = [dat2.stim dat2.resp(:,1)]; %4 columns
    
            
    for stim = 1:5
        if stim <5
            subdata = data(data(:,1)==stim,:); %pick the one stimulus pattern
        else
            subdata = data; %consider all stimulus patterns
        end

        %trials where LE > RE contrast
        LE_hc= subdata(subdata(:,2)>subdata(:,3),:);
        LE_hc(:,5) = max(LE_hc(:,2:3),[],2); %add the winner-take-all prediction
        
        %RE>LE
        RE_hc= subdata(subdata(:,2)<subdata(:,3),:);
        RE_hc(:,5) = max(RE_hc(:,2:3),[],2); %add the winner-take-all prediction
        
        %calculate RMSE, lower the better fit
        LE_RMSE = sqrt(mean((LE_hc(:,5)-LE_hc(:,4).^2)));
        RE_RMSE = sqrt(mean((RE_hc(:,5)-RE_hc(:,4).^2)));
        disp([LE_RMSE,RE_RMSE])
        isLEdom = LE_RMSE<RE_RMSE;
        eyedom(f,stim) = isLEdom;
    end
end

save('ARexpt_eyedom','eyedom');

%check to see how similar the result is based on which stimulus pattern was
%selected
sum(eyedom,1)
corr(eyedom)
