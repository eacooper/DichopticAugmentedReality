%comparing experiment 1 vs experiment 2 types of perceptual effect

clear all;
close all;

%load data saved from plot_propEffect_Expt*.m script
load('expt1_effect.mat');
expt1 = allsubj_ICR;
load('expt2_effect.mat');
expt2 = allsubj_ICR;

txt = {'contrast','brightness','luster','rivalry','depth'};
for ICR = [4 Inf]
    if ICR ==4
        disp('---------this is ICR ==4---------')
    else
        disp('---------this is monocular trial---------')
    end

    for Q = 1:5
        expt1_ans = expt1(expt1(:,1)==ICR & expt1(:,3)==Q,:);
        expt2_ans = expt2(expt2(:,1)==ICR & expt2(:,3)==Q,:);
        if ICR ==4
            figure(1);
        else
            figure(2);
        end

        %plot the individual responses
        subplot(2,3,Q); hold on;
        plot(0.7*ones(1,size(expt1_ans,1))+0.3*rand([1 size(expt1_ans,1)]), expt1_ans(:,4),'r.');
        plot(1,mean(expt1_ans(:,4)),'r*','MarkerSize',10);
        plot(1.7*ones(1,size(expt2_ans,1))+0.3*rand([1 size(expt2_ans,1)])+0.5, expt2_ans(:,4),'b.');
        plot(2,mean(expt2_ans(:,4)),'b*','MarkerSize',10);
        yticks([0 0.25 0.5 0.75 1]);
        xlim([0 3]);
        xticks([1 2]); 
        set(gca,'XTickLabel',{'expt1','expt2'});

        [h,p,ci,stats] = ttest2(expt1_ans(:,4),expt2_ans(:,4),'Vartype','equal'); %t test
        diff_mean = meanEffectSize(expt1_ans(:,4),expt2_ans(:,4),'Effect',"cohen"); %cohen's d
        diff_mean = diff_mean{1,1};

        if ICR ==4
            title('ICR = 4',txt{Q});
        else
            title('mono trials',txt{Q});
        end
        
        if p<0.05
            disp([txt{Q},' t(61)= ',num2str(stats.tstat), '  p=', num2str(p),'**  ', ' d=', num2str(diff_mean)])
        else
            disp([txt{Q},' t(61)= ',num2str(stats.tstat), '  p=', num2str(p),'    ', ' d=', num2str(diff_mean)])
        end
        
    end
end

