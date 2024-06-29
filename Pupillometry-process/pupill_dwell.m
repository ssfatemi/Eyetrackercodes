%% Dwell time percentage

%%
clc
clear all
load('pupillary_data_55.mat')

%% dtp level 1-6
clc
distraction = nan(55,24);
dtp = nan(55,24);
dtp_e = nan(55,24);
dtp_n = nan(55,24);
dtp2 = nan(55,24);
int = 5000;
strt = 1001;
for s= 1:55
    for i = 1:24
    total_time = sum(~isnan(emotional_aoi_indexes(s,i,1:6000))) + sum(~isnan(neutral_aoi_indexes(s,i,1:6000)));
    aoi_time = sum(emotional_aoi_indexes(s,i,strt:6000)==1) + sum(neutral_aoi_indexes(s,i,1000:6000)==1);
    aoi_time_e = sum(emotional_aoi_indexes(s,i,strt:6000)==1);
    aoi_time_n = sum(neutral_aoi_indexes(s,i,strt:6000)==1);
    
    not_aoi_time = sum(emotional_aoi_indexes(s,i,strt:6000)==0) + sum(neutral_aoi_indexes(s,i,strt:6000)==0);
    distraction_percent = round(not_aoi_time/6000,2);
    dwell_time_percent =  round(aoi_time/int,3);
    dwell_time_percent_e =  round(aoi_time_e/int,3);
    dwell_time_percent_n =  round(aoi_time_n/int,3);
    dwell_time_percent2 =  round((aoi_time+1000)/6000,2);
    distraction(s,i)= distraction_percent;
        if dwell_time_percent ==0
           dtp(s,i)= nan;
           dtp2(s,i)= nan;
        else
           dtp(s,i)= dwell_time_percent;
           dtp2(s,i)= dwell_time_percent2;
        end
        
        if dwell_time_percent_e ==0
           dtp_e(s,i)= nan;
        else
           dtp_e(s,i)= dwell_time_percent_e;
        end
        
        if dwell_time_percent_n ==0
           dtp_n(s,i)= nan;
        else
           dtp_n(s,i)= dwell_time_percent_n;
        end
        
    end
end
msg = 'distraction extracted'
% Distraction Samples plot
for i=1:100
    d_plot(i) = ((sum(sum(dtp>(i/100))))/1056)*100;
end
plot(d_plot)
%% violin
violinplot([reshape(dtp_e,1320,1), reshape(dtp_n,1320,1)]);
figure
violinplot([nanmean(dtp_e,1)', nanmean(dtp_n,1)']);
figure
violinplot([nanmean(dtp_e,2), nanmean(dtp_n,2)]);
[h_i_v,p_i_v,ci,stats_i_v] = ttest(nanmean(dtp_e,2), nanmean(dtp_n,2))
%hist(reshape(dtp,1,1320),20)
%%
violin_dtp_e = round(nanmean(dtp_e,2)*int/1000,2);
violin_dtp_n = round(nanmean(dtp_n,2)*int/1000,2);
violinplot(round([nanmean(dtp_e,2), nanmean(dtp_n,2)]*int/1000,2));
%% plot without distractions
clc
males = logical(subject_info(1:55,3));
plot_n=0;

plots_cell = cell(3,1);
figure
for dtp_i = [ 0.2 ];
    plot_n=plot_n+1;
    filtered_all_neutral_pupil_3d = nan(size(all_neutral_pupil_3d));
    filtered_all_emotional_pupil_3d = nan(size(all_emotional_pupil_3d));
    removed_subjects = [5,18,20,22,24,32,34,35,36,37,43];
%     for s = setdiff([1:sample_size],removed_subjects);
    for s= 1:55
        for i = 1:24
            if dtp(s,i)>= dtp_i && dtp(s,i) <= 0.999
                if ~isnan(all_neutral_pupil_3d(s,i))
                   filtered_all_neutral_pupil_3d(s,i,1:6000) = all_neutral_pupil_3d(s,i,1:6000);
                else
                   filtered_all_emotional_pupil_3d(s,i,1:6000) = all_emotional_pupil_3d(s,i,1:6000);
                end
            end
        end
    end

    i = 10;
    a = squeeze(nanmean(filtered_all_emotional_pupil_3d,2));
    b = squeeze(nanmean(filtered_all_neutral_pupil_3d,2));

    % figure
    % plot(a(i,:)','color',[1 0 0])
    % hold on
    % plot(b(i,:)','color',[0 0 1])
    % legend Emotional Neutral

    e = squeeze(nanmean(a(:,:),1));
    n = squeeze(nanmean(b(:,:),1));
    normal_e = round(sqrt(e)/10,3);
    normal_n = round(sqrt(n)/10,3);
    
    ploting_e = normal_e';
    ploting_n = normal_n';
    
    
    e_q_dif = quantile(round(sqrt(a)/10,3),0.75)- quantile(round(sqrt(a)/10,3),0.25);
    n_q_dif = quantile(round(sqrt(b)/10,3),0.75)- quantile(round(sqrt(b)/10,3),0.25);
    
%     ploting_shadow_e = nanstd(round(sqrt(a)/10,3));
%     ploting_shadow_n = nanstd(round(sqrt(b)/10,3));
    ploting_shadow_e = e_q_dif;
    ploting_shadow_n = n_q_dif;
    
    subplot(3,3,plot_n)
    plot(normal_e','color',[1 0 0]) %Emotional
    hold on
%   errorbar
%     y = ploting_e(;
%    err = std(ploting_shadow_e);
%     errorbar(y,err)
    
    hold on
    plot(normal_n','color',[0 0 1]) %Neutral
    hold on
%     plot(ploting_shadow_n','color',[0 0 0.3])
    legend Emotional Neutral
    plots_cell{plot_n} = [normal_e', normal_n'];   
    % Anova test mean
    a_2_6 = squeeze(nanmean(filtered_all_emotional_pupil_3d(:,:,2500:6000),2));
    b_2_6 = squeeze(nanmean(filtered_all_neutral_pupil_3d(:,:,2500:6000),2));
    e_2_6 = squeeze(nanmean(a_2_6,2));
    n_2_6 = squeeze(nanmean(b_2_6,2));
    normal_e = round(sqrt(e_2_6)/10,2);
    normal_n = round(sqrt(n_2_6)/10,2);
    %
    pupil_data_2_6 = [normal_e, normal_n];

    %anova1(pupil_data_2_6)

    % ttest
    [h,p,ci,stats] = ttest(normal_e, normal_n);
    l= {'Emotional Stimuli','Neutral Stimuli'};
    % figure
    % boxplot([normal_e, normal_n],'Labels',l)

    % violin plot
    % figure
    %violinplot([normal_e, normal_n]);
    title(['dtp >=' num2str(dtp_i) '--p:' num2str(round(p,3)) '--t:' num2str(round(stats.tstat,2))])
    disp(['dp>=' num2str(dtp_i) '--h:' num2str(round(h,3)) '--p:' num2str(round(p,3)) '--t:' num2str(round(stats.tstat,2))...
        ,'--Me:' num2str(round(nanmedian(normal_e),2)) '±' num2str(round(nanstd(normal_e),2))...
        ,'--Mn:' num2str(round(nanmedian(normal_n),2)) '±' num2str(round(nanstd(normal_n),2))])
    
end


%% All trials histogram
% a_2_6 = squeeze(nanmean(filtered_all_emotional_pupil_3d(:,:,2500:6000),3));
% b_2_6 = squeeze(nanmean(filtered_all_neutral_pupil_3d(:,:,2500:6000),3));
% normal_e = round(sqrt(a_2_6)/10,2);
% normal_n = round(sqrt(b_2_6)/10,2);
% normal_e = reshape(normal_e,1320,1);
% normal_n = reshape(normal_n,1320,1);

violinplot([normal_e,normal_n])
%% Standard deviation
clc
plot_n=0;
plots_cell = cell(3);
for dtp_i = [ 0.2]
    plot_n=plot_n+1;
    filtered_all_neutral_pupil_3d = nan(size(all_neutral_pupil_3d));
    filtered_all_emotional_pupil_3d = nan(size(all_emotional_pupil_3d));
    for s= 1:55
        for i = 1:24
            if dtp(s,i)>= dtp_i && dtp(s,i) <= 0.999
                if ~isnan(all_neutral_pupil_3d(s,i))
                   filtered_all_neutral_pupil_3d(s,i,1:6000) = all_neutral_pupil_3d(s,i,1:6000);
                else
                   filtered_all_emotional_pupil_3d(s,i,1:6000) = all_emotional_pupil_3d(s,i,1:6000);
                end
            end
        end
    end
    
    % Anova test std
    a_2_6 = squeeze(nanstd(filtered_all_emotional_pupil_3d(:,:,2500:6000),0,2));
    b_2_6 = squeeze(nanstd(filtered_all_neutral_pupil_3d(:,:,2500:6000),0,2));
    e_2_6 = squeeze(nanmean(a_2_6,2));
    n_2_6 = squeeze(nanmean(b_2_6,2));
    normal_e = round(sqrt(e_2_6)/10,2);
    normal_n = round(sqrt(n_2_6)/10,2);
    %
    pupil_data_2_6 = [normal_e, normal_n];

    %anova1(pupil_data_2_6)

    % ttest
    [h,p,ci,stats] = ttest(normal_e, normal_n);
    l= {'Emotional Stimuli','Neutral Stimuli'};
    % figure
    % boxplot([normal_e, normal_n],'Labels',l)
    
    % violin plot
    % figure
    %violinplot([normal_e, normal_n]);
    disp(['dtp>=' num2str(dtp_i) '--h:' num2str(round(h,3)) '--p:' num2str(round(p,3)) '--t:' num2str(round(stats.tstat,2))...
        ,'--Me:' num2str(round(nanmean(normal_e),2)) '±' num2str(round(nanstd(normal_e),2))...
        ,'--Mn:' num2str(round(nanmean(normal_n),2)) '±' num2str(round(nanstd(normal_n),2))])
    
end

%% plot without distractions --items--
plot_n=0;
figure
for dtp_i = [0.2];
    plot_n=plot_n+1;
    filtered_all_neutral_pupil_3d = nan(size(all_neutral_pupil_3d));
    filtered_all_emotional_pupil_3d = nan(size(all_emotional_pupil_3d));
    removed_subjects = [5,18,20,22,24,32,34,35,36,37,43];
%     for s = setdiff([1:sample_size],removed_subjects);
    for s= 1:55
        for i = 1:24
            if dtp(s,i)>= dtp_i && dtp(s,i) <= 0.999
                if ~isnan(all_neutral_pupil_3d(s,i))
                   filtered_all_neutral_pupil_3d(s,i,1:6000) = all_neutral_pupil_3d(s,i,1:6000);
                else
                   filtered_all_emotional_pupil_3d(s,i,1:6000) = all_emotional_pupil_3d(s,i,1:6000);
                end
            end
        end
    end

    i = 10;
    a = squeeze(nanmean(filtered_all_emotional_pupil_3d,1));
    b = squeeze(nanmean(filtered_all_neutral_pupil_3d,1));

    % figure
    % plot(a(i,:)','color',[1 0 0])
    % hold on
    % plot(b(i,:)','color',[0 0 1])
    % legend Emotional Neutral

    e = squeeze(nanmean(a(:,:),1));
    n = squeeze(nanmean(b(:,:),1));
    normal_e = round(sqrt(e)/10,3);
    normal_n = round(sqrt(n)/10,3);
    
    
    subplot(2,2,plot_n)
    plot(normal_e','color',[1 0 0]) %Emotional
    hold on
    plot(normal_n','color',[0 0 1]) %Neutral
    legend Emotional Neutral
        
    % Anova test mean
    a_2_6 = squeeze(nanmean(filtered_all_emotional_pupil_3d(:,:,2500:6000),1));
    b_2_6 = squeeze(nanmean(filtered_all_neutral_pupil_3d(:,:,2500:6000),1));
    e_2_6 = squeeze(nanmean(a_2_6,2));
    n_2_6 = squeeze(nanmean(b_2_6,2));
    normal_e = round(sqrt(e_2_6)/10,2);
    normal_n = round(sqrt(n_2_6)/10,2);
    %
    pupil_data_2_6 = [normal_e, normal_n];

    %anova1(pupil_data_2_6)

    % ttest
    [h,p,ci,stats] = ttest(normal_e, normal_n)
    l= {'Emotional Stimuli','Neutral Stimuli'};
    % figure
    % boxplot([normal_e, normal_n],'Labels',l)

    % violin plot
    % figure
    %violinplot([normal_e, normal_n]);
    title(['dp<=' num2str(dtp_i) '--p:' num2str(round(p,3)) '--t:' num2str(round(stats.tstat,2))])
    
end