clc
clear all
%% load data

path = 'samdata/';
listing = dir(path);
s = size(listing,1);

sams = cell(s-2,1);
saccades = cell(s-2,1);
trials = cell(s-2,1);

for i = 3:s-1
% for i = 57:57
    subject_char = num2str(i-2,'%03.f')
    file_name = [path,subject_char,'.xlsx'];
    T = readtable(file_name);
    sams{i-2} = T;
end

%%


%% load other data
load('data_56.mat')
load('subject_selection47_24.mat');
load('subject_info.mat');
corrected_l_r = load('corrected_l_r.mat');
corrected_l_r = corrected_l_r.d_rand_imgs;

%%



%% add 6-13 trials
clc
s = size(trials,1)-1;
% sample_size = s;
sample_size = 55;


% each variable saves emotional and neutral fix_duration, fix_label
neutral_sams =cell(sample_size,24);
emotional_sams =cell(sample_size,24);

for s = 31:sample_size
    T = sams{s};
    trial = 0;
    for i = 1:36
        if not(isnan(subject_selection(6,i))) 
            trial = trial + 1;
            if ((subject_selection(s,i))==1)
                x = T(table2array(T(:,'TRIAL_INDEX'))==trial,[40, 44]);
                emotional_sams{s,trial} = x;
            else
                x = T(table2array(T(:,'TRIAL_INDEX'))==trial,[40, 44]);
                neutral_sams{s,trial} = x;
            end
        end
    end
end

msg = 'done 1'

% add first 5 trials
% each variable saves emotional and neutral fix_duration, fix_label

for s = 1:5
    T = sams{s};
    trial = 0;
    for i = 1:36
        if not(isnan(subject_selection(6,i))) 
            trial = trial + 1;
            if ((subject_selection(s,i))==1)
                x = T(table2array(T(:,'TRIAL_INDEX'))==i,[40, 44]);
                emotional_sams{s,trial} = x;
            else
                x = T(table2array(T(:,'TRIAL_INDEX'))==i,[40, 44]);
                neutral_sams{s,trial} = x;
            end
        end
    end
end

% saccade columns: amp,start,end,ang,vel
% saccade columns: 4,44,17,5,6

% clear 'fixations'
% clear 'saccades'
% clear 'trials'

msg='done 2'
%%
pupil_data2.neutral_pupil = neutral_sams;
pupil_data2.emotional_pupil = emotional_sams;
neutral_sams2 = neutral_sams;
emotional_sams2 = emotional_sams;
%% save data
clc
save('sams_data_31_55.mat','neutral_sams2','emotional_sams2','pupil_data2')
msg = 'data saved'

%% Continue loading

clc
clear all
% load essentials

load('data_56.mat')
load('subject_selection47_24.mat');
load('subject_info.mat');
corrected_l_r = load('corrected_l_r.mat');
corrected_l_r = corrected_l_r.d_rand_imgs;

% load saved sams
load('sams_data_31_55.mat')
load('sams_data_55.mat')


%% combine smas data
clc
for i = 1:30
    for j = 1:24
        emotional_sams2{i,j} = emotional_sams{i,j};
        neutral_sams2{i,j} = neutral_sams{i,j};
    end
end
msg = "done"

%% saving new sams
pupil_data.neutral_pupil = neutral_sams2;
pupil_data.emotional_pupil = emotional_sams2;
neutral_sams = neutral_sams2;
emotional_sams = emotional_sams2;
%
save('sams_data_55.mat','neutral_sams','emotional_sams','pupil_data')
msg = 'data saved'

%% Start Processing
clc
clear all
% load essentials

load('data_56.mat')
load('subject_selection47_24.mat');
load('subject_info.mat');
corrected_l_r = load('corrected_l_r.mat');
corrected_l_r = corrected_l_r.d_rand_imgs;

% load sams
load('sams_data_55.mat')

%%
sample_size = 55;

all_emotional_pupil_3d = nan(sample_size,24,8000);
all_neutral_pupil_3d = nan(sample_size,24,8000);

all_emotional_pupil_3d_aoi = nan(sample_size,24,8000);
all_neutral_pupil_3d_aoi = nan(sample_size,24,8000);

emotional_aoi_indexes = nan(sample_size,24,8000);
neutral_aoi_indexes = nan(sample_size,24,8000);
aoi_indexes = nan(sample_size,24,8000);

removed_subjects = [5,18,24,34,35,36,37,43];
removed_subjects = [5,18,20,22,24,32,34,35,36,37,43];
for s = setdiff([1:sample_size],removed_subjects);
    trial = 0;
    for i = 1:24
        if (corrected_l_r(s,i)==1)
            lrbox_correction =1;
        else
            lrbox_correction =0;
        end
        
        if not(isempty(neutral_sams{s,i}))
            %pupil size
            T = neutral_sams{s,i};
            f = table2array(T(:,1));
            %index correction
            if (lrbox_correction ==1)
                indexes = find_index(f,'rand_img');
            else
                indexes = find_index(f,'stimulus_img');
            end
            %indexes = indexes | 1;
            x = table2array(T(indexes,2));
            x = str2double(x);
            if isnan(x)
                x = table2array(T(indexes,2));
            end
            all_neutral_pupil_3d_aoi(s,i,1:size(x,1)) = x;
            neutral_aoi_indexes(s,i,1:size(indexes,2)) = indexes;
            aoi_indexes(s,i,1:size(indexes,2)) = indexes;
            
            x = table2array(T(:,2));
            x = str2double(x);
            if isnan(x)
                x = table2array(T(:,2));
            end
            all_neutral_pupil_3d(s,i,1:size(x,1)) = x;
            
        else
            %pupil size
            T = emotional_sams{s,i};
            f = table2array(T(:,1));
            %index correction
            if (lrbox_correction ==1)
                indexes = find_index(f,'rand_img');
            else
                indexes = find_index(f,'stimulus_img');
            end
            %indexes = indexes | 1;
            x = table2array(T(indexes,2));
            x = str2double(x);
            if isnan(x)
                x = table2array(T(indexes,2));
            end
            all_emotional_pupil_3d_aoi(s,i,1:size(x,1)) = x;
            emotional_aoi_indexes(s,i,1:size(indexes,2)) = indexes;
            aoi_indexes(s,i,1:size(indexes,2)) = indexes;
            
            x = table2array(T(:,2));
            x = str2double(x);
            if isnan(x)
                x = table2array(T(:,2));
            end
            all_emotional_pupil_3d(s,i,1:size(x,1)) = x;
            
            
        end
    end
    s
end

clc
msg = "Features Assigned"
%%
save('pupillary_data_55.mat')
%%
clc
clear all
load('pupillary_data_55.mat')
%%
clc
distraction = nan(55,24);
for s= 1:55
    for i = 1:24
    total_time = sum(~isnan(emotional_aoi_indexes(s,i,1:6000))) + sum(~isnan(neutral_aoi_indexes(s,i,1:6000)));
    aoi_time = sum(emotional_aoi_indexes(s,i,1:6000)==1) + sum(neutral_aoi_indexes(s,i,1:6000)==1);
    not_aoi_time = sum(emotional_aoi_indexes(s,i,1:6000)==0) + sum(neutral_aoi_indexes(s,i,1:6000)==0);
    distraction_percent = round(not_aoi_time/total_time,2);
    
    distraction(s,i)= distraction_percent;
    end
end
msg = 'distraction extracted'
% Distraction Samples plot
for i=1:100
    d_plot(i) = sum(sum(distraction<(i/100)));
end
plot(d_plot)
%% distraction level 1-6
clc
distraction = nan(55,24);
for s= 1:55
    for i = 1:24
    total_time = sum(~isnan(emotional_aoi_indexes(s,i,1:6000))) + sum(~isnan(neutral_aoi_indexes(s,i,1:6000)));
    aoi_time = sum(emotional_aoi_indexes(s,i,1:6000)==1) + sum(neutral_aoi_indexes(s,i,1:6000)==1);
    not_aoi_time = sum(emotional_aoi_indexes(s,i,1000:6000)==0) + sum(neutral_aoi_indexes(s,i,1000:6000)==0);
    distraction_percent = round(not_aoi_time/6000,2);
    
    distraction(s,i)= distraction_percent;
    end
end
msg = 'distraction extracted'
% Distraction Samples plot
for i=1:100
    d_plot(i) = (sum(sum(distraction<(i/100)))-264)/1056;
end
plot(d_plot)

%% plot without distractions
clc
males = logical(subject_info(1:55,3));
plot_n=0;

plots_cell = cell(3,1);
figure
for dp = [ 0.3 0.5 0.8];
    plot_n=plot_n+1;
    filtered_all_neutral_pupil_3d = nan(size(all_neutral_pupil_3d));
    filtered_all_emotional_pupil_3d = nan(size(all_emotional_pupil_3d));
    removed_subjects = [5,18,20,22,24,32,34,35,36,37,43];
%     for s = setdiff([1:sample_size],removed_subjects);
    for s= 1:55
        for i = 1:24
            if distraction(s,i)<dp
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
    
    
    subplot(3,3,plot_n)
    plot(normal_e','color',[1 0 0]) %Emotional
    hold on
    plot(normal_n','color',[0 0 1]) %Neutral
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
    title(['dp<=' num2str(dp) '--p:' num2str(round(p,3)) '--t:' num2str(round(stats.tstat,2))])
    disp(['dp<=' num2str(dp) '--h:' num2str(round(h,3)) '--p:' num2str(round(p,3)) '--t:' num2str(round(stats.tstat,2))...
        ,'--Me:' num2str(round(nanmedian(normal_e),2)) '±' num2str(round(nanstd(normal_e),2))...
        ,'--Mn:' num2str(round(nanmedian(normal_n),2)) '±' num2str(round(nanstd(normal_n),2))])
    
end
%% Standard deviation
clc
plot_n=0;
plots_cell = cell(3);
for dp = [0.2 0.5 0.8];
    plot_n=plot_n+1;
    filtered_all_neutral_pupil_3d = nan(size(all_neutral_pupil_3d));
    filtered_all_emotional_pupil_3d = nan(size(all_emotional_pupil_3d));
    for s= 1:55
        for i = 1:24
            if distraction(s,i)<dp
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
    disp(['dp<=' num2str(dp) '--h:' num2str(round(h,3)) '--p:' num2str(round(p,3)) '--t:' num2str(round(stats.tstat,2))...
        ,'--Me:' num2str(round(nanmean(normal_e),2)) '±' num2str(round(nanstd(normal_e),2))...
        ,'--Mn:' num2str(round(nanmean(normal_n),2)) '±' num2str(round(nanstd(normal_n),2))])
    
end

%% plot without distractions --items--
plot_n=0;
figure
for dp = [0.9];
    plot_n=plot_n+1
    filtered_all_neutral_pupil_3d = nan(size(all_neutral_pupil_3d));
    filtered_all_emotional_pupil_3d = nan(size(all_emotional_pupil_3d));
    removed_subjects = [5,18,20,22,24,32,34,35,36,37,43];
%     for s = setdiff([1:sample_size],removed_subjects);
    for s= 1:55
        for i = 1:24
            if distraction(s,i)<dp
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
    a_2_6 = squeeze(nanmean(filtered_all_emotional_pupil_3d(males,:,2500:6000),1));
    b_2_6 = squeeze(nanmean(filtered_all_neutral_pupil_3d(males,:,2500:6000),1));
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
    title(['dp<=' num2str(dp) '--p:' num2str(round(p,3)) '--t:' num2str(round(stats.tstat,2))])
    
end