%% Pupil Statistics
clc

emotional_pupil_data = all_emotional_pupil_3d(:,:,1:6000);
neutral_pupil_data = all_neutral_pupil_3d(:,:,1:6000);
emotional_pupil_data_aoi = all_emotional_pupil_3d_aoi(:,:,1:6000);
neutral_pupil_data_aoi = all_neutral_pupil_3d_aoi(:,:,1:6000);

neutral_pupil_mean_aoi_0_6 = nanmean(neutral_pupil_data_aoi,3);
emotional_pupil_mean_aoi_0_6 = nanmean(emotional_pupil_data_aoi,3);
neutral_pupil_mean_aoi_1_6 = nanmean(neutral_pupil_data_aoi(:,:,1000:6000),3);
emotional_pupil_mean_aoi_1_6 = nanmean(emotional_pupil_data_aoi(:,:,1000:6000),3);

items_pupil_mean_aoi_1_6 = [nanmean(emotional_pupil_mean_aoi_1_6,1); nanmean(neutral_pupil_mean_aoi_1_6,1)];
subs_pupil_mean_aoi_1_6 = [nanmean(emotional_pupil_mean_aoi_1_6,2), nanmean(neutral_pupil_mean_aoi_1_6,2)];
items_pupil_mean_aoi_1_6 = items_pupil_mean_aoi_1_6';

items_pupil_mean_aoi_0_6 = [nanmean(emotional_pupil_mean_aoi_0_6,1); nanmean(neutral_pupil_mean_aoi_0_6,1)];
subs_pupil_mean_aoi_0_6 = [nanmean(emotional_pupil_mean_aoi_0_6,2), nanmean(neutral_pupil_mean_aoi_0_6,2)];
items_pupil_mean_aoi_0_6 = items_pupil_mean_aoi_0_6';

%% ttests
[h3,p3,~,~] = ttest(subs_pupil_mean_aoi_0_6(:,1),subs_pupil_mean_aoi_0_6(:,2))
l= {'Emotional Stimuli','Neutral Stimuli'};
boxplot([subs_pupil_mean_aoi_0_6(:,1),subs_pupil_mean_aoi_0_6(:,2)],'Labels',l)

%% items ttest
[h3,p3,~,~] = ttest(items_pupil_mean_aoi_1_6(:,1),items_pupil_mean_aoi_1_6(:,2))
l= {'Emotional Stimuli','Neutral Stimuli'};
boxplot([items_pupil_mean_aoi_1_6(:,1),items_pupil_mean_aoi_1_6(:,2)],'Labels',l)
%%
males2 = logical(subject_info(setdiff([3:55],removed_subjects),3));
%%

[h3,p3,~,stat] = ttest(neutral_subjects_v(males2,:),neutral_subjects_v(~males2,:))

%%
scatter(emotional_items_v,emotional_items_a,[],[1 0 0],'filled')
hold on;
scatter(neutral_items_v,neutral_items_a,[],[0 0 1],'filled')
