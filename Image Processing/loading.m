%% load neutral images
clc
clear all
%path = 'neutral_423/';
path = 'selection1/';

listing = dir(path);
s = size(listing,1);
all_gray_imgs = [];
all_imgs = [];
back = 242;
adjusted_imgs =[];

stimuli_cell = cell(35,4);
stimuli_info = cell(35,4);
stimuli_cell_labels = {'img','gray_img','item_mask','item_brightness_adjusted'};
stimuli_info_labels = {'brightness','adjusted_brightness','item_mask','item_brightness_adjusted'};

items=cell(1,1);

for i = 3:s
%for i = 3:3
    item = listing(i).name(1:end-4);
    img = imread([path, listing(i).name]);
    stimuli_cell{i-2,1} = img;

    img_gray = rgb2gray(img);
    stimuli_cell{i-2,2} = img_gray;
    
    w_bakcground = img_gray > back;
    stimuli_cell{i-2,3} = not(w_bakcground);
    
    a = -10;
    adjusted_img = img;
    adjusted_img(not(w_bakcground)) = img(not(w_bakcground)) + a;
    stimuli_cell{i-2,4} = adjusted_img;
    
    %all_imgs = cat(4,all_imgs,img);
    %all_gray_imgs = cat(3,all_gray_imgs,img_gray);
    %adjusted_imgs = cat(4,adjusted_imgs,adjusted_img);
    
    items{i-2,1}= item
end

% for j= 1:35
%     figure(j); imshow(stimuli_cell{j,3});
% end

%% Brightness
all_brightness = [];
for i = 1:size(stimuli_cell,1)
    nanback_all_gray_img = stimuli_cell{i,2};
    nanback_all_gray_img(not(stimuli_cell{i,3}))= NaN;
    luminance = squeeze(nanmean(nonzeros(nanback_all_gray_img),[1 2]));
    brightness = back - luminance;
    stimuli_info{i,1} = brightness;
    
    nanback_all_gray_img = rgb2gray(stimuli_cell{i,4});
    nanback_all_gray_img(not(stimuli_cell{i,3}))= NaN;
    luminance2 = squeeze(nanmean(nonzeros(nanback_all_gray_img),[1 2]));
    brightness2 = back - luminance2;
    stimuli_info{i,2} = brightness2;
    
    all_brightness=[all_brightness; brightness];
    
end

[mean(all_brightness) std(all_brightness)]

%% RGB Mean
clc
all_rgb_mean = [];
item_mean = [];
for i = 1:size(stimuli_cell,1)
    i
    item_mask = stimuli_cell{i,3};
    img = stimuli_cell{i,1};
    [m s] = color_proportion(img,item_mask);
    all_rgb_mean = [all_rgb_mean; m];
end
mean(all_rgb_mean,1)
std(all_rgb_mean,1)

%% img complexity
clc
all_complexity = [];
for i = 1:size(stimuli_cell,1)
    i
    item_mask = stimuli_cell{i,3};
    img = stimuli_cell{i,1};
    c = img_complexity(img,item_mask);
    all_complexity = [all_complexity; c];
end
mean(all_complexity,1)
std(all_complexity,1)

%% Contrast
clc
all_contrast = [];
for i = 1:size(stimuli_cell,1)
    i
    %item_mask = stimuli_cell{i,3};
    img = stimuli_cell{i,1};
    cont = img_contrast(img);
    all_contrast = [all_contrast; cont];
end
mean(all_contrast,1)
std(all_contrast,1)

%% Saving 
%save('neutral_info.mat','all_brightness','all_rgb_mean','all_contrast','all_complexity');
save('emotional_info.mat','all_brightness','all_rgb_mean','all_contrast','all_complexity');

