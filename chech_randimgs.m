%% 
clc
% what used in experiment:
used_rand_imgs = rand_imgs(:,1:24);
%the NaN values considered 0 (left img is rand) for instance trials 4,5
%always had left rnd img
used_rand_imgs(isnan(used_rand_imgs))= 0;

nn_rand_imgs = rand_imgs(:,~isnan(rand_imgs(6,:)));
d_rand_imgs = abs(used_rand_imgs - nn_rand_imgs);
sum(nansum(d_rand_imgs))