function [m,s] = color_proportion(img,item_mask)

 x = rgb2gray(img);
 
 r_channel = double(squeeze(img(:,:,1)));
 g_channel = double(squeeze(img(:,:,2)));
 b_channel = double(squeeze(img(:,:,3)));
 
 tg = NaN(size(r_channel));
 tg(item_mask) = x(item_mask);
 
 
 r=NaN(size(r_channel));
 g=NaN(size(r_channel));
 b=NaN(size(r_channel));
 
 r(item_mask) = 0.2989*r_channel(item_mask);
 g(item_mask) = 0.587*g_channel(item_mask);
 b(item_mask) = 0.1140*b_channel(item_mask);
 
 r_ratio=NaN(size(r_channel));
 g_ratio=NaN(size(r_channel));
 b_ratio=NaN(size(r_channel));
 
 tg(tg==0)=1;
  
 r_ratio(item_mask) = r(item_mask)./tg(item_mask);
 g_ratio(item_mask) = g(item_mask)./tg(item_mask);
 b_ratio(item_mask) = b(item_mask)./tg(item_mask);
 
%  if sum(isinf(r_ratio)) > 0
%     r_ratio = 1 - g_ratio - b_ratio;
%  elseif sum(isinf(g_ratio)) > 0
%     g_ratio = 1 - r_ratio - b_ratio;
%  elseif sum(isinf(b_ratio)) > 0
%     b_ratio = 1 - r_ratio - g_ratio;
%  else
%     
%  end
 
 rm = nanmean(r_ratio,[1 2]);
 gm = nanmean(g_ratio,[1 2]);
 bm = nanmean(b_ratio,[1 2]);
 
 rs = nanstd(reshape(r_ratio,[423*423,1]));
 gs = nanstd(reshape(g_ratio,[423*423,1]));
 bs = nanstd(reshape(b_ratio,[423*423,1]));
 
 
 
 %display([rm gm bm sum([rm gm bm])]);
 %display([rs gs bs])
 
 m = [rm gm bm];
 s = [rs gs bs];
 
end

