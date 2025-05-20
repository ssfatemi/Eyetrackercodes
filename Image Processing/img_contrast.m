function cont = img_contrast(img)

% %  tg = NaN(size(r_channel));
% %  tg(item_mask) = x(item_mask);
 
 x = rgb2gray(img);
 img_gr_double = double(x);
 cont = std2(img_gr_double);
end

