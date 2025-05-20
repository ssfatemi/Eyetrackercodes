function c = img_complexity(img,item_mask)

 
 x = rgb2gray(img);
 cedge = edge(x,'Canny');
 
 %c = sum(sum(cedge))/(size(img,1)*size(img,2));
 c = sum(sum(cedge))/sum(sum(item_mask));
end

