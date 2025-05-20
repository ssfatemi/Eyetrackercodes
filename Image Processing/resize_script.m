clc
clear all
path = 'imgs/';
%path = 'selection1/';

listing = dir(path);
s = size(listing,1);

for i = 3:s
    img = imread([path, listing(i).name]);
    img = imresize(img,[423 423])
    imwrite(img,listing(i).name)
end