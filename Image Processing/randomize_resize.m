function x = randomize_resize(img)
%input image
%output: 423x423 img 
    x = imresize(img,[423,30]);
    x = reshape(x,[900,3]);
    shuffled_index = randperm(900);
    x = x(shuffled_index,:);
    x = reshape(x,[30,30,3]);
    x = imresize(x,[423,423],'method','nearest');
    imshow
end

