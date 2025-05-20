function x = randomize_img(img)
%input image
%output: 423x423 randomized img with 30x30 pixels
    x = imresize(img,[23,23]);
    x = reshape(x,[529,3]);
    shuffled_index = randperm(529);
    x = x(shuffled_index,:);
    x = reshape(x,[23,23,3]);
    x = imresize(x,[423,423],'method','nearest');
end

