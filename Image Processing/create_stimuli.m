function img = create_stimuli(subject_index,neutral_path,emotional_path)

saving_path = ['subjects_stimuli\' , num2str(subject_index), '\'];
load('items.mat');
load('subject_selection50_edited.mat');


w3 = 455;
h2 = 384;
stimului_number =0;
    for i = 1:size(items,1)
        stimului_number = stimului_number + 1;

        item = items{i};
        item_file = [item, '.jpg'];

        if subject_selection(subject_index,i)==1
            path = emotional_path;
        else
            path = neutral_path;
        end

        main_img = imread([path, item_file]);
        rand_img = randomize_img(main_img);
        
        img = 128 * ones(768,1366,3);

        img_w =size(main_img,2);
        img_h =size(main_img,1);

        shift_right = -50;
        shift_left = 50;
        
        if rand_imgs(subject_index,i)==1
           r_img = rand_img;
           l_img = main_img;
        else
           r_img = main_img;
           l_img = rand_img;
        end
        
        img((h2-img_h/2):(h2+img_h/2-1),(w3-img_w/2+shift_right):(w3+img_w/2-1+shift_right),:) = double(l_img);
        img((h2-img_h/2):(h2+img_h/2-1),(2*w3-img_w/2+shift_left):(2*w3+img_w/2-1+shift_left),:) = double(r_img);
        
        img_saving_path =[saving_path, num2str(stimului_number),'.bmp'];
        imwrite(uint8(img), img_saving_path);
    end

end

