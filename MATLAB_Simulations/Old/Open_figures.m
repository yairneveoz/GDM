% Open_figures
clear
clc

% Read the JPG image
path1 = 'C:\Users\yairn\Downloads\2023_12_05\';
image_name = 'IMG_0001.tif';
% 
img = imread([path1,image_name]);

% Display the image
figure(1)
clf
imshow(img);

if 0
    for image_ind = 1:9
        image_name = ['IMG_000',num2str(image_ind),'.jpg'];
        img = imread([path1,image_name]);

        % Display the image
        figure(1)
        imshow(img);
        title(['Image ',num2str(image_ind),'.jpg'])
        pause(0.3)
        drawnow
    end
end

