clc; close all; clear;tic
%% read data
% fold_location = 'E:\作业2-数据';
% data_name = '\CSA-bridge.jpg';
fold_location = '/Users/weiyihai/中国科学院/中国科学院大学/国科大博一下学期资料数据/SAR信号处理与运动补偿22-23春季/作业/作业2-数据';
data_name = '/CSA-bridge.jpg';
filename = strcat(fold_location,data_name);
% 读取图像
img = imread(filename);
normalizedImg=rgb2gray(img);
%% 中值滤波去噪
% 将图像转换为灰度图像
% grayImg = rgb2gray(img);
grayImg = normalizedImg;
% 定义中值滤波器的大小
filterSize = 3;

% 应用中值滤波器
medianImg = medfilt2(grayImg, [filterSize filterSize]);
% % 显示原始图像和处理后的图像
% subplot(1,2,1), imshow(grayImg), title('原始图像');
% subplot(1,2,2), imshow(medianImg), title('中值滤波后的图像');
%%
% 对灰度图像进行直方图均衡化
equalizedImg = histeq(medianImg);
% % 计算灰度直方图
% histogram = imhist(equalizedImg);
% 
% % 可视化直方图
% bar(histogram);
% title('灰度直方图');
% xlabel('灰度级别');
% ylabel('像素数');
% 显示原始图像和均衡化后的图像
% figure;
% subplot(1, 2, 1);
% imshow(grayImg);
% title('原始图像');
% subplot(1, 2, 2);
% imshow(equalizedImg);
% title('直方图均衡化后的图像');
%%
% 显示原始图像和应用阈值后的图像
T = graythresh(equalizedImg);
sigm=otus(equalizedImg)/255;
thresholdImg = imbinarize(equalizedImg, sigm);
% figure;
% subplot(1, 2, 1);
% imshow(equalizedImg);
% title('原始图像');
% subplot(1, 2, 2);
% imshow(thresholdImg);
% title(['阈值为 ' num2str(T)]);

fgImg=[];
for i = 1:size(equalizedImg, 1)
    for j = 1:size(equalizedImg, 2)
        if thresholdImg(i, j) ==1
        else
            fgImg(end+1) = equalizedImg(i, j);
        end
    end
end

sigm1=otus(fgImg)/255;

finalImg = imbinarize(equalizedImg, sigm1);
% figure;
% subplot(1, 2, 1);
% imshow(thresholdImg);
% title(['T阈值为 ' num2str(T)]);
% subplot(1, 2, 2);
% imshow(finalImg);
% title(['T’阈值为 ' num2str(sigm1)]);

%膨胀
% se = strel('disk', 1); % 创建一个5像素半径的圆形结构元素
se = strel('rectangle', [1 1]);
eroded_img = imerode(finalImg, se);
%腐蚀
% se = strel('disk', 5); % 创建一个10像素半径的圆形结构元素
% se = strel('rectangle', [5 5]);
se = strel('line', 3, 0);
dilated_img = imdilate(eroded_img, se);
% 将图像求反并显示
neg_dilated_img = imcomplement(dilated_img);
figure;
imshow(neg_dilated_img);


% 进行连通域标记
[L, num] = bwlabel(neg_dilated_img, 8);

% 统计每个连通域的面积
stats = regionprops('table', L, 'Area');
areas = stats.Area;

