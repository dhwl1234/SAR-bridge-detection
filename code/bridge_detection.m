clc; close all; clear;tic
%% read data
fold_location = 'E:\作业2-数据';
data_name = '\CSA-bridge.jpg';
filename = strcat(fold_location,data_name);
% 读取图像
img = imread(filename);
normalizedImg=img;
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
thresholdImg = imbinarize(equalizedImg, T);
% figure;
% subplot(1, 2, 1);
% imshow(equalizedImg);
% title('原始图像');
% subplot(1, 2, 2);
% imshow(thresholdImg);
% title(['阈值为 ' num2str(T)]);
% bgImg=zeros(size(equalizedImg, 1),size(equalizedImg, 2));
% fgImg=zeros(size(equalizedImg, 1),size(equalizedImg, 2));
% 
% % 遍历图像，将像素值大于阈值的像素保存到数组 bgImg 中，将像素值小于等于阈值的像素保存到数组 fgImg 中
% for i = 1:size(equalizedImg, 1)
%     for j = 1:size(equalizedImg, 2)
%         if equalizedImg(i, j) > T*255
%             bgImg(i,j) = equalizedImg(i, j);
%         else
%             fgImg(i,j) = equalizedImg(i, j);
%         end
%     end
% end
% % 对于灰度级在 [0, T] 的像素，使用 Otsu 准则求出分割阈值 T'
% T1 = graythresh(fgImg);

fgImg=[];
for i = 1:size(equalizedImg, 1)
    for j = 1:size(equalizedImg, 2)
        if equalizedImg(i, j) > T*255
        else
            fgImg(end+1) = equalizedImg(i, j);
        end
    end
end

% 计算灰度图像的直方图
histCounts = imhist(fgImg/255);
% 计算灰度图像像素的总数
numPixels = numel(grayImg);
% 构造灰度级别向量
pixelVals = [0:255]';
% 计算直方图的累积分布函数
omega = cumsum(histCounts) / numPixels;
% 计算灰度级别的平均值
mu = cumsum(pixelVals .* histCounts) / numPixels;
% 计算灰度级别的总平均值
muT = mu(end);
% 计算前景和背景的方差
sigmaB = (muT * omega - mu).^2 ./ (omega .* (1 - omega));
% 选择最大方差对应的灰度级别作为阈值
[maxval, idx] = max(sigmaB);
T1 = pixelVals(idx) / 255;


finalImg = imbinarize(equalizedImg, T1);
figure;
subplot(1, 2, 1);
imshow(thresholdImg);
title(['T阈值为 ' num2str(T)]);
subplot(1, 2, 2);
imshow(finalImg);
title(['T’阈值为 ' num2str(T1)]);

