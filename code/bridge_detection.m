clc; close all; clear;tic
%% read data
fold_location = 'E:\作业2-数据';
data_name = '\CSA-bridge.jpg';
% fold_location = '/Users/weiyihai/中国科学院/中国科学院大学/国科大博一下学期资料数据/SAR信号处理与运动补偿22-23春季/作业/作业2-数据';
% data_name = '/CSA-bridge.jpg';
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

%%
% 显示原始图像和应用阈值后的图像
T = graythresh(equalizedImg);
thresholdImg = imbinarize(equalizedImg, T);

fgImg=[];
for i = 1:size(equalizedImg, 1)
    for j = 1:size(equalizedImg, 2)
        if thresholdImg(i, j) ==1
        else
            fgImg(end+1) = equalizedImg(i, j);
        end
    end
end

% sigm1=otus(fgImg)/255;
T2=graythresh(uint8(fgImg));

finalImg = imbinarize(equalizedImg , T2);
% 先腐蚀
se = strel('disk', 1); % 创建一个5像素半径的圆形结构元素
% se = strel('rectangle', [6 3]);
%  se = strel('line', 3, 0);
eroded_img = imerode(finalImg, se);

% 进行空洞填充操作
bw_filled = imfill(eroded_img,"holes");

% 进行孤立点去除操作
bw_without_isolated = bwareaopen(bw_filled, 1);

% 显示结果
% figure;
% subplot(2, 1, 1);
% imshow(bw_filled);
% title('填充空洞');
% subplot(2, 1, 2);
% imshow(bw_without_isolated);
% title('去除孤立点后的图像');
neg_dilated_img = imcomplement(bw_without_isolated);%图片取反
% 创建一个半径为 3 的圆形结构元素
se = strel('disk', 3);

% 对图像进行开运算处理，去除毛糙点
bw_opened = imopen(neg_dilated_img, se);
% figure;imshow(bw_opened);
% title('开运算处理，去除毛糙点');
% 
% %膨胀
% se = strel('disk', 2); % 创建一个10像素半径的圆形结构元素
% % se = strel('rectangle', [2 2]);
% % se = strel('line', 3, 0);
% dilated_img = imdilate(neg_dilated_img, se);
% % figure;
% % imshow(dilated_img);title('膨胀')
% %腐蚀
% % se = strel('disk', 2); % 创建一个5像素半径的圆形结构元素
% % se = strel('rectangle', [6 3]);
%  se = strel('line', 3, 0);
% eroded_img = imerode(dilated_img, se);
% figure;
% imshow(eroded_img);title('腐蚀')


% 创建一个半径为 5 的圆形结构元素
se = strel('disk', 2);
% 对图像进行闭运算处理
bw_closed = imclose(bw_opened, se);
% figure;
% imshow(bw_closed);title('闭运算')

% 进行连通域标记，得到各个独立的连接分量
[label, num] = bwlabel(bw_closed,4);


% 计算每个连接分量的面积
stats = regionprops(label, 'Area');

% 将所有面积值放入一个数组
areas = [stats.Area];

% 使用 Otsu 分析得到面积阈值
level = (graythresh(uint8(areas)))*max(areas);
delta  = 0.1*level;

% 根据阈值剔除面积过小的连接分量
bw_clean = ismember(label, find((areas>=level-3*delta)));

% 
% % 显示结果
% figure;
% imshow(bw_clean);
% figure;imshow(label);
% 进行逻辑与运算
bw_and = bw_clean & bw_closed;
figure;imshow(bw_and);title('与逻辑运算');
% 对二值图像进行形态学细化操作
bw_thinned = bwmorph(bw_clean, 'thin', Inf);

figure;imshow(bw_thinned);title('细化操作,提取主干线');
% 将主干线连接在一起
bw_filled = imfill(bw_thinned, 'holes');
bw_filled = bwmorph(bw_filled, 'bridge');

% 对河流边界线和河流主干线进行逻辑与运算，得到交点
edges = edge(bw_clean, 'canny');
figure;imshow(edges);title('提取河流边界线');

% 将两张图像重叠显示在一张图上
figure;
imshow(bw_thinned);
hold on;
h = imshow(edges);
set(h, 'AlphaData',0.4); % 设置边界线的透明度
title('河流主干线和边界线重叠显示');
intersection = bitand(edges, bw_thinned);

% 对交点进行配对处理，确定每个桥梁的中心位置
[rows, cols] = find(intersection);
d = pdist([rows, cols]);
idx = find(d < 20); % 这里的 20 是一个阈值，用于确定交点间的距离
pairs = nchoosek(idx, 2);
centers = (rows(pairs(:, 1)) + rows(pairs(:, 2))) / 2;








