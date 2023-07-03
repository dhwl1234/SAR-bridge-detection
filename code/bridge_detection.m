clc; close all; clear;tic
%% read data
% 获取主文件夹路径
mainFolderPath = fileparts(fileparts(mfilename('fullpath'))); % 当前脚本所在文件夹的上一级文件夹的上一级文件夹路径
% 构建image文件夹路径
imageFolderPath = fullfile(mainFolderPath, 'image');
% 构建图像文件路径
imageFilePath = fullfile(imageFolderPath, 'CSA-bridge.jpg');
% 读取图像
img = imread(imageFilePath);
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
% figure;imshow(medianImg), title('中值滤波后的图像',FontSize=30);
%%
% 对灰度图像进行直方图均衡化
equalizedImg = histeq(medianImg);
% figure;imshow(equalizedImg), title('直方图均衡化后的图像',FontSize=30);
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
% figure;imshow(finalImg), title(['两次otsu方法确定阈值，阈值=' num2str(T2),'范围0-1'],FontSize=30);
% 先腐蚀
se = strel('disk', 1); % 创建一个1像素半径的圆形结构元素
eroded_img = imerode(finalImg, se);
% 进行空洞填充操作
bw_filled = imfill(eroded_img,"holes");
% 进行孤立点去除操作
bw_without_isolated = bwareaopen(bw_filled, 1);
% figure;imshow(bw_without_isolated), title('腐蚀、空洞填充，孤立点去除后的图片',FontSize=30);
neg_dilated_img = imcomplement(bw_without_isolated);%图片取反
% figure;imshow(neg_dilated_img), title('将图片取反后的效果',FontSize=30);
% 创建一个半径为 3 的圆形结构元素
se = strel('disk', 3);
% 对图像进行开运算处理，去除毛糙点
bw_opened = imopen(neg_dilated_img, se);
% 创建一个半径为 5 的圆形结构元素
se = strel('disk', 2);
%  se = strel('rectangle', [15 15]);
% 对图像进行闭运算处理
bw_closed = imclose(bw_opened, se);

% figure;imshow(bw_closed);title('第一次开闭运算后的图片',FontSize=30)
%% 第二次开闭运算
% 定义结构元素

se = strel('rectangle',[10 14]); % 结构元素的大小可以根据具体情况调整
% 进行形态学操作
bw_smooth = imclose(bw_closed, se); % 闭运算
bw_smooth = imopen(bw_smooth, se); % 开运算

% % 显示图像
% figure;imshow(bw_smooth);title('第二次开闭运算后的图片',FontSize=30)

% 进行连通域标记，得到各个独立的连接分量
[label, num] = bwlabel(bw_smooth,4);


% 计算每个连接分量的面积
stats = regionprops(label, 'Area');

% 将所有面积值放入一个数组
areas = [stats.Area];

% 使用 otus分析得到面积阈值
level = (graythresh(uint8(areas)))*max(areas);
delta  = 0.1*level;
% 根据阈值剔除面积过小的连接分量
bw_clean2 = ismember(label, find((areas>=level-4*delta)));


% 
% % % 显示结果
% figure;imshow(bw_clean2);title('统计面积并且去除',FontSize=30);
% 进行逻辑与运算
bw_and = bw_clean2 & bw_closed;
% figure;imshow(bw_and);title('与逻辑运算',FontSize=30);
% figure;imshow(bw_and);title('与逻辑运算结果',FontSize=30);

% 提取河流主干线
bw_thinned = bwmorph(bw_clean2, 'thin', Inf);
% figure;imshow(bw_thinned);title('提取河流主干',FontSize=30);
%提取边界线
edges = edge(bw_clean2, 'canny');
% figure;imshow(edges);title('提取河流边界线',FontSize=30);
%连接河流   
bw_thinned_connect= ConnectRiver(bw_thinned);
% % % 将河流主干线和边界线两张图像重叠显示在一张图上
% figure;imshow(bw_thinned);hold on;
% h = imshow(edges);set(h, 'AlphaData', 0.5); % 设置边界线的透明度
% title('河流主干线和边界线一起显示',FontSize=30);
%  将重叠的河流主干线和边界线两张图像重叠显示在一张图上
% figure;imshow(bw_thinned_connect);hold on;
% h = imshow(edges);set(h, 'AlphaData', 0.5); % 设置边界线的透明度
% title('重叠的河流主干线和边界线一起显示',FontSize=30);

% 找出相交的点
intersect_points = bw_thinned_connect & edges;

% 获取相交点的横纵坐标
[rows, cols] = find(intersect_points);
coordinates = [cols, rows];
A_x=coordinates(1,1);
A_y=coordinates(1,2);
B_x=coordinates(end,1);
B_y=coordinates(end,2);
%计算距离
distance = sqrt((A_x-B_x)^2+(A_y-B_y));
% 计算中点坐标
x_M = (A_x + B_x) / 2;
y_M = (A_y + B_y) / 2;
% 计算斜率
k1 = -(B_y - A_y) / (B_x - A_x);%负号是为了弥补图片索引
k2 = -1 / k1;
% 构造方程

% 计算矩形的宽度和长度
width = 1.3 * distance;
length = 1.5 * distance;
% 计算矩形的顶点坐标
theta = (atan(k1));


xTopLeft = round(x_M -1/2*length*sin(theta)-1/2*width*cos(theta));
yTopLeft =round(y_M -1/2*length*cos(theta)+1/2*width*sin(theta)) ;
xTopRight = round(x_M -1/2*length*sin(theta)+1/2*width*cos(theta));
yTopRight = round(y_M -1/2*length*cos(theta)-1/2*width*sin(theta));

xBottomLeft = round(x_M +1/2*length*sin(theta)-1/2*width*cos(theta));
yBottomLeft = round(y_M +1/2*length*cos(theta)+1/2*width*sin(theta));
xBottomRight = round(x_M +1/2*length*sin(theta)+1/2*width*cos(theta));
yBottomRight = round(y_M +1/2*length*cos(theta)-1/2*width*sin(theta));

% % 在图像上绘制矩形
% angle = atan(k1)*180/pi;
% figure;imshow(edges);  % 假设已经加载了图像
% hold on;
% plot(x_M, y_M, 'ro', 'MarkerSize', 3);
% plot(xTopLeft, yTopLeft, 'ro', 'MarkerSize', 3);
% plot(xTopRight, yTopRight, 'ro', 'MarkerSize', 3);
% plot(xBottomLeft, yBottomLeft, 'ro', 'MarkerSize', 3);
% plot(xBottomRight, yBottomRight, 'ro', 'MarkerSize', 3);
% % rectangle('Position', [xBottomLeft, yBottomLeft, width, length], 'EdgeColor', 'red', 'LineWidth', 2);
% line([xTopLeft, xTopRight], [yTopLeft, yTopRight], 'Color', 'r', 'LineWidth', 1);
% line([xTopRight, xBottomRight], [yTopRight, yBottomRight], 'Color', 'r', 'LineWidth', 1);
% line([xBottomRight, xBottomLeft], [yBottomRight, yBottomLeft], 'Color', 'r', 'LineWidth',1);
% line([xBottomLeft, xTopLeft], [yBottomLeft, yTopLeft], 'Color', 'r', 'LineWidth', 1);
% title('搜索区域',FontSize=30);
% hold off;

% 创建掩膜
mask = poly2mask([xTopLeft, xTopRight, xBottomRight, xBottomLeft], [yTopLeft, yTopRight, yBottomRight, yBottomLeft], size(edges, 1), size(edges, 2));

% 将矩形外部的像素置为零
maskedImage = (edges) .* (mask);
figure;imshow(maskedImage);title('桥梁边界',FontSize=30);

%%
% 进行连通域标记，得到各个独立的连接分量
[label, num] = bwlabel(maskedImage,8);
figure;imshow((label));
% 获取连接分量的属性
properties = regionprops(label, 'PixelList');

% 存储索引为1和2的连接分量的横纵坐标
A_class = properties(1).PixelList;
B_class = properties(2).PixelList;

% 按照Y坐标升序对数据进行排序
sortedData = sortrows(A_class, 2);

% 初始化结果数组
A_final = [0,0];

% 遍历排序后的数据
for i = 1:size(sortedData, 1)
    currentY = sortedData(i, 2);
    
    % 检查当前Y值是否已经存在于结果数组中
    if ~ismember(currentY, A_final(:, 2))
        % 当前Y值不重复，将其添加到结果数组中
        A_final = [A_final; sortedData(i, :)];
    else
        % 当前Y值已经存在于结果数组中，更新对应的X值为较大值
        index = find(A_final(:, 2) == currentY);
        if sortedData(i, 1) > A_final(index, 1)
            A_final(index, 1) = sortedData(i, 1);
        end
    end
end
% 假设你的矩阵名为matrix
A_final(1, :) = [];  % 删除第一行

%% B 按照Y坐标升序对数据进行排序
sortedData = sortrows(B_class, 2);

% 初始化结果数组
B_final = [0,0];

% 遍历排序后的数据
for i = 1:size(sortedData, 1)
    currentY = sortedData(i, 2);
    
    % 检查当前Y值是否已经存在于结果数组中
    if ~ismember(currentY, B_final(:, 2))
        % 当前Y值不重复，将其添加到结果数组中
        B_final = [B_final; sortedData(i, :)];
    else
        % 当前Y值已经存在于结果数组中，更新对应的X值为较小值
        index = find(B_final(:, 2) == currentY);
        if sortedData(i, 1) < B_final(index, 1)
            B_final(index, 1) = sortedData(i, 1);
        end
    end
end
% 假设你的矩阵名为matrix
B_final(1, :) = [];  % 删除第一行
% 显示结果
figure;imshow((maskedImage));hold on;
plot(A_final(:, 1), A_final(:, 2), 'r.'); % 用红色点绘制结果
plot(B_final(:, 1), B_final(:, 2), 'r.'); % 用红色点绘制结果
hold off; % 关闭绘图保持模式
%% 求最小距离
% 假设 A_final 和 B_final 是二维矩阵，第一列是 X 坐标，第二列是 Y 坐标

% 初始化最小距离为一个较大的值
dmin = Inf;

% 初始化最小距离时的点坐标
pointA_min = [];
pointB_min = [];

% 遍历 A_final 中的点
for i = 1:size(A_final, 1)
    pointA = A_final(i, :);
    
    % 遍历 B_final 中的点
    for j = 1:size(B_final, 1)
        pointB = B_final(j, :);
        
        % 计算点 A 和点 B 之间的距离
        dist = sqrt((pointA(1) - pointB(1))^2 + (pointA(2) - pointB(2))^2);
        
        % 更新最小距离和对应的点坐标
        if dist < dmin
            dmin = dist;
            pointA_min = pointA;
            pointB_min = pointB;
        end
    end
end;
% 假设 A_final 和 B_final 是两个集合，分别包含 x 和 y 的坐标
%% 第一种直接拟合
% 对 A_final 进行线性拟合
coeff_A = polyfit(A_final(:, 1), A_final(:, 2), 1);
slope_A = coeff_A(1);     % 斜率
intercept_A = coeff_A(2); % 截距

% 对 B_final 进行线性拟合
coeff_B = polyfit(B_final(:, 1), B_final(:, 2), 1);
slope_B = coeff_B(1);     % 斜率
intercept_B = coeff_B(2); % 截距

% 绘制原始图片
figure;
imshow(normalizedImg);
hold on;

% 绘制 A_final 的拟合直线用y算x
yA = min(A_final(:, 2)):max(A_final(:, 2));
xA = (yA-intercept_A)/slope_A;
plot(xA, yA, 'r', 'LineWidth', 2);  % 以红色线条绘制直线

% 绘制 B_final 的拟合直线y算x
yB = min(B_final(:, 2)):max(B_final(:, 2));
xB = (yB-intercept_B)/slope_B;
plot(xB, yB, 'g', 'LineWidth', 2);  % 以绿色线条绘制直线

title('最小二乘法拟合桥梁边界',FontSize=30)
hold off;
%% 第二种方法选取集合中与平行线相差1-3个像素的集合
% %求A点
% C_a =A_y;%过A点的斜距
% A_Line_y = min(A_final(:,2)):max(A_final(:,2));
% A_Line_x = (A_Line_y-C_a)/k2+A_x;
% A_fit=[];
% for i = 1:size(A_Line_x,2)
%     if(abs(A_Line_x-A_final(i,1))<=3)
%         A_fit=[A_fit;[A_final(i,1),A_final(i,2)]];
%     end
% end
% %求B点
% C_b =B_y;%过B点的斜距
% B_Line_y = min(B_final(:,2)):max(B_final(:,2));
% B_Line_x = (B_Line_y-C_b)/k2+B_x;
% B_fit=[];
% for i = 1:size(B_Line_x,2)
%     if(abs(B_Line_x-B_final(i,1))<=5)
%         B_fit=[B_fit;[B_final(i,1),B_final(i,2)]];
%     end
% end
% % 对 A_fit 进行线性拟合
% coeff_A = polyfit(A_fit(:, 1), A_fit(:, 2), 1);
% slope_A = coeff_A(1);     % 斜率
% intercept_A = coeff_A(2); % 截距
% 
% % 对 B_fit 进行线性拟合
% coeff_B = polyfit(B_fit(:, 1), B_fit(:, 2), 1);
% slope_B = coeff_B(1);     % 斜率
% intercept_B = coeff_B(2); % 截距
% 
% % 绘制原始图片
% figure;
% imshow(normalizedImg);
% hold on;
% 
% % 绘制 A_final 的拟合直线用y算x
% yA = min(A_final(:, 2)):max(A_final(:, 2));
% xA = (yA-intercept_A)/slope_A;
% plot(A_Line_x, A_Line_y, 'r', 'LineWidth', 2);  % 以红色线条绘制直线
% 
% % 绘制 B_final 的拟合直线y算x
% yB = min(B_final(:, 2)):max(B_final(:, 2));
% xB = (yB-intercept_B)/slope_B;
% plot(B_Line_x, B_Line_y, 'g', 'LineWidth', 2);  % 以绿色线条绘制直线
% 
% hold off;