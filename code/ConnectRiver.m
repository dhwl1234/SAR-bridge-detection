function [final_image] = ConnectRiver(bw_thinned)
%% 链接桥梁的两端
% 寻找左边最后一个点
final_image = bw_thinned;
left_lastRow = 0;
left_lastCol = 0;
foundOne = false;
for j=1:size(bw_thinned,2)
     % 从左往右逐列查找
    if ~foundOne && any(bw_thinned(:, j))
        % 遇到1后开始记录
        foundOne = true;
    end
    if foundOne
        if sum(bw_thinned(:,j))>0
            % 当前列存在非零元素
            left_lastRow = find(bw_thinned(:, j), 1, 'last');
            left_lastCol = j;
        else
            % 当前列全为零，结束循环
            break;
        end
    end
end
%寻找右边对应的点
right_firsttRow = 0;
right_firstCol = 0;
right_downRow = 0;
foundOne = false;
for j=left_lastCol+1:size(bw_thinned,2)
     % 从左往右逐列查找
    if ~foundOne && any(bw_thinned(left_lastRow, j))
        % 遇到1后开始记录
        foundOne = true;
    end
    if foundOne
       right_firsttRow=left_lastRow;
       right_firstCol = j;
       break;
    end
end
%以右侧对应点往上下两侧搜索，每次的步长一样，列坐标更新为最大的（确保一直往右搜索）
% 当步长扩大到与相同列坐标下与其他点相交时，记录下差值，开始向右遍历，则一直向右遍历
% 直到差值为0，此时点为目标点

length=0;
for i=3:200
    if bw_thinned(right_firsttRow+i,right_firstCol)==1
        length = i;
        right_downRow = right_firsttRow+i;
            break;
    elseif bw_thinned(right_firsttRow-i,right_firstCol)==1
       length = i;      
       right_downRow = right_firsttRow-i;
           break;
    end
end
%往右移动
remp_right_firsttRow=right_firsttRow ;
temp_right_firstCol=right_firstCol ;
temp_right_downRow=right_downRow;
foundOne = false;
for j = 1:size(bw_thinned,2)
    % 求上点的位置
    for i=1:20
        if bw_thinned(right_firsttRow+i,right_firstCol+j)==1
            remp_right_firsttRow = right_firsttRow+i;
            temp_right_firstCol = right_firstCol+j;
            foundOne = true;
        end
        if  foundOne 
            if bw_thinned(right_firsttRow+i,right_firstCol+j)==0
                right_firsttRow =remp_right_firsttRow;

                foundOne = false;
                break;
            end
        end
    end
    %求下点的位置
    for i=1:20
        if bw_thinned(right_downRow-i,right_firstCol+j)==1
            temp_right_downRow = right_downRow-i;
            foundOne = true;
        end
        if foundOne 
            if bw_thinned(right_downRow-i,right_firstCol+j)==1%保持在最底端
                right_downRow = temp_right_downRow;
                foundOne = false;
                break;
            end
        end
    end
    %求距离
    length = right_downRow-right_firsttRow;
    if length ==0
        right_firstCol =temp_right_firstCol;
        break;
    end
end
%% 进行主轴连接
% 定义两个点的坐标
x1 = left_lastCol;
y1 = left_lastRow;
x2 = right_firstCol;
y2 = right_firsttRow;

% 计算直线的斜率和截距
m = (y2 - y1) / (x2 - x1);
c = y1 - m * x1;

% 对于图像中的每个 x 坐标（从 x1 到 x2），计算对应的 y 坐标
x = x1:x2;
y = round(m * x + c); % 使用 round 函数进行舍入

% 设置对应像素为 1
for i = 1:numel(x)
    final_image(y(i), x(i)) = 1;
end


end