function [indicatorX, indicatorY] = ddplot(X,Y)
%DDPLOT Summary of this function goes here
%   Detailed explanation goes here
Z = [X
    Y];

[szX1, ~] = size(X);
[szY1, ~] = size(Y);

indicatorX = zeros (szX1+szY1);
indicatorY = zeros (szX1+szY1);

for i=1:szX1+szY1
    for j = 1:szX1
        indicatorX(i) =  indicatorX(i) + norm(Z(i) - X(j));
    end
    for j = 1:szX1
        indicatorY(i) = indicatorY(i) + norm(Z(i) - Y(j));
    end

    plot(indicatorX,indicatorY, 'o')
end

