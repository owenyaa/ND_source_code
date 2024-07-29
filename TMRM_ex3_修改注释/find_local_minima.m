% function [Xmin] = getmin(y)
% a=length(y);
% j=1;
% for i=floor((a-1)/2):a
%     b=(y(i,1)-y(i-2,1))/2;
%     c=(y(i,1)+y(i-2,1))/2-y(i-1,1);
%     if y(i-2,1)>=y(i-1,1)&&y(i-1,1)<=y(i,1)&&c==0
%         Xmin(j)=y(i-1,1);
%         j=j+1;
%     elseif y(i-2,1)>=y(i-1,1)&&y(i-1,1)<=y(i,1)
%         Xmin(j)=y(i-1,1)+b^2/(4*c);
%         j=j+1;
%     end
% end
function [minimum_values] = find_local_minima(data)
% Find local minima in the input data vector.

% Get the length of the input data
data_length = length(data);

% Initialize an index variable and an array to store the minimum values
index = 1;
minimum_values = [];

% Iterate through the data to find local minima
for i = floor((data_length - 1) / 2):data_length
    % Calculate coefficients for a quadratic equation
    delta_x = (data(i, 1) - data(i - 2, 1)) / 2;
    delta_y = (data(i, 1) + data(i - 2, 1)) / 2 - data(i - 1, 1);
    
    % Check if the current point is a local minimum
    if data(i - 2, 1) >= data(i - 1, 1) && data(i - 1, 1) <= data(i, 1)
        if delta_y == 0
            % If the slope is zero, consider it a local minimum
            minimum_values(index) = data(i - 1, 1);
            index = index + 1;
        else
            % Calculate the x-coordinate of the local minimum using a quadratic equation
            minimum_values(index) = data(i - 1, 1) + delta_x^2 / (4 * delta_y);
            index = index + 1;
        end
    end
end
end
