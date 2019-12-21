%Berken Utku Demirel - 2166221
clear 
clc
close all
%% Loading image
addpath('C:\Users\berkenutku\Desktop\415_HW\P3_Sample Images');
[file,~,path] = uigetfile('*.mat','Workspace File');
[~,name,~] = fileparts(file);
data_struct = load(file);
data = data_struct.(name);
%% Taking input from users
Projection_angle_step_size = input("Enter projection angle step size: ");
Number_of_sampling_point = input("Enter the number of sampling points for each projection: ");
%% Calculate size of image
Size_of_image = size(data);
M = Size_of_image(1,1);
%% Calculating the t and theta values according to user inputs
t_value_first = -M * sqrt(2) / 2;
t_value_last = t_value_first * -1;
t_values = linspace( t_value_first,t_value_last, Number_of_sampling_point);

steps_for_theta = 0:floor(180/Projection_angle_step_size);
theta_values = Projection_angle_step_size .* steps_for_theta;

%% Calculatin relavant intersection points
projection_matrix = find_projections(theta_values,t_values,M,data);
%% Filtering

data_to_filter = projection_matrix.';

L = Number_of_sampling_point;

Triangle = triang(L);

fft_of_projection_data = fft2(data_to_filter);

result_of_multiplication = Triangle .* fft_of_projection_data.';

filtering_output = ifft2(result_of_multiplication);
%% Back Projection

back_projected_data = filtering_output;

back_projection_matrix_with_filter = find_back_projections(theta_values, t_values, M, filtering_output.');
back_projection_matrix = find_back_projections(theta_values, t_values, M, data_to_filter);

function y = find_address_distance(x, y, M, projection_matrix_value,back_projection_matrix)
points= [x.' y.'];
points = sortrows(points);
    for i=1:length(x)-1
        dist = sqrt( (points(i+1,1) - points(i,1))^2 + (points(i+1,2) - points(i,2))^2);
        mid_x = (points(i,1) + points(i+1,1)) / 2;
        mid_y = (points(i,2) + points(i+1,2)) / 2;
        rowdata = (M/2) - floor(mid_y);
        columndata = (M/2) + ceil(mid_x);
        back_projection_matrix(rowdata,columndata) = back_projection_matrix(rowdata,columndata) + (projection_matrix_value * dist); 
    end
    y = back_projection_matrix;
end

function y = find_back_projections(theta_values, t_values, M, projection_matrix)
X = linspace(-M/2,M/2,M+1); % X points in the border of image
Y = linspace(-M/2,M/2,M+1); % Y points in the border of image

back_projection_matrix = zeros(M, M);

theta_count = 1; % Counter for theta values
for theta = theta_values
    coeff_cos = cosd(theta); % Calculating cos values with respect to theta
    coeff_sin = sind(theta); % Calculating sin values with respect to theta
    t_count = 1; % Counter for t values
    for t = t_values
        bool_break = 1;
        if(coeff_sin == 0) % If projection angle is 0 degree
            if(abs(t) > M/2 || t == 0)
                p = 0;
                t_count = t_count + 1;
                bool_break = 0;
            else
                x(1,1:length(X)) = t;
                y = Y;
                back_projection_matrix = find_address_distance(x,y,M,projection_matrix(theta_count,t_count),back_projection_matrix);
                t_count = t_count + 1;
                bool_break = 0;
            end
        elseif(coeff_cos == 0) % If the projection angle is 90 degree
            if(abs(t) > M/2 || t == 0)
                p = 0;
                t_count = t_count + 1;
                bool_break = 0;
            else
                y(1,1:length(X)) = t;
                x = X;
                back_projection_matrix = find_address_distance(x,y,M,projection_matrix(theta_count,t_count),back_projection_matrix);
                t_count = t_count + 1;
                bool_break = 0;
            end
        else % Projection angles different than 0 and 90 degrees.
            x = ((t - Y * coeff_sin) / coeff_cos);
            y = ((t - X * coeff_cos) / coeff_sin);
        end
        if(bool_break)
            check_matrix_x = [x' Y.'];
            check_matrix_y = [y.' X.'];
            check_matrix_x_bool = abs(check_matrix_x(:,1)) <= M/2;
            bool_index_x = find(check_matrix_x_bool == 1);
            check_matrix_x = check_matrix_x(bool_index_x,:);
            check_matrix_y_bool = abs(check_matrix_y(:,1)) <= M/2;
            bool_index_y = find(check_matrix_y_bool == 1);
            check_matrix_y = check_matrix_y(bool_index_y,:);
            check_matrix = [check_matrix_x ; check_matrix_y(:,[2 1])];
            unique_matrix = unique(check_matrix,'rows');
            back_projection_matrix = find_address_distance(unique_matrix(:,1).', unique_matrix(:,2).', M,projection_matrix(theta_count,t_count),back_projection_matrix);
            t_count = t_count + 1;
        end
    end
    theta_count = theta_count + 1;
end

y = back_projection_matrix;
end

% Function to detect addresses and sum product of pixel value with distance
function y = find_attenu(x, y, M, data)
        points= [x.' y.'];
        points = sortrows(points);
        count = 0;
        for i=1:length(x)-1
            dist = sqrt( (points(i+1,1) - points(i,1))^2 + (points(i+1,2) - points(i,2))^2);
            mid_x = (points(i,1) + points(i+1,1)) / 2;
            mid_y = (points(i,2) + points(i+1,2)) / 2;
            rowdata = (M/2) - floor(mid_y);
            columndata = (M/2) + ceil(mid_x);
            count = count + data(rowdata,columndata) * dist;
        end
        y = count;
end


function y = find_projections(theta_values, t_values, M, data)
X = linspace(-M/2,M/2,M+1); % X points in the border of image
Y = linspace(-M/2,M/2,M+1); % Y points in the border of image

projection_matrix = zeros(length(t_values), length(theta_values));

theta_count = 1; % Counter for theta values
for theta = theta_values
    coeff_cos = cosd(theta); % Calculating cos values with respect to theta
    coeff_sin = sind(theta); % Calculating sin values with respect to theta
    t_count = 1; % Counter for t values
    p = zeros([1 length(t_values)]); % Projection vector for corresnponding theta value
    for t = t_values
        bool_break = 1;
        if(coeff_sin == 0) % If projection angle is 0 degree
            if(abs(t) > M/2 || t == 0)
                p(t_count) = 0;
                t_count = t_count + 1;
                bool_break = 0;
            else
                x(1,1:length(X)) = t;
                y = Y;
                p(t_count) = find_attenu(x,y,M,data);
                t_count = t_count + 1;
                bool_break = 0;
            end
        elseif(coeff_cos == 0) % If the projection angle is 90 degree
            if(abs(t) > M/2 || t == 0)
                p(t_count) = 0;
                t_count = t_count + 1;
                bool_break = 0;
            else
                y(1,1:length(X)) = t;
                x = X;
                p(t_count) = find_attenu(x,y,M,data);
                t_count = t_count + 1;
                bool_break = 0;
            end
        else % Projection angles different than 0 and 90 degrees.
            x = ((t - Y * coeff_sin) / coeff_cos);
            y = ((t - X * coeff_cos) / coeff_sin);
        end
        if(bool_break)
            check_matrix_x = [x' Y.'];
            check_matrix_y = [y.' X.'];
            check_matrix_x_bool = abs(check_matrix_x(:,1)) <= M/2;
            bool_index_x = find(check_matrix_x_bool == 1);
            check_matrix_x = check_matrix_x(bool_index_x,:);
            check_matrix_y_bool = abs(check_matrix_y(:,1)) <= M/2;
            bool_index_y = find(check_matrix_y_bool == 1);
            check_matrix_y = check_matrix_y(bool_index_y,:);
            check_matrix = [check_matrix_x ; check_matrix_y(:,[2 1])];
            unique_matrix = unique(check_matrix,'rows');
            p(t_count) = find_attenu(unique_matrix(:,1).', unique_matrix(:,2).', M, data);
            t_count = t_count + 1;
        end
    end
    projection_matrix(:,theta_count) = p;
    theta_count = theta_count + 1;
end

y = projection_matrix;

end

