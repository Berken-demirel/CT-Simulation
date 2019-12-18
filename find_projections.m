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


