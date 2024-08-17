function output = createModel(radius, center, J)

    % Create equally spaced vector of values between 0 and 2π 
    theta = linspace(0, 2*pi);
    
    % Calculate the x and y coordinates of the circle
    x = center(1) + radius * cos(theta);
    y = center(2) + radius * sin(theta);
    
    %*************************USER POINTS******************************
    % Generate random angles within the circle's range for users
    theta = 2*pi*rand(3, 1);
    
    % Generate random radii within the circle's range for users
    r = radius * sqrt(rand(3, 1));
    
    % Convert polar coordinates to Cartesian coordinates for users
    uX = center(1) + r .* cos(theta);
    uY = center(2) + r .* sin(theta);
    
    userDist = sqrt(uX.^2 + uY.^2);
    
    %*************************CLUTTER POINTS***************************
    % Generate random angles within the circle's range for users
    theta = 2*pi*rand(J, 1);
    
    % Generate random radii within the circle's range for users
    r = radius * sqrt(rand(J, 1));
    
    % Convert polar coordinates to Cartesian coordinates for users
    cX = center(1) + r .* cos(theta);
    cY = center(2) + r .* sin(theta);
    
    clutterDist = sqrt(cX.^2 + cY.^2);
    
    %*************************PLOTTING GRAPH***************************
    % Plotting the circle
    % plot(x, y);
    % %Adding points to figure
    % hold on; 
    % plot(0, 0, "^"); %plotting the BS
    % plot(uX, uY, "square"); %plotting the 3 users
    % plot(cX, cY, "o"); %plotting the 3 users
    % hold off; 

    combinedVector = vertcat(userDist, clutterDist); %concatenate vectors vertically
    output = combinedVector;
end