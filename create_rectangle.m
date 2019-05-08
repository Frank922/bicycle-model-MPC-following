function [vx,vy] = create_rectangle(x,y,theta,s)
    % s is the length of the vehicle
    %the width of the vehicle is 0.75s
    x_d_front = x + s*cos(theta);
    y_d_front = y + s*sin(theta);
    x_d_back = x - s*cos(theta);
    y_d_back = y - s*sin(theta);
    %m = (y_d_front - y)/(x_d_front - x);
    %m_perp = -1/m; %in perpendicular line equation 
    width = 0.75*s;
    angle_val = pi/2 - theta; 
    len_perp = 0.5*width;
    vx(1,1) = x_d_back - len_perp*cos(angle_val);
    vy(1,1) = y_d_back + len_perp*sin(angle_val);
    vx(1,2) = x_d_back + len_perp*cos(angle_val);
    vy(1,2) = y_d_back - len_perp*sin(angle_val);
    vx(1,3) = x_d_front + len_perp*sin(theta);
    vy(1,3) = y_d_front - len_perp*cos(theta);
    vx(1,4) = x_d_front - len_perp*cos(angle_val);
    vy(1,4) = y_d_front + len_perp*sin(angle_val);
end
    
    
    
    