%SLMPC high way scenario following a high way 
clf;
clear;

%the overall delta time for each step and the obstacle constant velocity is
%0.5
Ts = 0.1; %time for each step
v_obstacle = 0.02;
counter_1 = 1; %counter for to plot the obstacle vehicle
counter_2 = 1; %the reference counter for obtaining the reference values

%the below code is to obtain a rectangle representing a car driving along a
%line
x_ob = 0.1;
y_ob = 4;
theta_ob = 0;
s_ob = 2;
% [vx, vy] = create_rectangle(x_ob,y_ob,theta_ob,s_ob);
% fill(vx, vy, [0.25,0.41,0.88])
% axis equal
% hold on

x_map = [0:0.1:30];
num = 4;
y_map = num*ones(size(x_map));
plot(x_map, y_map, 'g--')
axis equal
hold on

% while x_ob <= x_map(1,end)    
%     x_ob = x_ob + Ts*v_obstacle;
%     y_ob = 4;
%     theta_ob = 0;
%     s_ob = 2.5;
%     if mod(counter_1,3000) == 0
%         [vx, vy] = create_rectangle(x_ob,y_ob,theta_ob,s_ob);
%         fill(vx, vy, [0.25,0.41,0.88])
%         axis equal
%         hold on
%     end
%     if mod(counter_1, 3001) == 0
%         fill(vx, vy, [1,1,1])
%         axis equal
%         hold on
%         length_ob = vx(end) - vx(1);
%         width_ob = vy(end) - vy(1);
%         rectangle('position',[vx(2),vy(2),length_ob,width_ob],'EdgeColor','b','LineStyle','--');
%         hold on
%     end
%     counter_1 = counter_1 + 1; 
% end

%the specifications of the jackel robot
length_front_axle = 0.150; %this is given in meters
length_rear_axle = 0.150;  %this is given in meters
mass_jackel = 17;   %this is given as 17kg
radius_of_jackel = 0.125; %in meters
length_of_jackel = 0.51; %in meters
yaw_inertia_jackel = (1/4)*(mass_jackel)*(radius_of_jackel^2) + (1/12)*(mass_jackel)*(length_of_jackel^2);
%constraints for the Jackel
heading_max = pi/6;
heading_min = -pi/6;
a_max = 0.5;
a_min = -0.5;
x_min = -2;
x_max = 105;
y_min = -20;
y_max = 20;
eta_min = -pi;
eta_max = pi;
vel_max = 0.1;
vel_min = -0.1;

%system parameters
n_prediction = 10; %prediction window size
n_x = 4; %number of states
n_u = 2; %number of inputs
n_c = 4; %number of active constraints

x_con_max = [x_max; y_max; eta_max; vel_max];
x_con_min = [x_min; y_min; eta_min; vel_min];
u_constraint_high = [heading_max;a_max];
u_constraint_low = [heading_min;a_min];

%input constraints are shown as upper and lower matrices for the
%quadratic programming
[lb_matrix, ub_matrix] = SLMPC_functions.input_constraints(u_constraint_high,u_constraint_low,n_u,n_prediction);
%state constraints formation
[x_constraint_high, x_constraint_low] =SLMPC_functions.formulate_constraints(x_con_max,x_con_min,n_c,n_prediction); 


%initial values of the robot
x_robot = 1;
y_robot = 0;
eta_robot = 1e-3;
robot_vel = 1e-3;
x_state = [x_robot;y_robot;eta_robot;robot_vel]; %the states of the robot
u_in = [0;0]; %the input of the robot with delta and acceleration

%kinematic equations for the bicycle model of the jackel robot and assume it is
%frontal drive, eta is the inertial heading, which is the angle between the
%x-axis and the center of gravity of the robot, there is to input values, a
%and heading_ang
C_lin = eye(4);  %the output C matrix is identity
%the nonlin_eq(x,u) function giving in the form dx/dt that is the model of the  
sys.lf = length_front_axle;
sys.lr = length_rear_axle; 
model = @(x,delta,a) SLMPC_functions.nonlin_eq(x,delta,a,sys); 
Cd = C_lin;
state_matrix = [x_robot;y_robot;eta_robot;robot_vel]; %use this matrix to store each state value at each time index k


%the weighting matrices and the set up for the MPC controller 
weight_error_x = 1;
weight_error_y = 1;
weight_error_eta = 1;
weight_error_vel = 1;
weight_input_heading = 1;
weight_input_a = 1;

q_small = [weight_error_x 0 0 0; 0 weight_error_y 0 0; 0 0 weight_error_eta 0; 0 0 0 weight_error_vel];
r_small = [weight_input_heading 0; 0 weight_input_a];

state_weights = [weight_error_x weight_error_y weight_error_eta weight_error_vel];
input_weights = [weight_input_heading weight_input_a];
Q = diag(repmat(state_weights, 1, n_prediction)); 
R = diag(repmat(input_weights, 1, n_prediction));

[vx, vy] = create_rectangle(x_robot,y_robot,eta_robot,0.5);
fill(vx, vy, 'r')
axis equal
hold on
drawnow

%symbolic function used to do the derivative of the B matrix for the

%linearized model of the robot
%[df_B_1,df_B_2,df_B_3] = SLMPC_functions.linear_B_matrix();

% Setting the quadprog with 200 iterations at maximum
opts = optimoptions('quadprog', 'MaxIter', n_prediction, 'Display','off');
tic;
while x_robot*1.05 <= x_map(1,end) 
    if counter_2 == 100
        j = counter_2
    end
    %obtain the references on the path 
    [x_ref,y_ref,eta_ref,index_ref] = SLMPC_functions.obtain_reference(x_robot,y_robot,eta_robot,x_map,y_map,counter_2); %for only constant line function
    plot(x_ref,y_ref,'r*');
    axis equal
    hold on
    vel_ref = 0;
    if index_ref + n_prediction >= length(x_map)
        r_val = SLMPC_functions.arrange_ref(x_map,y_ref,eta_ref,vel_ref,length(x_map)-n_prediction,n_prediction,n_x,n_u);
    else
        r_val = SLMPC_functions.arrange_ref(x_map,y_ref,eta_ref,vel_ref,index_ref,n_prediction,n_x,n_u);
    end
    
    %r_val = [x_ref;y_ref;eta_ref;0];
    [x_state, dx_state] = SLMPC_functions.RK4(x_state,u_in,Ts,model);     % i.e. simulate one step forward with Runge-Kutta 4 order integrator
    %linearize the bicycle model
    [A_lin,B_lin,K_lin] = SLMPC_functions.linearize_model(x_state,dx_state,u_in,sys);
    A_lin = double(A_lin);
    B_lin = double(B_lin);
    K_lin = double(K_lin);
    [Ad,Bd,Kd] = SLMPC_functions.discretize(A_lin,B_lin,K_lin,Ts);
    [G, f] = SLMPC_functions.grad_n_hess(R, Q, Ad, Bd, Cd, [], Kd, r_val, n_prediction, x_state);
    [Px,Hx,Ex,V] = SLMPC_functions.state_constraints(Ad, Bd, Kd, n_prediction, n_x, n_u, x_state,n_c);
    %state constraints below are shown up as the inequality constraints for
    %the quadratic programming
    Aineq = [V*Hx;-V*Hx];
    bineq = [x_constraint_high - V*Px*x_state-V*Ex*Kd;-x_constraint_low+V*Px*x_state+V*Ex*Kd];
    u_input = quadprog(G, f, Aineq, bineq, [], [], lb_matrix, ub_matrix, [],opts);
    heading_input = u_input(1);
    a_input = u_input(2);
    
    %update of the robot state each time
    robot_vel = robot_vel + a_input*Ts;
    beta_robot = atan((length_rear_axle/(length_front_axle+length_rear_axle))*tan(heading_input));
    eta_robot = eta_robot + (robot_vel/length_rear_axle)*sin(beta_robot)*Ts;
    y_robot = y_robot + robot_vel*sin(eta_robot+beta_robot)*Ts;
    x_robot = x_robot + robot_vel*cos(eta_robot+beta_robot)*Ts;
    state_matrix = [state_matrix [x_robot;y_robot;eta_robot;robot_vel]];
    if x_robot == 0 
        x_robot = 1e-3;
    end
    if y_robot == 0 
        y_robot = 1e-3;
    end
    if eta_robot == 0 
        eta_robot = 1e-3;
    end
    if robot_vel == 0 
        robot_vel = 1e-3;
    end
    x_state = [x_robot;y_robot;eta_robot;robot_vel];
    
    if mod(counter_2,20) == 0
        [vx, vy] = create_rectangle(x_robot,y_robot,eta_robot,0.5);
        fill(vx, vy, 'r')
        axis equal
        hold on
        drawnow
    end
    counter_2 = counter_2 + 1;
    toc
end

plot(x_map, y_map, 'g--')
axis equal
hold on