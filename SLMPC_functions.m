classdef SLMPC_functions
    methods(Static)
        function [x_ref,y_ref,eta_ref,index_ref] = obtain_reference(x_robot,y_robot,eta_robot,x_map,y_map, counter_2)
            I = find(x_map >= x_robot + 0.5);
            index_ref = I(1);
            x_ref = x_map(I(1));
            y_ref = y_map(1,end);
            eta_ref = atan((y_map(1,end)-y_map(1,1))/(x_map(1,end)-x_map(1,1)));
        end
        function r_val = arrange_ref(x_map,y_ref,eta_ref,vel_ref,index_ref,np,nx,nu)
            r_val = zeros(np*nx,1);
            r_index = 1;
            for i = index_ref:index_ref+np-1
                r_val(r_index,1) = x_map(i);
                r_index = r_index + 1;
                r_val(r_index,1) = y_ref;
                r_index = r_index + 1;
                r_val(r_index,1) = eta_ref;
                r_index = r_index + 1;
                r_val(r_index,1) = vel_ref;
                r_index = r_index + 1;
            end
        end
        function dx = nonlin_eq(x,a,delta,sys)
            %This function represents nonlinear model of simple mathematical 
            lf = sys.lf;
            lr = sys.lr;
            beta = atan((lr/(lf+lr))*tan(delta));
            dx =zeros(length(x),1); 
            dx(1) = x(4)*cos(x(3)+beta);
            dx(2) = x(4)*sin(x(3)+beta);
            dx(3) = (x(4)/lr)*sin(beta);
            dx(4) = a;
        end
        function [x_state, dx_state] = RK4(x1,u,Ts,k)
        % This function calculates next position, speed in small time Ts 
        % with initial position x and input u,which is considered to be constant in 
        % time Ts. In order to find position RK4 integration method is used 
        
        %u(1) is the input for the heading angle of the vehicle and u(2) is
        %the longitudinal acceleration 
                %Algorithm of RK4
            k1 = k(x1,u(1),u(2)); 
            k2 = k(x1 + (Ts/2)*k1,u(1),u(2));         
            k3 = k(x1 + (Ts/2)*k2,u(1),u(2));        
            k4 = k(x1+Ts*k3,u(1),u(2));        
            dx_state=k1;
            x_state = x1 + Ts*(k1 + 2*k2 + 2*k3 + k4)/6;
        end
        function [df_B_1,df_B_2,df_B_3] = linear_B_matrix()
            syms f(x) v eta lr lf
            f(x) = v*cos(eta+atan((lr/(lf + lr))*tan(x)));
            df_B_1 = diff(f,x);
            f(x) = v*sin(eta+atan((lr/(lf + lr))*tan(x)));
            df_B_2 = diff(f,x);
            f(x) = (v/lr)*sin(atan((lr/(lf + lr))*tan(x)));
            df_B_3 = diff(f,x);
        end
        function [A_lin,B_lin,K_lin] = linearize_model(x_state,dx_state,u_in,sys)
            %the below functions are used to directly sub with values
            %instead of using subs functions to increase the speed. 
            %x_state = [x_robot,y_robot,eta_robot,robot_vel]
            lf_val = sys.lf;
            lr_val = sys.lr;
            B_lin_1 = -(lr_val*x_state(4)*sin(x_state(3) + atan((lr_val*tan(u_in(1)))/(lf_val + lr_val)))*(tan(u_in(1))^2 + 1))/((lf_val + lr_val)*((lr_val^2*tan(u_in(1))^2)/(lf_val + lr_val)^2 + 1));
            B_lin_2 = (lr_val*x_state(4)*cos(x_state(3) + atan((lr_val*tan(u_in(1)))/(lf_val + lr_val)))*(tan(u_in(1))^2 + 1))/((lf_val + lr_val)*((lr_val^2*tan(u_in(1))^2)/(lf_val + lr_val)^2 + 1)); 
            B_lin_3 = (x_state(4)*(tan(u_in(1))^2 + 1))/((lf_val + lr_val)*((lr_val^2*tan(u_in(1))^2)/(lf_val + lr_val)^2 + 1)^(1/2)) - (lr_val^2*x_state(4)*tan(u_in(1))^2*(tan(u_in(1))^2 + 1))/((lf_val + lr_val)^3*((lr_val^2*tan(u_in(1))^2)/(lf_val + lr_val)^2 + 1)^(3/2));
            beta = atan((lr_val/(lf_val+lr_val))*tan(u_in(1)));
            A_lin = [0 0 -x_state(4)*sin(x_state(3)+beta) cos(x_state(3)+beta); 0 0 x_state(4)*cos(x_state(3)+beta) sin(x_state(3)+beta); 0 0 0 sin(beta)/lr_val; 0 0 0 0];
            B_lin = [B_lin_1 0; B_lin_2 0; B_lin_3 0; 0 1];
            %Kc is offcet created by linear approximation of the system and initial
            %parameters.
            K_lin=dx_state-A_lin*x_state-B_lin*u_in;
        end
        function [Ad,Bd,Kd] = discretize(A_lin,B_lin,K_lin,Ts)
            %Applies discretization on continious system of the form
            %dx/dt = Ac*x +Bc*u+Kc and calculates such A,B,K that provides
            %analogical discrete system of form x(k+1)=A*x(k)+B*u(k)+K
            %more info in https://en.wikipedia.org/wiki/Discretization 
            mat_val = [A_lin B_lin; [0 0 0 0; 0 0 0 0] [0 0; 0 0]]*Ts;
            matrix_AB = expm(mat_val);                      
            Ad = matrix_AB(1:4,1:4);
            Bd = matrix_AB(1:4,5:6);
            %approximation of the discrete systems: dx = A_Lx+B_Lu+K_op and
            %this can be shown as x[k+1] = (I+A_L*T)x+(B_L*T)u+k*T
            Kd = K_lin*Ts;
        end
        function [G, f] = grad_n_hess(R, Q, Ad, Bd, C, D, Kd, rr, np, x)
            %This function returns gradient and hessian matrices of the cost
            %function using function calc_hp to obtain horizon matrices.
            [Hx,Px,Km] = SLMPC_functions.calc_hp(Ad,Bd,C,D,np);
            K1 = (Km*Kd)-rr;
            G = Hx'*Q*Hx + R;
            f = Hx'*Q*(Px*x+K1);
        end
        function [Hx,Px,Km] = calc_hp(A,B,C,D,np)
            %This function calculates prediction matrices for vector x and output
            %vector y with prediction horizon np
    
            %Initialization
    
            % number of states
            nx = size(A, 1);
            %number of inputes 
            nu = size(B, 2);
            %number of outputs 
            no=size(C,1);

            %zero initialization  
            Px=zeros(np*nx,size(A,2));
            Hx=zeros(np*size(B));
            H_x=zeros(np*size(B));
            Km=zeros(no*np,size(C,2));
            S=zeros(size(C));
            
            %start of the main loop 
            for ind1=1:np

                % Filling Matrices Px,P,Km recursively 
                Px((1+(ind1-1)*nx):ind1*nx,1:nx)=A^(ind1-1);
                Km(1+no*(ind1-1):no*ind1,:)=S;
                S=S+C*A^(ind1-1);

                %Filling Marices Hx, H recurcively  
                for ind2=1:np            
                    if(ind1>=ind2)                
                        H_x((1+(ind1-1)*nx):(ind1)*nx,(1+(ind2-1)*nu):(ind2)*nu)=A^(ind1-ind2)*B;                
                    end
                end
            %End of the main loop    
            end
        
        %fix Hx 
        
        j = 1;
        i = 1;
        while j <= np*length(B)
            if j ~= 1
                Hx(j:j+(length(B)-1),:) = H_x(i:i+(length(B)-1),:);
                i = i + length(B);
            end
            j = j + length(B);
        end
        %End of function  
        end
        function [Px,Hx,Ex,V] = state_constraints(Ad, Bd, Kd, np, nx, nu, xd, nc)
            %to find the state constraints x[k] = Px*xd + Hx*u[k] + Ex*Kd
            Px = zeros((np+1)*nx,nx);
            Hx = zeros((np+1)*nx,np*nu);
            Ex = zeros((np+1)*nx,nx);
            
            %create the Px matrix 
            i = 1;
            power_val = 1;
            while i <= (np+1)*nx
                if i == 1
                    Px(i:i+nx-1,:) = eye(nx);
                else
                    Px(i:i+nx-1,:) = Ad^(power_val);
                    power_val = power_val+1;
                end
                i = i + nx;
            end
            %create the Hx matrix 
            i = 1; 
            for j = 1:np+1
                if j == 1
                    Hx(i:i+nx-1,:) = zeros(nx,np*nu);
                else
                    counter = 1;
                    Hx_buff = Bd;
                    while counter <= j - 2 
                        Hx_buff = [(Ad^counter)*Bd Hx_buff];
                        counter = counter + 1;
                    end
                    Hx(i:i+nx-1,:) = [Hx_buff zeros(nx,np*nu-size(Hx_buff,2))];
                end
                i = i + nx;
            end
            %create the Ex matrix 
            S = zeros(nx,nx);
            i = 1;
            for  j = 1:np+1
                if j == 1
                    Ex(i:i+nx-1,:) = S;
                elseif j == 2
                    S = S + eye(nx);
                    Ex(i:i+nx-1,:) = S;
                else
                    S = S + Ad^(j-2);
                    Ex(i:i+nx-1,:) = S;
                end
                i = i + nx;
            end
            
            %for bicycle model only
            %the below part is for calculating V without obstacle avoidance
            %and pure line following
            V = eye((np+1)*nc);
        end
        
        function [lb_matrix, ub_matrix] = input_constraints(u_constraint_high,u_constraint_low,nu,np)
            lb_matrix = zeros(nu*np,1);
            ub_matrix = zeros(nu*np,1);
            j = 1;
            while j < nu*np
                lb_matrix(j:j+(nu-1),:) = u_constraint_low;
                ub_matrix(j:j+(nu-1),:) = u_constraint_high;
                j = j + nu;
            end
        end
        
        function [x_constraint_high, x_constraint_low] = formulate_constraints(x_con_max,x_con_min,nc,np)
            x_constraint_high = zeros((np+1)*nc,1);
            x_constraint_low = zeros((np+1)*nc,1);
            j = 1;
            while j < nc*(np+1)
                x_constraint_high(j:j+(nc-1),:) = x_con_max;
                x_constraint_low(j:j+(nc-1),:) = x_con_min;
                j = j + nc;
            end
        end
        
        %the below functions are used to bypass an obstacle in a high way
        %scenario
        function [x_constraint_high, x_constraint_low] = modify_state_constraints(x_constraint_high, x_constraint_low,x_ob_high,x_ob_low,y_ob_high,y_ob_low, nc)
            index = 1;
            while index <= length(x_constraint_high)
                for j = 1:nc+2
                    if j == 2
                        x_constraint_high(index+j-1,1) = x_ob_low;
                        x_constraint_low(index+j-1,1) = x_ob_high;
                    elseif j == 4 
                        x_constraint_high(index+j-1,1) = y_ob_low;
                        x_constraint_low(index+j-1,1) = y_ob_high;
                    end
                end
                index = index + (nc+2);
            end
        end
        
        function [a_matrix, b_constraint,x_limit,y_limit] = find_constraint_up(ref_robot,x_robot,y_robot,robot_front_up,robot_front_down,robot_back_up,robot_back_down,ob_up_left,ob_up_right,ob_down_left,ob_down_right,back_mid)
            %now locate where the robot is wihtin respect to the obstacle
            if robot_back_down(2) < ob_up_left(2) && (robot_back_down(1) + 1.0) < ob_down_right(1)
                constraint_slope = (ob_up_left(2) - robot_front_down(2))/(ob_up_left(1) - (robot_front_down(1) + 0.5));
                constraint_b = robot_front_down(2) - constraint_slope*(robot_front_down(1)+0.5);
                if constraint_slope*robot_back_down(1) + constraint_b <= robot_back_down(2) %the bottom back corner of the robot is above the constraint vector line
                    a_matrix = [-constraint_slope ; 1];
                    b_constraint = constraint_b;
                    x_limit = [(robot_front_down(1)+0.5) ob_up_left(1)];
                    y_limit = [robot_front_down(2) ob_up_left(2)];
                else %the bottom back corner of the robot is below the constraint vector line
                    constraint_slope = (ob_up_left(2) - robot_back_down(2))/(ob_up_left(1) - (robot_back_down(1) + 0.5));
                    a_matrix = [-constraint_slope ; 1];
                    b_constraint = robot_back_down(2) - constraint_slope*(robot_back_down(1)+0.5); 
                    x_limit = [(robot_back_down(1) + 0.5) ob_up_left(1)];
                    y_limit = [robot_back_down(2) ob_up_left(2)];
                end
            elseif robot_back_down(2) >= ob_up_left(2) && (robot_back_down(1) + 1.0) < ob_down_right(1)
                a_matrix = [0 ; 1];
                b_constraint = y_robot;
                x_limit = [ob_up_left(1) ob_up_right(1)];
                y_limit = [ob_up_left(2) ob_up_right(2)];
            elseif robot_back_down(2) >= ob_up_left(2) && (robot_back_down(1) + 1.0) >= ob_down_right(1)
                constraint_slope = (ref_robot(2) - robot_back_down(2))/(ref_robot(1) - (robot_back_down(1) - 0.5));
                a_matrix = [-constraint_slope ; 1];
                b_constraint = robot_back_down(2) - constraint_slope*(robot_back_down(1)-0.5); 
                x_limit = [(robot_back_down(1) - 0.5) ref_robot(1)];
                y_limit = [robot_back_down(2) ref_robot(2)];
            end
        end
        
        %the below function is to modify the V vector for the constraint
        %avoidance
        function [Px,Hx,Ex,V_high,V_low] = state_constraints_modified(a_matrix,b_constraint,Ad, Bd, Kd, np, nx, nu, xd,nc)
            %to find the state constraints x[k] = Px*xd + Hx*u[k] + Ex*Kd
            Px = zeros((np+1)*nx,nx);
            Hx = zeros((np+1)*nx,np*nu);
            Ex = zeros((np+1)*nx,nx);
            
            %create the Px matrix 
            i = 1;
            power_val = 1;
            while i <= (np+1)*nx
                if i == 1
                    Px(i:i+nx-1,:) = eye(nx);
                else
                    Px(i:i+nx-1,:) = Ad^(power_val);
                    power_val = power_val+1;
                end
                i = i + nx;
            end
            %create the Hx matrix 
            i = 1; 
            for j = 1:np+1
                if j == 1
                    Hx(i:i+nx-1,:) = zeros(nx,np*nu);
                else
                    counter = 1;
                    Hx_buff = Bd;
                    while counter <= j - 2 
                        Hx_buff = [(Ad^counter)*Bd Hx_buff];
                        counter = counter + 1;
                    end
                    Hx(i:i+nx-1,:) = [Hx_buff zeros(nx,np*nu-size(Hx_buff,2))];
                end
                i = i + nx;
            end
            %create the Ex matrix 
            S = zeros(nx,nx);
            i = 1;
            for  j = 1:np+1
                if j == 1
                    Ex(i:i+nx-1,:) = S;
                elseif j == 2
                    S = S + eye(nx);
                    Ex(i:i+nx-1,:) = S;
                else
                    S = S + Ad^(j-2);
                    Ex(i:i+nx-1,:) = S;
                end
                i = i + nx;
            end
            
            %for bicycle model only
            %the below part is for calculating V without obstacle avoidance
            V_high = eye((np+1)*nc);
            %below is to obtain V_low matrix
            %v_small = [1 0 0 0 0; 0 1 0 0 0; -a_matrix(1,1) -a_matrix(2,1) 0 0 0; 0 0 1 0 0; 0 0 0 1 0];
            V_low = zeros((np+1)*(nc+1),(np+1)*nc); %allocate memory for V_low
            i = 1;
            j = 1;
            while j <= (np+1)*nc
                V_low(j,i) = 1;
                V_low(j+1,i+1) = 1;
                V_low(j+2,i+2) = 1;
                V_low(j+3,i+3) = 1;
                V_low(j+4,i) = a_matrix(1);
                V_low(j+4,i+1) = a_matrix(2);
                i = i + nc;
                j = j + (nc+1);
            end
        end
        function [x_constraint_high_modified, x_constraint_low_modified] = constraint_obstacle(x_constraint_high, x_constraint_low, a_matrix,  b_constraint, nc, np) 
            x_constraint_high_modified = x_constraint_high;
            x_con_min = [x_constraint_low(1); x_constraint_low(2); x_constraint_low(3); x_constraint_low(4); b_constraint];
            
            x_constraint_low_modified = zeros((np+1)*(nc+1),1);
            j = 1;
            while j <= length(x_constraint_low_modified)
                x_constraint_low_modified(j:j+nc,:) = x_con_min;
                j = j + (nc + 1);
            end
        end
        
        function [H,P,E] = calc_hp_HPE(A,B,C,D,np)
            %This function calculates prediction matrices for vector x and output
            %vector y with prediction horizon np
    
            %Initialization
    
            % number of states
            nx = size(A, 1);
            %number of inputes 
            nu = size(B, 2);
            %number of outputs 
            no=size(C,1);

            %zero initialization  
            P=zeros(np*nx,size(A,2));
            H=zeros(np*size(B));
            H_x=zeros(np*size(B));
            E=zeros(no*np,size(C,2));
            S=zeros(size(C));
            
            %start of the main loop 
            for ind1=1:np

                % Filling Matrices Px,P,Km recursively 
                P((1+(ind1-1)*nx):ind1*nx,1:nx)=A^(ind1-1);
                E(1+no*(ind1-1):no*ind1,:)=S;
                S=S+C*A^(ind1-1);

                %Filling Marices Hx, H recurcively  
                for ind2=1:np            
                    if(ind1>=ind2)                
                        H_x((1+(ind1-1)*nx):(ind1)*nx,(1+(ind2-1)*nu):(ind2)*nu)=A^(ind1-ind2)*B;                
                    end
                end
            %End of the main loop    
            end
        
        %fix Hx 
        
        j = 1;
        i = 1;
        while j <= np*length(B)
            if j ~= 1
                H(j:j+(length(B)-1),:) = H_x(i:i+(length(B)-1),:);
                i = i + length(B);
            end
            j = j + length(B);
        end
        %End of function  
        end
        function [r_val] = modify_constraint(a_matrix,b_constraint,x_robot,y_robot,np,nx,y_high,int_num)
            if int_num == 0
                y_max = y_high; %the max value for y
                k = -a_matrix(1); %slope of the constraint
                b = b_constraint;
                x_max = (y_max - b)/k;
                x = [x_robot:0.1:x_max];
                y = k*x+b;
                I = find(y >= y_robot + 0.2);
                index_ref = I(1);
                x_ref = x(I(1));
                y_ref = y(I(1));
                eta_ref = atan((y(1,end)-y(1,1))/(x(1,end)-x(1,1)));
                r_val = zeros(np*nx,1);
            end
            j = 1;
            i = 0;
            while j <= np*nx
                r_val(j,1) = x(index_ref + i);
                r_val(j+1,1) = y(index_ref + i);
                r_val(j+2,1) = eta_ref;
                r_val(j+2,1) = 0;
                i = i + 1;
                j = j + nx;
            end
        end
        %obstacle bypass function
        function [x_constraint_high, x_constraint_low, r_val] = bypass_obstacle(x_map,y_map,x_con_high, x_con_low, r_val_old,x_robot,y_robot,robot_front_up,robot_front_down,robot_back_up,robot_back_down,ob_up_left,ob_up_right,ob_down_left,ob_down_right,back_mid,nx,np)
            if robot_back_down(2) < ob_up_left(2) && robot_back_down(1) < ob_up_right(1) + 2 %scenario 1 when the robot is below and behind the obstacle
                x_constraint_low = x_con_low;
                j = 1;
                while j <= (np+1)*nx
                    x_con_high(j,1) = ob_up_left(1);
                    j = j + nx;
                end
                x_constraint_high = x_con_high;
                y_ref = [ob_up_left(2):0.2:x_con_high(2)];
                x_ref = ob_up_left(1)*ones(size(y_ref));
                eta_ref = pi/2;
                I = find(y_ref >= robot_back_down(2) + 0.5);
                index_ref = I(1);
                i = 0;
                j = 1;
                r_val = zeros(np*nx,1);
                while j <= np*nx
                    r_val(j,1) = x_ref(index_ref + i);
                    r_val(j+1,1) = y_ref(index_ref + i);
                    r_val(j+2,1) = eta_ref;
                    r_val(j+3,1) = 0;
                    i = i + 1;
                    j = j + nx;
                end
            elseif robot_back_down(2) >= ob_up_left(2) && robot_back_down(1) < ob_up_right(1) + 2 %scenario 2 when the robot is above and behind the obstacle
                x_constraint_high = x_con_high;
                j = 1;
                while j <= (np+1)*nx
                    x_con_low(j+1,1) = ob_up_left(2)+0.3;
                    j = j + nx;
                end
                x_constraint_low = x_con_low;
                x_ref = [ob_up_right(1):0.2:x_con_high(1)];
                y_ref = ob_up_right(2)*ones(size(x_ref));
                eta_ref = 0;
                I = find(x_ref >= x_robot + 0.5);
                index_ref = I(1);
                i = 0;
                j = 1;
                r_val = zeros(np*nx,1);
                while j <= np*nx
                    r_val(j,1) = x_ref(index_ref + i);
                    r_val(j+1,1) = y_ref(index_ref + i);
                    r_val(j+2,1) = eta_ref;
                    r_val(j+3,1) = 0;
                    i = i + 1;
                    j = j + nx;
                end
            else %scenario 3 when the robot passes to the front of the obstacle
                x_constraint_high = x_con_high;
                j = 1;
                while j <= (np+1)*nx
                    x_con_low(j,1) = robot_back_down(1);
                    j = j + nx;
                end
                x_constraint_low = x_con_low;
                I = find(x_map >= robot_front_down(1) + 0.5);
                index_ref = I(1);
                i = 0;
                j = 1;
                r_val = zeros(np*nx,1);
                while j <= np*nx
                    r_val(j,1) = x_map(index_ref + i);
                    r_val(j+1,1) = y_map(index_ref + i);
                    r_val(j+2,1) = 0;
                    r_val(j+3,1) = 0;
                    i = i + 1;
                    j = j + nx;
                end
            end
        end
        
        function [vx_1,vy_1] = create_rectangle(x,y,theta,s)
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
            vx_1(1,1) = x_d_back - len_perp*cos(angle_val);
            vy_1(1,1) = y_d_back + len_perp*sin(angle_val);
            vx_1(1,2) = x_d_back + len_perp*cos(angle_val);
            vy_1(1,2) = y_d_back - len_perp*sin(angle_val);
            vx_1(1,3) = x_d_front + len_perp*sin(theta);
            vy_1(1,3) = y_d_front - len_perp*cos(theta);
            vx_1(1,4) = x_d_front - len_perp*cos(angle_val);
            vy_1(1,4) = y_d_front + len_perp*sin(angle_val);
        end
    end
end