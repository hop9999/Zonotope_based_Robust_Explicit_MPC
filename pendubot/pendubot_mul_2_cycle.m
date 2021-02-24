clear all;
clc;
close all;
cvx_clear;

n = 4; m = 1;
z = 250;

dt = 10e-3;
%number of iteration
Count = 99;

%%%%%%%%%%%%%%%%%%%
%model params
traj = readmatrix('traj_cycle_fast.csv');
% 
% t1 = linspace(0,3,300);
% t2 = linspace(0,3,Count);
% % % 
% traj = interp1(t1,traj',t2)';

x_des = traj(1:4,:);
u_des = traj(5,:);

A1 = zeros(n,n,Count);
B1 = zeros(n,m,Count);

A2 = zeros(n,n,Count);
B2 = zeros(n,m,Count);

A3 = zeros(n,n,Count);
B3 = zeros(n,m,Count);

A4 = zeros(n,n,Count);
B4 = zeros(n,m,Count);

p = [0.0329
    0.0060
    0.0038
    0.1265
    0.0197
    0.0038
    0.0004
    0.0466
    0.0000];

me1 = 0.04;
le1 = 0.2;
me2 = 0.01;
le2 = 0.25;

for i = 1:Count
    x_i = x_des(:,i);
    u_i = u_des(:,i);

    A1(:,:,i) = A_mtrx(x_i(1:2), x_i(3:4),u_i,p,-0.01,le1,-me2,le2)*dt + eye(n);
    B1(:,:,i) = B_mtrx(x_i(1:2), x_i(3:4),u_i,p,-0.01,le1,-me2,le2)*dt;
    
    A2(:,:,i) = A_mtrx(x_i(1:2), x_i(3:4),u_i,p,me1,le1,-me2,le2)*dt + eye(n);
    B2(:,:,i) = B_mtrx(x_i(1:2), x_i(3:4),u_i,p,me1,le1,-me2,le2)*dt;
    
    A3(:,:,i) = A_mtrx(x_i(1:2), x_i(3:4),u_i,p,-0.01,le1,me2,le2)*dt + eye(n);
    B3(:,:,i) = B_mtrx(x_i(1:2), x_i(3:4),u_i,p,-0.01,le1,me2,le2)*dt;

    A4(:,:,i) = A_mtrx(x_i(1:2), x_i(3:4),u_i,p,me1,le1,me2,le2)*dt + eye(n);
    B4(:,:,i) = B_mtrx(x_i(1:2), x_i(3:4),u_i,p,me1,le1,me2,le2)*dt;


end

x_0 = [0.0;0.0;0.0;0.0];
x_g = [0.0;0.0;0.0;0.0];

G_start = diag([0.2 0.2 2.0 2.0]);

G_end = diag([0.2 0.2 1.2 1.2]);

u_lim = 10;
T_max = [u_lim];

%%%%%%%%%%%%%%%%%%%
%disturbance
%dist zonotop size
w = 4;
W = 0.95*diag([0.001 0.001 0.001 0.001]);
%%%%%%%%%%%%%%%%%%%
y = 1;

cvx_solver gurobi_2
cvx_begin
    cvx_solver_settings( 'BarHomogeneous', 1 )
    cvx_solver_settings( 'CrossoverBasis', 1 )
    cvx_solver_settings( 'NumericFocus', 3 )
    cvx_solver_settings( 'Method', 3 )
    cvx_precision best 
    variable G(n, z+n, (Count+1))
    variable T(m, z+n, Count)
    
   % variable x(n, (Count+1))
   % variable u(m, (Count+1))
    
   % variable Gamma_end(n, z+n)
   % variable betta_end(z+n, 1)
    
   %variable Gamma_u(m, z+n, Count)
   %variable betta_u(m, 1, Count)
    
    variable bounding_box_margin(n,Count) nonnegative
 
    Tvectorized = reshape(T(:, :, :), m*(z+n )*Count, []);
    Gvectorized_x1 = reshape(G(1, :, 2:end),  (z+n)*Count, []);
    Gvectorized_dx1 = reshape(G(3, :, 2:end), (z+n)*Count, []);
    Gvectorized_x2 = reshape(G(2, :, 2:end),  (z+n)*Count, []);
    Gvectorized_dx2 = reshape(G(4, :, 2:end), (z+n)*Count, []);
    bounding_box_margin_vect = reshape(bounding_box_margin,n*Count,[]);
    
    cost_norm = 2;
    
    minimize(1*norm(bounding_box_margin_vect,cost_norm) + 1*norm(Tvectorized,cost_norm) + 100*norm(Gvectorized_x1,cost_norm) + 100*norm(Gvectorized_dx1,cost_norm) + 100*norm(Gvectorized_x2,cost_norm) + 100*norm(Gvectorized_dx2,cost_norm)  )

    subject to

        %init
        G(:, :, 1) == G(:, :, Count+1)
        %goal
        %G(:, :, Count+1) == G_end*Gamma_end
        %x_g - x(:,Count+1) == G_start*betta_end
        %norm(Gamma_end,inf) <= 1
        
        for i = 1:Count

            %control limit
            %T(:, :, i) == T_max*Gamma_u(:,:,i)
            %0 - u(:,i) == T_max*betta_u(:,:,i)
            %norm([Gamma_u(:,:,i)],Inf) <= 1

            G1 = A1(:,:,i)*G(:, :, i) + B1(:,:,i)*T(:,:,i);
            G2 = A2(:,:,i)*G(:, :, i) + B2(:,:,i)*T(:,:,i);
            G3 = A3(:,:,i)*G(:, :, i) + B3(:,:,i)*T(:,:,i);
            G4 = A4(:,:,i)*G(:, :, i) + B4(:,:,i)*T(:,:,i);
            
            F1 = [G1+G2 G1-G2]/2;
            F2 = [G3+G4 G3-G4]/2;
            
            F = [F1+F2 F1-F2]/2;
            
            F = [F(:,1:4) F(:,9:end) F(:,5:8)];
            F_1 = F(:,1:z-n);
            F_2 = F(:,z-n+1:end);
            for j = 1:n
                norm(F_2(j,:),1) <= bounding_box_margin(j,i)
            end
            G(:, :, i+1) == [F_1 diag(bounding_box_margin(:,i)) W]

        end
cvx_end


if strcmp(cvx_status, 'Inaccurate/Solved') || strcmp(cvx_status, 'Solved')
    Z_red = zeros(n,n,Count);
    Z1 = zeros(n,n,Count);
    Z2 = zeros(1,n,Count);
    for i = 1:Count
        c = x_des(:,i)';
        X = [G(:,:, i) -G(:,:, i)]';
        C0 = X'*X;
        [U,S,V] = svd(C0);
        Z_ad = V*G(:,:, i);
        Z_red(:,:,i) = V'*diag([norm(Z_ad(1,:),1) norm(Z_ad(2,:),1) norm(Z_ad(3,:),1) norm(Z_ad(4,:),1)]);
        Z1(:,:,i) = Z_red(:,:,i)*diag([1/norm(Z_red(:,1,i),2) 1/norm(Z_red(:,2,i),2)  1/norm(Z_red(:,3,i),2)  1/norm(Z_red(:,4,i),2)]);
        Z2(:,:,i) = [norm(Z_red(:,1,i),2) norm(Z_red(:,2,i),2) norm(Z_red(:,3,i),2) norm(Z_red(:,3,i),2)];
        draw_zonotope([Z_red(1,:,i);Z_red(3,:,i)], [x_des(1,i);x_des(3,i)]);
        draw_zonotope([Z_red(2,:,i);Z_red(4,:,i)], [x_des(2,i);x_des(4,i)]);
    end
    %plot(x_des(1,:),x_des(3,:),'*')
    %plot(x_des(2,:),x_des(4,:),'*')
    figure
    T;
    
    %%%%%%%%%%%%%%%%%
    
    Point_Count = 150;
    
    Points_gen = 2*rand(z+n, Point_Count) - ones(z+n, Point_Count); %for {-1, 1} convention
    
    Xset = zeros(n, Point_Count, Count);
    Uset = zeros(m, Point_Count, Count);
    Policy = zeros(m,n,Count);
    % Xset = cell(Count+1, 1);
    % Uset = cell(Count, 1);
    
    Xset(:, :, 1) = G(:, :, 1)*Points_gen;
    for i = 1:Count
        T(:, :, i) * pinv(G(:,:, i))
        Policy(:,:,i) = T(:, :, i) * pinv(G(:,:, i));
        %Uset(:, :, i) = u(i) + T(:, :, i) * pinv(G(:, :, i)) * (Xset(:, :, i) - x(:,i));
        Uset(:, :, i) = T(:, :, i) * pinv(G(:, :, i)) * (Xset(:, :, i));
        Dist_gen = 2*rand(w, Point_Count) - ones(w, Point_Count); %for {-1, 1} convention
        Wset = W*Dist_gen;
        Xset(:, :, i+1) = (A1(:, :, i) + A2(:, :, i))/2*Xset(:, :, i) + (B1(:, :, i) + B2(:, :, i))/2*Uset(:, :, i) + Wset;
        %         if mod(i,2)==1
        %             draw_zonotope(G_max(:, :, 1), x(:,i),'FaceColor', [0.1, 0.4, 0.1]);
        %         end
        
    end
    
    %%%%%%%%%%%%
    for i = 1:Point_Count
        plot(reshape(Xset(1, i, :), [], 1), ...
            reshape(Xset(3, i, :), [], 1), ...
            'LineWidth', 3); hold on;
    end
    
    for i = 1:Point_Count
        plot(reshape(Xset(1, i, 1), [], 1), ...
            reshape(Xset(3, i, 1), [], 1), ...
            'o', 'MarkerFaceColor', 'b'); hold on;
    end
    
    for i = 1:Point_Count
        plot(reshape(Xset(1, i, end), [], 1), ...
            reshape(Xset(3, i, end), [], 1), ...
            '^', 'MarkerFaceColor', 'g', 'MarkerSize', 10); hold on;
    end
    
    %%%%%%%%%%%%%%
    figure()
    for i = 1:Point_Count
        plot(reshape(Xset(2, i, :), [], 1), ...
            reshape(Xset(4, i, :), [], 1), ...
            'LineWidth', 3); hold on;
    end
    
    for i = 1:Point_Count
        plot(reshape(Xset(2, i, 1), [], 1), ...
            reshape(Xset(4, i, 1), [], 1), ...
            'o', 'MarkerFaceColor', 'b'); hold on;
    end
    
    for i = 1:Point_Count
        plot(reshape(Xset(2, i, end), [], 1), ...
            reshape(Xset(4, i, end), [], 1), ...
            '^', 'MarkerFaceColor', 'g', 'MarkerSize', 10); hold on;
    end
    
    policy = reshape(Policy,[4,Count])';
    csvwrite('G_cycle_opt.csv',G)
    csvwrite('T_cycle_opt.csv',T)
    t1 = linspace(0,3,Count);
	t2 = linspace(0,3,300);
	new_pol = interp1(t1,policy,t2)';
	csvwrite('pol_cycle_opt.csv',policy);
    Z_rede = reshape(Z_red,n,Count*n);
    Z_1e = reshape(Z1,n,Count*n);
    Z_2e = reshape(Z2,1,Count*n);
    xe = reshape(x_des,1,Count*n);
    Z_data = [Z_rede;Z_1e;Z_2e;xe];
    csvwrite('Z_cycle_opt.csv',Z_data);
    
end
