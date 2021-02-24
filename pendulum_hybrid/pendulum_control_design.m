clear all;
clc;
close all;
cvx_clear;

n = 2; m = 1;
z = 16;

dt = 4e-3;
%number of iteration
Count = 90;

%%%%%%%%%%%%%%%%%%%
%model params
k_i = 0.0306;
mgl = 0.1852 + 0.46*9.8*0.15;
I = 0.0045 + 0.46*0.15^2 + 0.46*0.055^2/2;
k_qf = 0.00109;
k_qc = 0.46;
K = 129;
q_c = 3.254-pi;
q_lin = pi;

%A nominal
Af1 = [0                1
       -1*0.8*mgl*cos(q_lin)/(0.7*I)  -1*k_qf/(0.7*I)]*dt + eye(2);
    
Af2 = [0                1
       -1*0.8*mgl*cos(q_lin)/(1.6*I)  -1*k_qf/(1.6*I)]*dt + eye(2);
   
Af3 = [0                1
       -1*1.2*mgl*cos(q_lin)/(0.7*I)  -1*k_qf/(0.7*I)]*dt + eye(2);
    
Af4 = [0                1
       -1*1.2*mgl*cos(q_lin)/(1.6*I)  -1*k_qf/(1.6*I)]*dt + eye(2);

   
Ac1 = [0                1
       (-mgl*cos(q_lin) - 0.8*K)/I  -0.8*k_qc/I]*dt + eye(2);
  
Ac2 = [0                1
       (-mgl*cos(q_lin) - 0.8*K)/I  -1.2*k_qc/I]*dt + eye(2);
   
Ac3 = [0                1
       (-mgl*cos(q_lin) - 1.2*K)/I  -0.8*k_qc/I]*dt + eye(2);
  
Ac4 = [0                1
       (-mgl*cos(q_lin) - 1.2*K)/I  -1.2*k_qc/I]*dt + eye(2);
   
%%%%%%%%%%%%%%%%%%%%
%B nominal
B = [0
     k_i/I]*dt;
%%%%%%%%%%%%%%%%%%%%
%C nominal
Cf1 = [0
       -0.0*I/I]*dt;

Cf2 = [0
       -0.0*I/I]*dt;
  
Cc1 = [0
      0.8*K*q_c/I]*dt;
  
Cc2 = [0
      1.2*K*q_c/I]*dt;
%%%%%%%%%%%%%%%%%%%


x_0 = [0.0;0.0];
x_g = [0.0;0.0];

G_start = [diag([0.02 0.2]) zeros(n,z)];

G_end = [diag([0.02 0.2]) zeros(n,z)];

u_lim = 8;
T_max = [u_lim];

%limits
x_c = [q_c+0.2;0];
G_c = diag([0.21 10]);

x_f = [q_c-0.3;0];
G_f = diag([0.31 10]);

%%%%%%%%%%%%%%%%%%%
%disturbance
%dist zonotop size
w = 2;
W = diag([0.0001,0.001]);
%%%%%%%%%%%%%%%%%%%
y = 1;

cvx_solver Gurobi_2
cvx_begin
    cvx_precision best
    variable G(n, z+n, (Count+1))
    variable T(m, z+n, Count)
    
    variable x(n, (Count+1))
    variable u(m, (Count+1))
    
    variable Gamma_end(z+n, z+w)
    variable betta_end(z+n, 1)
    
    variable Gamma_f(n, z+w, Count)
    variable betta_f(n, 1, Count)
    
    variable Gamma_c(n, z+w, Count)
    variable betta_c(n, 1, Count)
    
    variable Gamma_u(m, z+n, Count)
    variable betta_u(m, z+n, Count)
    
    variable af(n,Count) nonnegative
    variable ac(n,Count) nonnegative
    variable y nonnegative
    variable z1(Count) binary
    variable du(m, Count);
    M = 100;

    Tvectorized = reshape(T(:, :, :), m*(z+n )*Count, []);
    Gvectorized_x = reshape(G(1, :, 2:end),  (z+n)*Count, []);
    Gvectorized_dx = reshape(G(2, :, 2:end), (z+n)*Count, []);

    af_vect = reshape(af,n*Count,[]);
    ac_vect = reshape(ac,n*Count,[]);
    cost_norm = 1;
    n_mid = floor(Count/2);

    minimize(1*norm(af_vect,1) + 1*norm(ac_vect,1) + 10*norm(u,cost_norm) + 0*norm(x(1,:),cost_norm) + 0*norm(x(2,:),inf) + 1*norm(Tvectorized,cost_norm) + 100*norm(Gvectorized_x,cost_norm) + 100*norm(Gvectorized_dx,cost_norm) + 0*(x(1,n_mid) - (q_c + 0.01))^2 - 00*y)
%наложить доп конмстреинт на зонотопы
    subject to
        
        
        y == 16
        %init
%         G(:, :, 1) == G_start
        x_0 == x(:,1)
        G(:, :, 1) == G(:, :, Count+1)
        x(:,Count+1) == x(:,1)
        %x(:,1) == [0;0]
         %goal
        G(:, :, Count+1) == G_end*Gamma_end
        x_g - x(:,Count+1) == G_end*betta_end
        norm([Gamma_end,betta_end],inf) <= 1
        
        for i = 1:Count
            du(:,i) == u(:,i+1) - u(:,i);
            abs(du(:,i)) <= 0.1
              (x(1,i) - q_c) <= M*z1(i)
              (x(1,i) - q_c) >= -M*(1-z1(i))
            x(1,n_mid) >= q_c + 0.005
            %control limit
            T(:, :, i) == T_max*Gamma_u(:,:,i)
            0 - u(:,i) == T_max*betta_u(:,:,i)
            norm([Gamma_u(:,:,i),betta_u(:,:,i)],Inf) <= 1
            
            xf1 = Af1*x(:,i) + B*u(:,i) + Cf1;
            xf2 = Af2*x(:,i) + B*u(:,i) + Cf1;
            xf3 = Af3*x(:,i) + B*u(:,i) + Cf2;
            xf4 = Af4*x(:,i) + B*u(:,i) + Cf2;
            
            Gf1 = Af1*G(:, :, i) + B*T(:,:,i);
            Gf2 = Af2*G(:, :, i) + B*T(:,:,i);
            Gf3 = Af3*G(:, :, i) + B*T(:,:,i);
            Gf4 = Af4*G(:, :, i) + B*T(:,:,i);
            
            Ff1 = [Gf1 + Gf2 (xf1 - xf2) Gf1 - Gf2]/2;
            Ff2 = [Gf3 + Gf4 (xf3 - xf4) Gf3 - Gf4]/2;
            cf1 = (xf1 + xf2)/2;
            cf2 = (xf3 + xf4)/2;
            
            Ff = [Ff1 + Ff2 (cf1 - cf2) Ff1 - Ff2]/2;
            cf = (cf1 + cf2)/2;
            
            Ff = [Ff(:,1:2) Ff(:,4:end) Ff(:,3)];
            Ff_1 = Ff(:,1:z-n);
            Ff_2 = Ff(:,z-n+1:end);
            for j = 1:n
                norm(Ff_2(j,:),1) <= af(j,i)
            end
            
            xc1 = Ac1*x(:,i) + B*u(:,i) + Cc1;
            xc2 = Ac2*x(:,i) + B*u(:,i) + Cc1;
            xc3 = Ac3*x(:,i) + B*u(:,i) + Cc2;
            xc4 = Ac4*x(:,i) + B*u(:,i) + Cc2;
            
            Gc1 = Ac1*G(:, :, i) + B*T(:,:,i);
            Gc2 = Ac2*G(:, :, i) + B*T(:,:,i);
            Gc3 = Ac3*G(:, :, i) + B*T(:,:,i);
            Gc4 = Ac4*G(:, :, i) + B*T(:,:,i);
            
            Fc1 = [Gc1 + Gc2 (xc1 - xc2) Gc1 - Gc2]/2;
            Fc2 = [Gc3 + Gc4 (xc3 - xc4) Gc3 - Gc4]/2;
            cc1 = (xc1 + xc2)/2;
            cc2 = (xc3 + xc4)/2;
            
            Fc = [Fc1 + Fc2 (cc1 - cc2) Fc1 - Fc2]/2;
            cc = (cc1 + cc2)/2;
            
            Fc = [Fc(:,1:2) Fc(:,5:end) Fc(:,3:4)];
            Fc_1 = Fc(:,1:z-n);
            Fc_2 = Fc(:,z-n+1:end);
            for j = 1:n
                norm(Fc_2(j,:),1) <= ac(j,i)
            end
            
            %hydrid_constr
            G(:, :, i) == G_f*Gamma_f(:, :, i)
            x_f - x(:,i) == G_f*betta_f(:, :, i)
            norm([Gamma_f(:, :, i),betta_f(:, :, i)],inf) <= 1 + 100*z1(i)*M
            
            G(:, :, i) == G_c*Gamma_c(:, :, i)
            x_c - x(:,i) == G_c*betta_c(:, :, i)
            norm([Gamma_c(:, :, i),betta_c(:, :, i)],inf) <= 1 + 100*(1 - z1(i))*M
            
            norm(reshape(G(:, :, i+1) - [Ff_1 diag(af(:,i)),y*W],n*(z+n),[]),1) <= z1(i)*M
            norm(x(:,i+1) - cf,1) <= z1(i)*M
            
            norm(reshape(G(:, :, i+1) - [Fc_1 diag(ac(:,i)),y*W],n*(z+n),[]),1) <= (1 - z1(i))*M
            norm(x(:,i+1) - cc,1) <= (1 - z1(i))*M

        end
cvx_end
 
if strcmp(cvx_status, 'Inaccurate/Solved') || strcmp(cvx_status, 'Solved') || strcmp(cvx_status, 'Suboptimal')
    T;
    
    %%%%%%%%%%%%%%%%%
    
   Point_Count = 40;
    
    Points_gen = 2*rand(z+n, Point_Count) - ones(z+n, Point_Count); %for {-1, 1} convention

    Xset = zeros(n, Point_Count, Count+1);
    Uset = zeros(m, Point_Count, Count);
    Policy = zeros(Count,5);
    Policy_tikh = zeros(Count,5);
    % Xset = cell(Count+1, 1);
    % Uset = cell(Count, 1);
    
    Xset(:, :, 1) = x_0 + G(:, :, 1)*Points_gen;
    for i = 1:Count
        Policy(i,:) = [x(:,i)',u(i), T(:, :, i) * pinv(G(:, :, i))];
        lambda = 0.0001;
        G_inv = inv(G(:, :, i)'*G(:, :, i) + lambda*eye(z+n))*G(:, :, i)';
        Policy_tikh(i,:) = [x(:,i)',u(i), T(:, :, i) * G_inv];
        Uset(:, :, i) = u(i) + T(:, :, i) * G_inv * (Xset(:, :, i) - x(:,i));
        Dist_gen = 2*rand(2, Point_Count) - ones(2, Point_Count); %for {-1, 1} convention
        Wset = y*W*Dist_gen;
%         Ff = [Af*G(:, :, i) + B*T(:,:,i) + Af*G(:, :, i) + B*T(:,:,i) Af*G(:, :, i) + B*T(:,:,i) - Af*G(:, :, i) - B*T(:,:,i)]/2;
%         Ff1 = Ff(:,1:z-n);
%             
%         Fc = [Ac*G(:, :, i) + B*T(:,:,i) + Ac*G(:, :, i) + B*T(:,:,i) Ac*G(:, :, i) + B*T(:,:,i) - Ac*G(:, :, i) - B*T(:,:,i)]/2;
%         Fc1 = Fc(:,1:z-n);
%             %i
        if z1(i) == 0
            A1 = Af1;
            A2 = Af2;
            A3 = Af3;
            A4 = Af4;
            C1 = Cf1;
            C2 = Cf1;
            C3 = Cf2;
            C4 = Cf2;
            %z1(i)*M 
            %norm(reshape(G(:, :, i+1) - [Ff1 diag(af(:,i)),W],n*(z+n),[]),1) ;
        else
            A1 = Ac1;
            A2 = Ac2;
            A3 = Ac3;
            A4 = Ac4;
            C1 = Cc1;
            C2 = Cc1;
            C3 = Cc2;
            C4 = Cc2;
            
            %(1 - z1(i))*M 
            %norm(reshape(G(:, :, i+1) - [Fc1 diag(ac(:,i)),W],n*(z+n),[]),1) ;
        end

        %disp([z1(i),A,C])
        draw_zonotope(G(:, :, i), x(:,i));
        Xset(:, 1:10, i+1) = A1*Xset(:, 1:10, i) + B*Uset(:, 1:10, i) + C1 + Wset(:, 1:10);
        Xset(:, 11:20, i+1) = A2*Xset(:, 11:20, i) + B*Uset(:, 11:20, i) + C2 + Wset(:, 11:20);
        Xset(:, 21:30, i+1) = A3*Xset(:, 21:30, i) + B*Uset(:, 21:30, i) + C3 + Wset(:, 21:30);
        Xset(:, 31:40, i+1) = A4*Xset(:, 31:40, i) + B*Uset(:, 31:40, i) + C4 + Wset(:, 31:40);
    end
    
%     
%     for i = 1:Point_Count
%         plot(reshape(Xset(1, i, :), [], 1), ...
%             reshape(Xset(2, i, :), [], 1), ...
%             'LineWidth', 3); hold on;
%     end
%     
%     for i = 1:Point_Count
%         plot(reshape(Xset(1, i, 1), [], 1), ...
%             reshape(Xset(2, i, 1), [], 1), ...
%             'o', 'MarkerFaceColor', 'b'); hold on;
%     end
%     
%     for i = 1:Point_Count
%         plot(reshape(Xset(1, i, end), [], 1), ...
%             reshape(Xset(2, i, end), [], 1), ...
%             '^', 'MarkerFaceColor', 'g', 'MarkerSize', 10); hold on;
%     end
    Policy(:,1) = Policy(:,1) + q_lin;
    Policy_tikh(:,1) = Policy_tikh(:,1) + q_lin;
    draw_zonotope(G(:, 1:z, 1), x(:,1));
    %draw_zonotope(G_max(:, :, 1), x(:,1));
    %draw_zonotope(G_end(:, :, 1), x(:,end));
    %plot(x(1,:),x(2,:),'r*')
    plot([q_c q_c],[-2.5 2.5],'k','LineWidth',2)
    %xlim([-0.15 0.15])
    %ylim([-1 1])
    xlabel('x')
    ylabel('dx')
    grid on
    csvwrite('policy_point.csv',Policy_tikh)
    csvwrite('G_pend.csv',G)
    csvwrite('T_pend.csv',T)
    t1 = linspace(0,dt*Count,Count);
	t2 = linspace(0,dt*Count,2*Count);
	new_pol = interp1(t1,Policy,t2)';
	csvwrite('pol_point_500.csv',new_pol');
end
Z_red = zeros(2,2,Count);
Z1 = zeros(2,2,Count);
Z2 = zeros(1,2,Count);
for i = 1:Count
    c = x(:,i)';%[policy(i,1) - pi; policy(i,2)];
    X = [G(:,:,i) -G(:,:,i)]';
    C0 = X'*X;
    [U,S,V] = svd(C0);
    Z_ad = V*G(:,:,i);
    Z_red(:,:,i) = V'*diag([norm(Z_ad(1,:),1) norm(Z_ad(2,:),1)]);
    Z1(:,:,i) = Z_red(:,:,i)*diag([1/norm(Z_red(:,1,i),2) 1/norm(Z_red(:,2,i),2)]);
    Z2(:,:,i) = [norm(Z_red(:,1,i),2) norm(Z_red(:,2,i),2)];
    %draw_zonotope(Z_red(:,:,i), c');
    %draw_zonotope(G, x(i,:)');
end
x_r = [x(1,1:end-1) + pi; x(2,1:end-1)]
Z_data = [reshape(Z_red,2,Count*2);reshape(Z1,2,Count*2);reshape(Z2,1,Count*2);reshape(x_r,1,Count*2)];
csvwrite('Z_data.csv',Z_data);