classdef BicycleModel_v2 < handle
    % state is 
    % X = [x_cm,y_cm,z_cm,  % position in ecef frame (m)
    %       theta,          % heading in E-N plane (rad)
    %       v_xy,v_z,        % velocity in enu E-N plane and vertical (m/s)
    %       b_clk, b_dot]'         % gps receiver clock bias and rate(m)
    % control input is 
    % u = [psi,             % steering angle, assumed to be normally distributed zero mean with sigma proportional to 1/v
    %       F]'           % input longitudinal force 
    
    % model assumptions:
    %   no slip dynamics
    %   negligible grade
    %   aerodynamic resistance
    properties (Constant)
        l1 = 0.355; % (m) long. dist from CoM to rear wheel contact patch
        l2 = 0.686; % (m) long. dist from CoM to front wheel contact patch
        l = BicycleModel_v2.l1+BicycleModel_v2.l2;
        m = 70; % (kg) mass 
        r_w = 0.683/2; % (m) radius of rear wheel
        CdA = 0*0.34; % average from  https://link.springer.com/article/10.1007/s12283-017-0234-1
        rho = 1.225; % kg/m^3 density of air
    end
    methods(Static)

        function X_k1 = f_dyn(X_k,u_k,dt)
            % linearized dynamics at given state 
            X_k1 = X_k + dt*BicycleModel_v2.f_dyn_nonlinear(X_k,u_k);
            % wrap theta
            %X_k1(4) = wrapToPi(X_k1(4));
        end

        function X_dot = f_dyn_nonlinear(X,u)
            % full nonlinear, continuous time dynamics model
            x_ecef = X(1);y_ecef = X(2);z_ecef = X(3);theta = X(4);
            v_en = X(5);v_u = X(6);b_clk = X(7);b_dot = X(8);
            psi = u(1); F = u(2);

            %%%%%%%%%%%%
            v_u = 0;

            l = BicycleModel_v2.l;
            r_w = BicycleModel_v2.r_w;
            m = BicycleModel_v2.m;
            rho = BicycleModel_v2.rho;
            CdA = BicycleModel_v2.CdA;
            R_enu = BicycleModel_v2.R_enu([x_ecef;y_ecef;z_ecef]);
            % define a vector of zeros for cleanliness
            z_vec = @(n) zeros(n,1);

            V_enu = [v_en*cos(theta);v_en*sin(theta);v_u];

            X_ecef_dot = R_enu' * V_enu;

            theta_dot = v_en*tan(psi)/l;

            ven_dot = 1/m*(F-0.5*CdA*rho*v_en^2*sign(v_en));
            b_clk_dot = b_dot;
            b_clk_dot = 0;
            
            vu_dot = 0; bdot_dot = 0;

            X_dot = [X_ecef_dot', theta_dot,...
                     ven_dot, vu_dot, b_clk_dot, bdot_dot]';
        end
        

        function jacobian = F_jac(X,u,dt)
            % linearized dynamics
            x_ecef = X(1);y_ecef = X(2);z_ecef = X(3);theta = X(4);
            v_en = X(5);v_u = X(6);b_clk = X(7);b_dot = X(8);
            %%%%%%%%%
            v_u = 0;

            psi = u(1); F = u(2);
            l = BicycleModel_v2.l;
            r_w = BicycleModel_v2.r_w;
            m = BicycleModel_v2.m;
            rho = BicycleModel_v2.rho;
            CdA = BicycleModel_v2.CdA;
            R_enu = BicycleModel_v2.R_enu([x_ecef;y_ecef;z_ecef]);

            % define a vector of zeros for cleanliness
            z_vec = @(n) zeros(n,1);
            % coefficient term in velocity derivative

            % partial derivatives wrt individual states
            p_px = z_vec(8);
            p_py = z_vec(8);
            p_pz = z_vec(8);
            p_pth = [ R_enu'*[-v_en*sin(theta);v_en*cos(theta); v_u]; z_vec(5)];
            p_pven = [ R_enu'*[cos(theta);sin(theta);0]; tan(psi)/l; CdA*rho/m*v_en*sign(v_en); z_vec(3)];
            p_pvu = [R_enu'*[0;0;1]; z_vec(5)];
            %%%%%%%%%%%%%%%
            p_pvu = [R_enu'*[0;0;0]; z_vec(5)];
            p_pb = z_vec(8);
            p_pbdot = [z_vec(6);1;0];
            jacobian = eye(8) + dt*[p_px,p_py,p_pz,p_pth,p_pven,p_pvu,p_pb,p_pbdot];
    %                       % dividing by the velocity squared makes the partial derivatives very clean 
    %                       % dividing by the velocity squared makes the partial derivatives very clean 
        end
        
        function B = B_jac(X,u,dt)
            % partial derivative of dynamics w.r.t control inputs 
            % used to model control input noise as process noise
            x_ecef = X(1);y_ecef = X(2);z_ecef = X(3);theta = X(4);
            v_en = X(5);v_u = X(6);b_clk = X(7);b_dot = X(8);
            psi = u(1); F = u(2);
            l = BicycleModel_v2.l;
            r_w = BicycleModel_v2.r_w;
            m = BicycleModel_v2.m;
            rho = BicycleModel_v2.rho;
            CdA = BicycleModel_v2.CdA;
            R_enu = BicycleModel_v2.R_enu([x_ecef;y_ecef;z_ecef]);
            % define a vector of zeros for cleanliness
            z_vec = @(n) zeros(n,1);

            pth_ppsi = v_en*sec(psi)^2/l;
            pvxy_pF = 1/m;
            pvxy_ppsi = 0; % i am a dumb

            pXdot_ppsi = [z_vec(3);pth_ppsi;pvxy_ppsi;z_vec(3)];
            pXdot_pF = [z_vec(4);pvxy_pF;z_vec(3)];

            B = dt*[pXdot_ppsi,pXdot_pF];
        end

        function Q = Q_dyn(X,u,dt)
            % generates the process and control input noise covariance matrix
            v_en = X(5);

            Q_process = diag([zeros(1,3),... % assigned later
                    .4^2,... % theta
                    .5,... % ven
                    .5,... % vu
                    1, ... % clock bias
                    5, ... % clock bias rate
                    ])*dt;
            Q_process(1:3,1:3) = BicycleModel_v2.R_enu(X(1:3))'*diag([1.5^2,1.5^2,1^2]*dt) * R_enu(X(1:3));


            Q = Q_process;%+ B*diag([s_sq_psi,s_sq_F])*B';
        end
    end

    methods(Static)
        function R = R_enu(p_ecef)
            % gets the rotation matrix from ecef to enu at a point p_ecef in meters
            geod = Toolbox.geocTogeod(Toolbox.ECEFTogeoc(p_ecef/1000));
            R = Toolbox.ECEFToENU(eye(3),geod(2),geod(3));
        end
    end
end
