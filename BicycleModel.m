classdef BicycleModel < handle
    % state is 
    % X = [x_cm,y_cm,z_cm,  % position in ecef frame (m)
    %       theta,          % heading in plane (rad)
    %       v_x,v_y,        % velocity in plane (m/s)
    %       omega,          % wheel angular rate (rad/s)
    %       b_clk]'         % gps receiver clock bias (m)
    % control input is 
    % u = [psi,             % steering angle, assumed to be normally distributed zero mean with sigma proportional to 1/v
    %       tau]'           % input torque, assumed to be normally distributed zero mean with constant sigma
    % process noise in z_cm,b_clk], normally distributed i.i.d 
    % w = [0,0,w_z,
    %       0,0,0,0,w_b]';
    %
    % model assumptions:
    %   longitudinal slip in driven tire
    %   no side slip
    %   negligible grade
    %   rear wheel normal force is constant
    properties (Constant)
        l1 = 0.355; % (m) long. dist from CoM to rear wheel contact patch
        l2 = 0.686; % (m) long. dist from CoM to front wheel contact patch
        m = 70; % (kg) mass 
        r_w = 0.683/2; % (m) radius of rear wheel
        I_w = 1.470/BicycleModel.r_w^2; % moment of inertia of wheel 
        slip_curve = [0,.0625,.125,.1875,.25,.375,.5,.625,.75,.875,1]; % slip parameters of mu-slip curve
        mu_curve = [0,800,900,925,890,875,812,800,750,730,650]/1000; % mu parameters of mu-slip curve
    end
    methods(Static)

        function jacobian = F_jac(X,u,dt)
            % linearized dynamics
            x_cm = X(1);y_cm = X(2);z_cm = X(3);theta = X(4);
            v_x = X(5);v_y = X(6);omega = X(7);b_clk = X(8);
            psi = u(1); tau = u(2);
            l1 = BicycleModel.l1;
            l2 = BicycleModel.l2;
            r_w = BicycleModel.r_w;
            m = BicycleModel.m;
            I_w = BicycleModel.I_w;
            v = norm([v_x;v_y]);
            s = BicycleModel.calc_slip(v,omega);
            mu = BicycleModel.mu_slip(s);
            g = 9.81;
            v_dot = mu*g*l1/(l1+l2);
            v_dot = sign(v_dot)*min(abs(v_dot),5);
            
            X_dot = BicycleModel.f_dyn_nonlinear(X,u);
            theta_dot = X_dot(4);
            if abs(v)<1e-2
                pv_pvx = 0;
                pv_pvy = 0;
                ps_pom = 0;
                ps_pv = 0;
            else
                pv_pvx = v_x/v;
                pv_pvy = v_y/v;
                ps_pom = r_w/v;
                ps_pv = -r_w/v^2;
            end
            pmu_ps = BicycleModel.mu_slip_ds(s);
            pvdot_pmu = g*l1/(l1+l2);


            %X = [x,y,z,th,vx,vy,om,b]'
            % partial derivatives of f_dyn w.r.t individual states
            p_px = zeros(8,1);
            p_py = zeros(8,1);
            p_pz = zeros(8,1);

            p_pth = [0,0,0,0, -v_dot*sin(theta)-v*cos(theta)*theta_dot,...
                                v_dot*cos(theta)-v*sin(theta)*theta_dot,...
                                0,0]';
            p_pvx = [1,0,0, tan(psi)/(l1+l2)*pv_pvx,...
                            cos(theta)*pvdot_pmu*pmu_ps*ps_pv*pv_pvx - sin(theta)*tan(psi)/(l1+l2)*2*v*pv_pvx,...
                            sin(theta)*pvdot_pmu*pmu_ps*ps_pv*pv_pvx + cos(theta)*tan(psi)/(l1+l2)*2*v*pv_pvx,...
                            1/I_w*m*r_w* pvdot_pmu*pmu_ps*ps_pv*pv_pvx,...
                            0]';
            p_pvy = [0,1,0, tan(psi)/(l1+l2)*pv_pvy,...
                            cos(theta)*pvdot_pmu*pmu_ps*ps_pv*pv_pvy - sin(theta)*tan(psi)/(l1+l2)*2*v*pv_pvy,...
                            sin(theta)*pvdot_pmu*pmu_ps*ps_pv*pv_pvy + cos(theta)*tan(psi)/(l1+l2)*2*v*pv_pvy,...
                            1/I_w*m*r_w* pvdot_pmu*pmu_ps*ps_pv*pv_pvy,...
                            0]';
            p_om = [0,0,0,0,...
                            cos(theta)*pvdot_pmu*pmu_ps*ps_pom,...
                            sin(theta)*pvdot_pmu*pmu_ps*ps_pom,...
                            1/I_w*m*r_w* pvdot_pmu*pmu_ps*ps_pom,...
                            0]';
            p_b = zeros(8,1);

            jacobian = eye(8) + dt*[p_px,p_py,p_pz,p_pth,p_pvx,p_pvy,p_om,p_b];
        end
        
        function B = B_jac(X,u,dt)
            % partial derivative of dynamics w.r.t control inputs 
            % used to model control input noise as process noise
            x_cm = X(1);y_cm = X(2);z_cm = X(3);theta = X(4);
            v_x = X(5);v_y = X(6);omega = X(7);b_clk = X(8);
            psi = u(1); tau = u(2);
            l1 = BicycleModel.l1;
            l2 = BicycleModel.l2;
            r_w = BicycleModel.r_w;
            m = BicycleModel.m;
            I_w = BicycleModel.I_w;
            v = norm([v_x;v_y]);
            
            ptheta_ppsi = v*sec(psi)^2/(l1+l2);
            pxdot_ppsi = [0,0,0,ptheta_ppsi,...
                    -v*sin(theta)*ptheta_ppsi,...
                    v*cos(theta)*ptheta_ppsi,...
                    0,0]';
            pxdot_ptau = [zeros(6,1);1/I_w;0];
            B = dt*[pxdot_ppsi,pxdot_ptau];
        end

        function X_k1 = f_dyn(X_k,u_k,dt)
            % linearized dynamics at given state 
            X_k1 = X_k + dt*BicycleModel.f_dyn_nonlinear(X_k,u_k);
            if X_k1(7)<1
                X_k1(7) = 0;
            end
        end

        function X_dot = f_dyn_nonlinear(X,u)
            % full nonlinear, continuous time dynamics model
            x_cm = X(1);y_cm = X(2);z_cm = X(3);theta = X(4);
            v_x = X(5);v_y = X(6);omega = X(7);b_clk = X(8);
            psi = u(1); tau = u(2);
            l1 = BicycleModel.l1;
            l2 = BicycleModel.l2;
            r_w = BicycleModel.r_w;
            m = BicycleModel.m;
            I_w = BicycleModel.I_w;

            x_dot = v_x;
            y_dot = v_y;

            v = norm([v_x;v_y]);

            s = BicycleModel.calc_slip(v,omega);
            mu = BicycleModel.mu_slip(s);
            g = 9.81;
            v_dot = .5*mu*g*l1/(l1+l2);
            v_dot = sign(v_dot)*min(abs(v_dot),5);

            theta_dot = v*tan(psi)/(l1+l2);
            v_x_dot = v_dot*cos(theta) - v*sin(theta)*theta_dot;
            v_y_dot = v_dot*sin(theta) + v*cos(theta)*theta_dot;
            omega_dot = .5*1/I_w*( tau-v_dot*m*r_w );

            % unmodeled dynamics
            z_dot = 0; b_clk_dot = 0;
            
            X_dot = [x_dot, y_dot, z_dot, theta_dot,...
                     v_x_dot, v_y_dot, omega_dot, b_clk_dot]';
        end
        
        function Q = Q_dyn(X,u,dt)
            % generates the process and control input noise covariance matrix
            x_cm = X(1);y_cm = X(2);z_cm = X(3);theta = X(4);
            v_x = X(5);v_y = X(6);omega = X(7);b_clk = X(8);
            psi = u(1); tau = u(2);
            v = norm([v_x;v_y]);
        

            % TODO: need to tune
            Sigma_process = diag([.4*ones(1,2),.4,... % xyz
                    .75/4,... % theta
                    .6/2,.6/2,... % vxy
                    .7/3,... % omega
                    1.3, ... % clock bias
                    ])*dt;
            B = BicycleModel.B_jac(X,u,dt);
            % assume steer angle standard deviation is 
            % inversely proportional to speed
            sigma_psi = (1.25/(v+.01))^2;
            %sigma_psi = (0)^2;
            sigma_tau = (.75)^2;
            Q = B*diag([sigma_psi,sigma_tau])*B' + Sigma_process;
        end
    end
    methods(Static)

        function s = calc_slip(v,om)
            if abs(v)<1e-2
                s = 0.0000;
            else
                s = (BicycleModel.r_w*om-v)./v;
            end
        end

        function mu = mu_slip(s)
            % mu-slip curve described by a piecewise function
            % adapted from data in https://www.researchgate.net/profile/Oliver-Maier-6/publication/308612161_Vertical_and_Longitudinal_Characteristics_of_a_Bicycle_Tire/links/57e8dc0808ae9e5e4558d590/Vertical-and-Longitudinal-Characteristics-of-a-Bicycle-Tire.pdf
            % using the 6 bar inflation pressure, 
            % Fz = 1.0 kN chart
            mu = interp1(BicycleModel.slip_curve,...
                BicycleModel.mu_curve,...
                abs(s),'linear','extrap');
            mu = mu.*sign(s);
        end

        function dmu = mu_slip_ds(s)
            % takes the piecewise derivative of the mu-slip curve at s
            deriv = diff(BicycleModel.mu_curve)./diff(BicycleModel.slip_curve);
            dmu = interp1(BicycleModel.slip_curve(1:end-1),...
                    deriv,abs(s),'previous','extrap');
        end
    end
end
