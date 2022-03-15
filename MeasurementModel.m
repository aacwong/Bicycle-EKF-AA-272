classdef MeasurementModel < handle
    % measurements are input as structs with fields dependent on type
    % measurements are either 
    %
    % type=='wheelspeed':
    % z = [om]              % measured wheel speed of driven wheel
    % meas.z = om;
    %
    % type=='gnss':
    % z = [rho;             % pseudorange
    %       rho_dot]        % pseudorange_rate
    % meas.z = [rho;rho_dot];
    % meas.sigma = [sigma_rho;sigma_rho_dot];
    % meas.sv_pos = [x_sv;y_sv;z_sv];   % position of satellite in ecef 
    % meas.sv_vel = [vx_sv;vy_sv;vz_sv];   % vel of satellite in ecef 
    % meas.sv_B = sv_clock_err; % clock bias
    % meas.sv_B_dot = sv_clock_err_rate;    % clock bias rate
    %
    %
    % state is 
    % X = [x_cm,y_cm,z_cm,  % position in ecef frame (m)
    %       theta,          % heading in plane (rad)
    %       v_x,v_y,        % velocity in plane (m/s)
    %       omega,          % wheel angular rate (rad/s)
    %       b_clk]'         % gps receiver clock bias (m)
    properties (Constant)
        c = 299792458; % speed of light
        motion_model = 2; % set which motion model we are using
    end
    methods(Static)
        function h = h_meas(X,meas)
        % returns predicted measurement of type given state and measurement model params
            if strcmp(meas.type,'wheelspeed')
                if MeasurementModel.motion_model==2
                    if abs(X(5)) < BicycleModel_v2.v_thresh
                        h = 0;
                    else
                        h = abs(X(5))/BicycleModel_v2.r_w ; 
                    end
                    % wheel speed measurements are always positive
                else
                    h = X(7);
                end
            elseif strcmp(meas.type,'gnss')

                sat_pos = meas.sv_pos;
                sat_vel = meas.sv_vel;
                tof = meas.tof;
                % correct time of flight
                for i=1:size(sat_pos,2)
                    sat_pos(:,i) = FlightTimeCorrection(sat_pos(:,i),tof(i))';
                end
                B = meas.sv_B;
                B_dot = meas.sv_B_dot;
                c = MeasurementModel.c;
                pos = X(1:3);
                if MeasurementModel.motion_model==2
                    R_enu = BicycleModel_v2.R_enu(pos);
                    vel = [X(5)*cos(X(4));X(5)*sin(X(4));X(6)];
                    b = X(7);
                    bdot = X(8);
                else
                    R_enu = eye(3);
                    vel = [X(5:6);0];
                    b = X(8);
                    bdot = 0;
                end

                e = sat_pos-pos; e = e./vecnorm(e);

                rho = vecnorm(sat_pos-pos) + c*(-B)+b;
                rho_dot = dot(sat_vel-R_enu'*vel,e) +bdot + -c*B_dot;
                h = [rho,rho_dot]';
            else
                error('unsupported measurement type!')
            end
        end

        function H = H_jac(X,meas)
            % jacobian of measurement w.r.t states
            if strcmp(meas.type,'wheelspeed')
                if MeasurementModel.motion_model==2
                    H = zeros(1,8);
                    %if abs(X(5))>=BicycleModel_v2.v_thresh
                        H(5) = sign(X(5))*1/BicycleModel_v2.r_w;
                    %end
                    % wheel speed is always nonzero, so have to account for sign flip
                else
                    H = zeros(1,8);
                    H(7) = 1;
                end
            elseif strcmp(meas.type,'gnss')
                %h = MeasurementModel.h_meas(X,meas);
                sat_pos = meas.sv_pos;
                tof = meas.tof;
                % correct time of flight
                for i=1:size(sat_pos,2)
                    sat_pos(:,i) = FlightTimeCorrection(sat_pos(:,i),tof(i))';
                end
                sat_vel = meas.sv_vel;
                B = meas.sv_B;
                B_dot = meas.sv_B_dot;
                c = MeasurementModel.c;
                pos = X(1:3);
                if MeasurementModel.motion_model==2
                    theta = X(4);
                    v_en = X(5);
                    v_u = X(6);
                    vel_enu = [v_en*cos(theta);v_en*sin(theta);v_u]; % in enu frame
                    bu = X(7);
                    bdot = X(8);
                    e = sat_pos-pos; e = e/norm(e);
                    R_enu = BicycleModel_v2.R_enu(pos);

                    rho_true = vecnorm(sat_pos-pos);
                    prho_pxyz = [(pos-sat_pos)./rho_true]';
                    prho_pb= ones(size(meas.sv_id))';

                    v_rel = sat_vel-R_enu'*vel_enu;

                    % jacobian of unit vector wrt position
                    prho_dot_pxyz = zeros(size(sat_pos))';
                    for i=1:length(meas.sv_id)
                        pe_pxyz = (sat_pos(:,i)-pos)*(sat_pos(:,i)-pos)'/rho_true(i)^3 - 1/rho_true(i)*eye(3);
                        prho_dot_pxyz(i,:) = (v_rel(:,i)'*pe_pxyz);
                    end
                    prho_dot_pven = (-R_enu*e)'*[cos(theta);sin(theta);0];
                    %if abs(v_en)<BicycleModel_v2.v_thresh
                    %    %prho_dot_pven = 0*prho_dot_pven;
                    %else
                    %    v_e = v_en*cos(theta); v_n = v_en*sin(theta);
                    %    dve_dven = 1/cos(theta);dvn_dven = 1/sin(theta);
                    %    if isinf(dve_dven); dve_dven = 0;end
                    %    if isinf(dvn_dven); dvn_dven = 0;end
                    %    prho_dot_pven = (R_enu*e)' * [dve_dven;dvn_dven;0];
                    %end
                    %prho_dot_pven = (-R_enu*e)'*[cos(theta);sin(theta);0];
                    prho_dot_pvu = (-R_enu*e)' *[0;0;1];
                    prho_dot_pth = (-R_enu*e)' *[-v_en*sin(theta);v_en*sin(theta);0];
                else
                    vel = [X(5:6);0];
                    bu = X(8);
                    bdot = 0;
                    e = sat_pos-pos; e = e/norm(e);

                    rho_true = norm(sat_pos-pos);
                    prho_pxyz = [(pos-sat_pos)./rho_true]';
                    prho_pb= 1;

                    v_rel = sat_vel-vel;
                    % jacobian of unit vector wrt position
                    pe_pxyz = (sat_pos-pos)*(sat_pos-pos)'/rho_true^3 - 1/rho_true*eye(3);

                    % really hope this is good math cuz its cool
                    prho_dot_pxyz = (v_rel'*pe_pxyz);
                    prho_dot_pvxyz = -eye(3)* e;
                    prho_dot_pvxy = prho_dot_pvxyz(1:2);
                    prho_dot_pvz = prho_dot_pvxyz(3);

                end
                if MeasurementModel.motion_model==2
                    n = numel(meas.sv_id);
                    H = [prho_pxyz, zeros(n,3), prho_pb, zeros(n,1);
                        prho_dot_pxyz, prho_dot_pth, prho_dot_pven, prho_dot_pvu, zeros(n,1), ones(n,1)];
                else
                    H = [prho_pxyz, zeros(1,4), prho_pb;
                         prho_dot_pxyz, 0, prho_dot_pvxy', 0,0];
                end
            else
                error('unsupported measurement type!')
            end
        end

        function R = R_meas(X,meas)
            % returns the covariance matrix of a particular measurement
            if strcmp(meas.type,'wheelspeed')
                R = 30; 
            elseif strcmp(meas.type,'gnss')
                % estimated value gnss log
                sigma_sq = [35,30]'.*meas.sigma;
                R = diag(reshape(sigma_sq',[],1));
            end
        end
    end         
end
