clear
addpath gps-measurement-tools/opensource/
load_save = false
set(groot,'defaultAxesXGrid','on')
set(groot,'defaultAxesYGrid','on')

if load_save
    load('save_sim.mat')
else

    load('initial_state.mat')
    px = initial_state.pos_ecef(1);
    py = initial_state.pos_ecef(2);
    pz = initial_state.pos_ecef(3);
    theta = 0;
    vx = 0;
    vy = 0;
    omega = 0;
    b = initial_state.b_sec*MeasurementModel.c;
    psi = 0;
    tau = 0;


    X0 = [px,py,pz,theta,vx,vy,omega,b]';
    P0 = diag(ones(8,1));
    % minimum time step before updating state estimate
    dt = .05;

    load('wheelspeed_meas.mat')
    wheelspeed_measurements = struct('timestamp',timestamp,'omega',omega,'type','wheelspeed');
    load('gnss_measurements.mat')
    gnss_measurements.type = 'gnss';

    t0 = max([wheelspeed_measurements.timestamp(1),
        gnss_measurements.timestamp(1)]);
    %tf = t0+600*2;
    tf = min([wheelspeed_measurements.timestamp(end),
        gnss_measurements.timestamp(end)]);


    [X_hat,P_hat,t_vec,innovs] = simulate(X0,P0,dt,t0,tf,wheelspeed_measurements,gnss_measurements);
end

figure(3); clf; hold on;
[whl_innov,whl_cov,rho_innov,rho_cov,rhor_innov,rhor_cov] = plot_innovations(innovs);

[X,P,t] = condense(X_hat,P_hat,t_vec);

R_enu = Toolbox.ECEFToENU(eye(3),deg2rad(37.475),deg2rad(-122.1666));
t = t-t0;
load('wlt_soln.mat')
wls_soln = wlt_soln';
wls_soln_enu = R_enu*wls_soln;
wls = wls_soln-X(1:3,1);
pos_enu = R_enu*X(1:3,:);
X = X-X(:,1);
% approximate rotation matrix from ecef to enu
X_no_whl = load('X_no_whl.mat');
X_no_whl = X_no_whl.X;
pos_no_whl_enu = R_enu*(X_no_whl(1:3,:)+X_hat(1:3,1));

figure(1); clf; hold on;
plot(X(1,:),X(2,:),'DisplayName','EKF+wheel speed')
plot(X_no_whl(1,:),X_no_whl(2,:),'DisplayName','EKF')
scatter(wls(1,:),wls(2,:),15,'DisplayName','WLS');
title('XY trajectory relative to initial Position')
xlabel('X ECEF (m)')
ylabel('Y EcEF (m)')
legend
axis equal
grid on

% approximate rotation matrix from ecef to enu
R_enu = Toolbox.ECEFToENU(eye(3),deg2rad(37.475),deg2rad(-122.1666));

% bridge comparison
s_bridge = 13000;
e_bridge = 18000;
p_bridge = X(1:3,s_bridge:e_bridge);
p_no_whl = X_no_whl(1:3,s_bridge:e_bridge);
[~,s_wls] = min(vecnorm(wls-p_bridge(:,1)));
[~,e_wls] = min(vecnorm(wls-p_bridge(:,end)));
p_bridge = R_enu*p_bridge;
p_no_whl = R_enu*p_no_whl;
p_wls = R_enu*wls(:,s_wls:e_wls);
det_P = arrayfun(@(k) det(P(:,:,k)), 1:size(P,3));

figure(4); clf; hold on;
plot(p_bridge(1,:),p_bridge(2,:),'DisplayName','EKF+wheel speed');
scatter(p_wls(1,:),p_wls(2,:),'DisplayName','WLS');

title('XY Trajectory near Bridge over Hwy 101')
xlabel('E (m)')
ylabel('N (m)')
legend
axis equal
grid on


if ~load_save
    %save('save_sim')
end

function [whl_innov,whl_innov_cov,rho_innov,rho_innov_cov,rhor_innov,rhor_innov_cov] = plot_innovations(meas_innov)
    whl_innov = [];
    whl_innov_cov = [];
    rho_innov = [];
    rhor_innov = [];
    rho_innov_cov = [];
    rhor_innov_cov = [];

    for meas=meas_innov
        innov = meas{1}{1}{1};
        innov_cov = meas{1}{1}{2};
        if size(innov_cov,1)==1
            whl_innov(end+1) = innov;
            whl_innov_cov(end+1) = innov_cov^2;
        else
            rho_innov(end+1) = innov(1);
            rho_innov_cov(end+1) = innov_cov(1,1)^2;
            rhor_innov(end+1) = innov(2);
            rhor_innov_cov(end+1) = innov_cov(2,2)^2;
        end
    end

    plot_bar = @(w,p) plot(-w*[-1,1],diff(ylim)*(1-p)*[1,1],'lineWidth',3);
    p_1 = .66;
    p_2 = .95;

    a=subplot(3,2,2); 
    plot(abs(whl_innov_cov));
    title('Innovation Covariances')
    a.YScale='log';
    ylabel('wheel speed sensor (rad^2/s^2)')
    subplot(3,2,1);
    h = histogram(whl_innov,30);
    h.HandleVisibility='off';
    xlabel('\delta rad/s')
    ylabel('wheel speed sensor')
    title('Innovations')
    hold on;
    [~,w_66] = percentile(whl_innov,p_1);
    [~,w_95] = percentile(whl_innov,p_2);
    %plot_bar(w_66,p_1);
    plot_bar(w_95,p_2);
    fprintf('wheel speed sensor: %.2i%% percentile +- %.2f, %.2i%% percentile +- %.2f\n',p_1*100,w_66,p_2*100,w_95);
    legend(sprintf('95%%, \\sigma=%.1f',w_95/2))

    a=subplot(3,2,4); 
    plot(abs(rho_innov_cov));
    a.YScale='log';
    ylabel('pseudorange (m^2)')
    subplot(3,2,3);
    xlabel('\delta m')
    ylabel('pseudorange')
    h=histogram(rho_innov,30);
    h.HandleVisibility='off';
    hold on;
    [~,r_66] = percentile(rho_innov,p_1);
    [~,r_95] = percentile(rho_innov,p_2);
    %plot_bar(r_66,p_1);
    plot_bar(r_95,p_2);
    fprintf('pseudorange: %.2i%% percentile +- %.2f, %.2i%% percentile +- %.2f\n',p_1*100,r_66,p_2*100,r_95);
    legend(sprintf('95%%, \\sigma=%.1f',r_95/2))

    a=subplot(3,2,6); 
    plot(abs(rhor_innov));
    a.YScale='log';
    ylabel('pseudorange rate (m^2/s^2)')
    subplot(3,2,5);
    h=histogram(rhor_innov,30);
    h.HandleVisibility='off';
    xlabel('\delta m/s')
    ylabel('doppler')
    hold on;
    [~,rr_66] = percentile(rhor_innov,p_1);
    [~,rr_95] = percentile(rhor_innov,p_2);
    %plot_bar(rr_66,p_1);
    plot_bar(rr_95,p_2);
    fprintf('pseudorange rate: %.2i%% percentile +- %.2f, %.2i%% percentile +- %.2f\n',p_1*100,rr_66,p_2*100,rr_95);
    legend(sprintf('95%%, \\sigma=%.1f',rr_95/2))
end
function [X,P,t] = condense(X_vec,P_vec,t_vec)
    % removes repeated time steps in the ekf output due to simulatneous update steps at the same
    [t,idx,~] = unique(t_vec,'last');
    X = X_vec(:,idx);
    P = P_vec(:,:,idx);
end

function [X_hat_sim,P_hat_sim,t_vec,innov_cov] = ...
        simulate(X0,P0,dt,t0,tf,wheelspeed_meas,gnss_meas)
    t = t0;
    X_hat_sim = [];
    P_hat_sim = [];
    innov_cov = {};
    t_vec = t;
    X_hat_sim(:,1) = X0;
    P_hat_sim(:,:,1) = P0;
    % assume constant zero control input
    u = zeros(2,1);

    % X0 is the state estimate
    % X0_hat is the intermediate estimate before measurement

    ts_wheel = wheelspeed_meas.timestamp;
    ts_wheel = ts_wheel(ts_wheel>=t0);
    idx_wheel = min(find(ts_wheel>=t0));
    ts_gnss = gnss_meas.timestamp;
    idx_gnss = min(find(ts_gnss>=t0));
    last_disp_time = 0;
    while t <= tf
        % calculate next time step [measurement or update]
        if idx_wheel> length(ts_wheel)
            dt_wheel = inf;
        else
            dt_wheel = ts_wheel(idx_wheel)-t;
        end
        if idx_gnss> length(ts_gnss)
            dt_gnss = inf;
        else
            dt_gnss = ts_gnss(idx_gnss)-t;
        end
        if all([dt_wheel,dt_gnss] > dt)
            % time to next measurement is too big, need to update state estimate without measurement
            [X0,P0] = ekf_predict(X0,P0,u,dt);
            t = t+dt;
        elseif dt_wheel<dt_gnss
            [X0,P0] = ekf_predict(X0,P0,u,dt);
            % next measurement is a wheel speed meas
            % predict
            [X0_hat,P0_hat] = ekf_predict(X0,P0,u,dt_wheel);
            % pack measurement
            z = wheelspeed_meas.omega(idx_wheel);
            meas = struct('type','wheelspeed','z',z);
            % update
            [X0,P0,i_cov] = ekf_update(X0_hat,P0_hat,meas,dt);
            innov_cov{end+1} = {i_cov};
            % for estimation without wheel speed 
            %X0 = X0_hat; P0 = P0_hat;
            % advance time
            idx_wheel = idx_wheel+1;
            t = t+dt_wheel;
        elseif dt_gnss<=dt_wheel
            % next measurement is gnss
            % predict
            [X0_hat,P0_hat] = ekf_predict(X0,P0,u,dt_gnss);
            % pack measurement
            z = [gnss_meas.pseudorange_m(idx_gnss);...
                gnss_meas.pseudorange_rate_mps(idx_gnss)];
            meas = struct('type','gnss','z',z,...
            'sigma',[gnss_meas.pseudorange_sigma(idx_gnss);...
                gnss_meas.pseudorange_rate_sigma(idx_gnss)],...
            'sv_pos',gnss_meas.sv_pos_ecef(:,idx_gnss),...
            'sv_vel',gnss_meas.sv_vel_ecef(:,idx_gnss),...
            'sv_B',gnss_meas.sv_clock_err_(idx_gnss),... %lol
            'sv_B_dot',gnss_meas.sv_clock_err_rate(idx_gnss), ...
            'sv_id',gnss_meas.sv_ids(idx_gnss)...
            );
            % update
            [X0,P0,i_cov] = ekf_update(X0_hat,P0_hat,meas,dt);
            innov_cov{end+1} = {i_cov};
            % advance time
            idx_gnss = idx_gnss+1;
            t = t+dt_gnss;
        end
        X_hat_sim(:,end+1) = X0;
        P_hat_sim(:,:,end+1) = P0;
        t_vec(end+1) = t;

        if mod(round(t-t0),100)==0 && last_disp_time ~= round(t-t0)
            fprintf('simulated elapsed time: %i\n',round(t-t0))
            last_disp_time = round(t-t0);
        end
    end
end


% predict
function [mu_hat,P_hat] = ekf_predict(mu,P,u,dt)
    mu_hat = BicycleModel.f_dyn(mu,u,dt);
    F = BicycleModel.F_jac(mu,u,dt);
    Q = BicycleModel.Q_dyn(mu,u,dt);
    P_hat = F*P*F'+Q;
end
% update
function [mu,P,innov] = ekf_update(mu_hat,P_hat,meas,dt)
    if strcmp(meas.type,'wheelspeed')
        [mu,P,innov] = ekf_update_step(mu_hat,P_hat,meas,dt);
        return
    else
        dtof = inf;
        tof_old = 0;
        meas.tof = 0;

        while dtof > 1e-5
            [mu,P,innov] = ekf_update_step(mu_hat,P_hat,meas,dt);
            pos = mu(1:3);

            tof = norm(meas.sv_pos-pos)/MeasurementModel.c;
            dtof = abs(tof-tof_old);
            tof_old = tof;
            meas.tof = tof;
        end
    end
end
function [mu,P,innov] = ekf_update_step(mu_hat,P_hat,meas,dt)
    %%mu = mu_hat; P = P_hat; return
    z = meas.z;
    y_tilde = z-MeasurementModel.h_meas(mu_hat,meas);
    H = MeasurementModel.H_jac(mu_hat,meas);
    R = MeasurementModel.R_meas(mu_hat,meas);
    innov_cov = R+H*P_hat*H';
    innov = {y_tilde,innov_cov,meas};
    K = P_hat*H'*(R+H*P_hat*H')^-1;
    if abs(trace(innov_cov))>1e6
        disp('very high innovation!')
        innov{1} = 0*y_tilde;
        innov{2} = 0*innov_cov;
        if strcmp(meas.type,'gnss')
            fprintf('sv: %i\n',meas.sv_id) ;
        end
    end
    mu = mu_hat+K*y_tilde;
    P = (eye(size(P_hat)) - K*H)*P_hat;
end

function [avg,dev] = percentile(vals,percent)
    % gets the +- bounds from mean of vals that percent of vals are within
    avg = mean(vals);
    abs_vals = sort(abs(vals-avg));
    n = length(abs_vals);
    nperc = round(n*percent);
    dev = abs_vals(nperc);
end

function [s_wls,e_wls] = wls_idx(path,wls)
    [~,s_wls] = min(vecnorm(wls-path(:,1)));
    [~,e_wls] = min(vecnorm(wls-path(:,end)));
end

