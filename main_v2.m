clear
addpath gps-measurement-tools/opensource/
load_save = false;
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
    vxy = 0;
    vz = 0;
    omega = 0;
    b = initial_state.b_sec*MeasurementModel.c;
    b_dot = 0;
    psi = 0;
    tau = 0;


    X0 = [px,py,pz,theta,vxy,vz,b,b_dot]';
    P0 = diag([3,3,3,10,1,0,5,5]);
    % minimum time step before updating state estimate
    dt = .05;

    load('wheelspeed_meas.mat')
    wheelspeed_measurements = struct('timestamp',timestamp,'omega',omega,'type','wheelspeed');
    load('gnss_measurements.mat')
    % add sat el and az data relative to initial position
    lla_init = ecef2lla(X0(1:3,1)');
    sv_pos = gnss_measurements.sv_pos_ecef;
    [sv_n,sv_e,sv_d] = ecef2ned(sv_pos(1,:),sv_pos(2,:),sv_pos(3,:),lla_init(1),lla_init(2),lla_init(3),wgs84Ellipsoid('meter'));
    sv_u = -sv_d;
    sv_el = atan2d(sv_u,vecnorm([sv_n;sv_e]));
    sv_az = atan2d(sv_e,sv_n);
    gnss_measurements.sv_el = sv_el;
    gnss_measurements.sv_az = sv_az;
    gnss_measurements.type = 'gnss';

    t0 = max([wheelspeed_measurements.timestamp(1),
        gnss_measurements.timestamp(1)]);
    gnss_t_offset = t0-gnss_measurements.timestamp(1);
    tf = t0+600*5;
    % exp: 
    %wheelspeed_measurements.timestamp = wheelspeed_measurements.timestamp+gnss_t_offset;
    tf = min([wheelspeed_measurements.timestamp(end),
        gnss_measurements.timestamp(end)]);

    [X_hat,P_hat,t_vec,meas_innovs] = simulate(X0,P0,dt,t0,tf,wheelspeed_measurements,gnss_measurements);
end

innovs = plot_innovations(meas_innovs);

check_consistency(innovs);


X = X_hat;
P = P_hat;
t = t_vec;
%[X,P,t] = condense(X_hat,P_hat,t_vec); % this does nothing now

plot_states(t,X,P);


R_enu = BicycleModel_v2.R_enu(X(1:3,1));
lla_init = ecef2lla(X_hat(1:3,1)');
sv_pos = gnss_measurements.sv_pos_ecef;
[sv_n,sv_e,sv_d] = ecef2ned(sv_pos(1,:),sv_pos(2,:),sv_pos(3,:),lla_init(1),lla_init(2),lla_init(3),wgs84Ellipsoid('meter'));
sv_u = -sv_d;
sv_el = atan2(sv_u,vecnorm([sv_n;sv_e]));
t = t-t0;
load('wlt_soln.mat')
load('wls_pvt.mat')
wls_vel = R_enu*wls_pvt.allVelMps';
wls_ven = vecnorm(wls_vel(1:2,:));

wls_soln = wlt_soln(5:end,:)';
wls_enu = R_enu*wls_soln;
wls = wls_soln-X(1:3,1);
pos_enu = R_enu*X(1:3,:);
pos_lla = ecef2lla(X_hat(1:3,:)')';

figure(1); clf;
geoscatter(wls_pvt.allLlaDegDegM(:,1),wls_pvt.allLlaDegDegM(:,2),'DisplayName','WLS')
hold on;
gplot = geoplot(pos_lla(1,:),pos_lla(2,:),'.','DisplayName','EKF');
t_row = [dataTipTextRow('time',t)];
gplot.DataTipTemplate.DataTipRows(end+1) = t_row;
%pos_enu = pos_enu-pos_enu(:,1);
%wls_enu = wls_enu-wls_enu(:,1);
%plot3(pos_enu(1,:),pos_enu(2,:),pos_enu(3,:),'.')
%scatter3(wls_enu(1,:),wls_enu(2,:),wls_enu(3,:))
%zlim([min(pos_enu(3,:)),max(pos_enu(3,:))]);
%ylim([min(pos_enu(2,:)),max(pos_enu(2,:))]);
%xlim([min(pos_enu(1,:)),max(pos_enu(1,:))]);
%xlabel('E (m)'); ylabel('N (m)'); zlabel('U (m)')
%title('ENU trajectory')
%axis equal
legend


X = X-X(:,1);

%return
% calibration of wheel speed?
wls_vel_interp = interp1(wls_pvt.FctSeconds,vecnorm(wls_vel),(0:1:(tf-t0)) + wls_pvt.FctSeconds(1)+gnss_t_offset,'spline');
whl_interp = interp1(wheelspeed_measurements.timestamp,wheelspeed_measurements.omega,t0:1:tf,'spline');
%figure(10); clf; hold on;
%plot(whl_interp*BicycleModel_v2.r_w,wls_vel_interp,'.');
%plot(0:10,0:10,'LineWidth',2);
%xlabel('wheelspeed * r')
%ylabel('snapshot')
%axis equal; axis square;
%ylim([-10,10])


% approximate rotation matrix from ecef to enu
%X_no_whl = load('X_no_whl.mat');
%X_no_whl = X_no_whl.X;
%pos_no_whl_enu = R_enu*(X_no_whl(1:3,:)+X_hat(1:3,1));
%
%figure(1); clf; hold on;
%plot(pos_enu(1,:),pos_enu(2,:),'DisplayName','EKF+wheel speed')
%%plot(X_no_whl(1,:),X_no_whl(2,:),'DisplayName','EKF')
%scatter(wls_enu(1,:),wls_enu(2,:),15,'DisplayName','WLS');
%title('XY trajectory relative to initial Position')
%xlabel('X ECEF (m)')
%ylabel('Y EcEF (m)')
%legend
%axis equal
%grid on
%figure(2); clf; hold on;
%plot(t+t0,abs(X(5,:)));
%scatter(wheelspeed_measurements.timestamp,wheelspeed_measurements.omega*BicycleModel.r_w);
%
ax = gca;
ax(1) = [];
figure(6); clf; hold on;
ax(end+1)=subplot(4,2,[1,2]); hold on;
plot(t,sqrt(squeeze(P(1,1,:))),'.','DisplayName','\sigma x')
plot(t,sqrt(squeeze(P(2,2,1:end))),'.','DisplayName','\sigma y')
plot(t,sqrt(squeeze(P(3,3,:))),'.','DisplayName','\sigma z')
title('X,Y,Z standard deviations'); ylabel('(m)'); 
legend
mean_sigma = @(i) mean(sqrt(squeeze(P(i,i,:))));
prc = @(i,p) prctile(sqrt(squeeze(P(i,i,:))),p);

fprintf('sigma x: 50pctile: %.2f, 95pctile: %.2f m \n',prc(1,50),prc(1,95));
fprintf('sigma y: 50pctile: %.2f, 95pctile: %.2f m \n',prc(2,50),prc(2,95));
fprintf('sigma z: 50pctile: %.2f, 95pctile: %.2f m \n',prc(3,50),prc(3,95));
fprintf('sigma th: 50pctile: %.2f, 95pctile: %.2f rad \n',prc(4,50),prc(4,95));
fprintf('sigma v_en: 50pctile: %.2f, 95pctile: %.2f m/s \n',prc(5,50),prc(5,95));
fprintf('sigma b_u: 50pctile: %.2f, 95pctile: %.2f m \n',prc(7,50),prc(7,95));
fprintf('sigma b_dot: 50pctile: %.2f, 95pctile: %.2f m/s \n',prc(8,50),prc(8,95));


ax(end+1)=subplot(4,2,[3,4]); hold on;
plot(t,sqrt(squeeze(P(4,4,:))),'.','DisplayName','\sigma \theta')
plot(t,sqrt(squeeze(P(5,5,:))),'.','DisplayName','\sigma v en')
%plot(t,sqrt(squeeze(P(6,6,:))),'.','DisplayName','\sigma v u')
title('\theta, V-EN, standard deviations'); 
ylabel('rad, m/s'); 
legend

ax(end+1)=subplot(4,2,[5,6]); hold on;
plot(wheelspeed_measurements.timestamp-t0, wheelspeed_measurements.omega*BicycleModel_v2.r_w,'DisplayName','speed sensor');
plot(t,abs(X(5,:)),'.','DisplayName','est V_{en}')
plot(wls_pvt.FctSeconds-wls_pvt.FctSeconds(1)-gnss_t_offset,wls_ven,'DisplayName','wls')
title('speed estimates'); ylabel('m/s')
ax(end).XLim(1) = 0;
legend

ax(end+1)=subplot(4,2,7); hold on;
s=scatter(gnss_measurements.timestamp-t0,gnss_measurements.sv_el,30,gnss_measurements.sv_ids,'.','DisplayName','el');
s.DataTipTemplate.DataTipRows(3).Label='sv id';
ylabel('el (deg)')
title('elevation')
ax(end+1)=subplot(4,2,8); hold on;
s = scatter(gnss_measurements.timestamp-t0,gnss_measurements.sv_az,30,gnss_measurements.sv_ids,'.','DisplayName','az');
s.DataTipTemplate.DataTipRows(3).Label='sv id';
svid_row = [dataTipTextRow('time',t)];
gplot.DataTipTemplate.DataTipRows(end+1) = t_row;
ylabel('az (deg)')
title('azimuth')
ax(end).YTick=-180:45:180;



xlabel('time (s)')

linkaxes(ax,'x')
ax(end).XLim(2) = t(end);
%
%% approximate rotation matrix from ecef to enu
%R_enu = Toolbox.ECEFToENU(eye(3),deg2rad(37.475),deg2rad(-122.1666));
%
%%% bridge comparison
%s_bridge = 13000;
%e_bridge = 18000;
%p_bridge = pos_enu(1:3,s_bridge:e_bridge);
%%p_no_whl = X_no_whl(1:3,s_bridge:e_bridge);
%s_wls = min(find( vecnorm(wls_enu-p_bridge(:,1)) < 50 ));
%e_wls = min(find( vecnorm(wls_enu-p_bridge(:,end)) < 50 ));
%%p_no_whl = R_enu*p_no_whl;
%%
%%
%figure(4); clf; hold on;
%plot3(p_bridge(1,:),p_bridge(2,:),p_bridge(3,:),'DisplayName','EKF+wheel speed');
%scatter(wls_enu(1,s_wls:e_wls),wls_enu(2,s_wls:e_wls),'DisplayName','WLS');
%scatter3(wls_enu(1,s_wls:e_wls),wls_enu(2,s_wls:e_wls),wls_enu(3,s_wls:e_wls),'DisplayName','WLS');
%
%title('XY Trajectory near Bridge over Hwy 101')
%xlabel('E (m)')
%ylabel('N (m)')
%legend
%axis equal
%grid on

det_P = arrayfun(@(k) det(P(:,:,k)), 1:size(P,3));

if ~load_save
    %save('save_sim')
end

function innovs = plot_innovations(meas_innov)
    whl_innov = [];
    whl_innov_cov = [];
    whl_meas = {};
    whl_filter = {};
    gnss_innov = {};
    gnss_innov_cov = {};
    gnss_meas = {};
    gnss_filter = {};


    for meas=meas_innov
        innov = meas{1}{1}{1};
        innov_cov = meas{1}{1}{2};
        m = meas{1}{1}{3};
        filter = meas{1}{1}{4};
        if size(innov_cov,1)==1
            whl_innov(end+1) = innov;
            whl_innov_cov(end+1) = innov_cov^2;
            whl_meas{end+1} = m;
            whl_filter{end+1} = filter;
        else
            gnss_innov{end+1} = innov;
            gnss_innov_cov{end+1} = innov_cov;
            gnss_meas{end+1} = m;
            gnss_filter{end+1} = filter;
        end
    end
    innovs = struct();
    innovs.whl.y = whl_innov;
    innovs.whl.S = whl_innov_cov;
    innovs.whl.meas = whl_meas;
    innovs.whl.filter = whl_filter;
    innovs.gnss.y = gnss_innov;
    innovs.gnss.S = gnss_innov_cov;
    innovs.gnss.meas = gnss_meas;
    innovs.gnss.filter = gnss_filter;

    %p_1 = .66;
    %p_2 = .95;
    
    %a=subplot(3,2,2); 
    %histogram(sqrt(whl_innov_cov(2:end)));
    %title('Innovation sqrt(cov)')
    %%a.YScale='log';
    %ylabel('wheel speed sensor (rad/s)')
    %subplot(3,2,1);
    %h = histogram(whl_innov,30);
    %h.HandleVisibility='off';
    %xlabel('\delta rad/s')
    %ylabel('wheel speed sensor')
    %title('Innovations')
    %hold on;
    %[~,w_66] = percentile(whl_innov,p_1);
    %[~,w_95] = percentile(whl_innov,p_2);
    %%plot_bar(w_66,p_1);
    %%plot_bar(w_95,p_2);
    %fprintf('wheel speed sensor: %.2i%% percentile +- %.2f, %.2i%% percentile +- %.2f\n',p_1*100,w_66,p_2*100,w_95);
    %%legend(sprintf('95%%, \\sigma=%.1f',w_95/2))
    
    %a=subplot(3,2,4); 
    %histogram(sqrt(rho_innov_cov));
    %%a.YScale='log';
    %ylabel('pseudorange (m)')
    %subplot(3,2,3);
    %xlabel('\delta m')
    %ylabel('pseudorange')
    %h=histogram(rho_innov,30);
    %h.HandleVisibility='off';
    %hold on;
    %[~,r_66] = percentile(rho_innov,p_1);
    %[~,r_95] = percentile(rho_innov,p_2);
    %%plot_bar(r_66,p_1);
    %%plot_bar(r_95,p_2);
    %fprintf('pseudorange: %.2i%% percentile +- %.2f, %.2i%% percentile +- %.2f\n',p_1*100,r_66,p_2*100,r_95);
    %%legend(sprintf('95%%, \\sigma=%.1f',r_95/2))
    
    %a=subplot(3,2,6); 
    %histogram(sqrt(rhor_innov_cov));
    %a.YScale='log';
    %ylabel('pseudorange rate (m/s)')
    %subplot(3,2,5);
    %h=histogram(rhor_innov,30);
    %h.HandleVisibility='off';
    %xlabel('\delta m/s')
    %ylabel('doppler')
    %hold on;
    %[~,rr_66] = percentile(rhor_innov,p_1);
    %[~,rr_95] = percentile(rhor_innov,p_2);
    %%plot_bar(rr_66,p_1);
    %%plot_bar(rr_95,p_2);
    %fprintf('pseudorange rate: %.2i%% percentile +- %.2f, %.2i%% percentile +- %.2f\n',p_1*100,rr_66,p_2*100,rr_95);
    %%legend(sprintf('95%%, \\sigma=%.1f',rr_95/2))
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
    
    % elevation mask in degrees
    el_mask = gnss_meas.sv_el>10;
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
            meas = struct('type','wheelspeed','z',z,'t',wheelspeed_meas.timestamp(idx_wheel));

            % update
            [X0,P0,i_cov] = ekf_update(X0_hat,P0_hat,meas,dt,t-t0);
            innov_cov{end+1} = {i_cov};
        % exp: for estimation without wheel speed 
        %X0 = X0_hat; P0 = P0_hat;
        % advance time
        idx_wheel = idx_wheel+1;
        t = t+dt_wheel;
    elseif dt_gnss<=dt_wheel
        % next measurement is gnss
        % predict
        [X0_hat,P0_hat] = ekf_predict(X0,P0,u,dt_gnss);
        % pack measurement % batch measurements together
        idxs = abs(gnss_meas.timestamp - t+dt_gnss)<.1;

        % exp: 
        % elevation mask
        %idxs = idxs & el_mask;
        
        %% filter by Cn0
        %Cn0_min = 25;
        %idxs = idxs & (gnss_meas.Cn0 > Cn0_min | isnan(gnss_meas.Cn0));
        z = [gnss_meas.pseudorange_m(idxs),...
            gnss_meas.pseudorange_rate_mps(idxs)];
        z = reshape(z,[],1);
        sigma = [gnss_meas.pseudorange_sigma(idxs);...
            gnss_meas.pseudorange_rate_sigma(idxs)];
        meas = struct('type','gnss','z',z,...
            'sigma',sigma,...
            'sv_pos',gnss_meas.sv_pos_ecef(:,idxs),...
            'sv_vel',gnss_meas.sv_vel_ecef(:,idxs),...
            'sv_B',gnss_meas.sv_clock_err_(idxs),... %lol
            'sv_B_dot',gnss_meas.sv_clock_err_rate(idxs), ...
            'sv_id',gnss_meas.sv_ids(idxs),...
            't',gnss_meas.timestamp(idx_gnss)...
            );
        % update
        [X0,P0,i_cov] = ekf_update(X0_hat,P0_hat,meas,dt,t-t0);
        innov_cov{end+1} = {i_cov};
    % advance time
    idx_gnss = idx_gnss+sum(idxs);
    t = t+dt_gnss;
end
%keyboard
X_hat_sim(:,end+1) = X0;
P_hat_sim(:,:,end+1) = P0;
t_vec(end+1) = t;

if mod(round(t-t0),200)==0 && last_disp_time ~= round(t-t0)
    fprintf('simulated elapsed time: %i\n',round(t-t0))
    last_disp_time = round(t-t0);
end
end
end

function [] = check_consistency(innovs)
    nu_whl = innovs.whl.S.^-.5.*innovs.whl.y;
    nis_whl = innovs.whl.y.*innovs.whl.S.^-1.*innovs.whl.y;
    meas_whl = [innovs.whl.meas{:}];
    t_whl = [meas_whl.t];
    n_gnss = numel(innovs.gnss.y);
    nis_rho = zeros(1,n_gnss);
    nis_rhor = zeros(1,n_gnss);
    nis_gnss_dof = zeros(1,n_gnss);
    t_gnss = zeros(1,n_gnss);
    filter_percent = zeros(1,n_gnss);
    for i=1:n_gnss
        y = innovs.gnss.y{i};
        n_y = numel(y)/2;
        nis_gnss_dof(i) = n_y;
        y_rho = y(1:n_y);
        y_rhor = y(n_y+1:end);
        S_rho = innovs.gnss.S{i}(1:n_y,1:n_y);
        S_rhor = innovs.gnss.S{i}(n_y+1:end,n_y+1:end);
        nis_rho(i) = y_rho'*S_rho^-1*y_rho;
        nis_rhor(i) = y_rhor'*S_rhor^-1*y_rhor;
        t_gnss(i) = innovs.gnss.meas{i}.t;
        filter_percent(i) = sum(innovs.gnss.filter{i})/numel(innovs.gnss.filter{i});
    end
    ax = [];
    t_start = min([t_gnss(1),t_whl(1)]);
    figure(3); clf; hold on;
    ax(end+1) = subplot(3,1,1); hold on;
    plot(t_gnss-t_start,nis_rho,'.')
    plot(t_gnss-t_start,chi2inv(.99,nis_gnss_dof),'--')
    title('nis pseudorange')
    legend('NIS','99% \chi^2 test')
    ax(end+1) = subplot(3,1,2); hold on;
    plot(t_gnss-t_start,nis_rhor,'.')
    plot(t_gnss-t_start,chi2inv(.99,nis_gnss_dof),'--')
    title('nis doppler')
    ax(end+1) = subplot(3,1,3); hold on;
    plot(t_whl-t_start,nis_whl,'.')
    plot(t_whl-t_start,chi2inv(.99,ones(size(t_whl))),'--')
    title('nis wheel speed')
    xlabel('t since start (s)')
    linkaxes(ax,'x');

    fprintf('%.2f%% of pseudoranges are within 99%% expectation\n', sum(nis_rho<chi2inv(.99,nis_gnss_dof))/n_gnss * 100)
    fprintf('%.2f%% of pseudoranges are within 65%% expectation\n', sum(nis_rho<chi2inv(.65,nis_gnss_dof))/n_gnss * 100)
    fprintf('%.2f%% of doppler are within 99%% expectation\n', sum(nis_rhor<chi2inv(.99,nis_gnss_dof))/n_gnss * 100)
    fprintf('%.2f%% of doppler are within 65%% expectation\n', sum(nis_rhor<chi2inv(.65,nis_gnss_dof))/n_gnss * 100)
    fprintf('%.2f%% of wheel speeds are within 99%% expectation\n', sum(nis_whl<chi2inv(.99,1))/numel(nis_whl) * 100)
    fprintf('%.2f%% of wheel speeds are within 65%% expectation\n', sum(nis_whl<chi2inv(.65,1))/numel(nis_whl) * 100)
end

function plot_states(t,X,P)
    labels={'X','Y','Z','theta','V_en','V_u','b','b dot'};
    ylabels={'m','m','m','rad','m/s','m/s','m','m/s'};
    figure(5); clf;
    ax = [];
    for i=1:8
        if i==6
            continue
        end
        ax(end+1) = subplot(4,2,i);  hold on;
        plot(t-t(1),X(i,:)-X(i,1),'.');
        plot(t-t(1),3*sqrt(squeeze(P(i,i,:)))' + X(i,:)-X(i,1),'.');
        plot(t-t(1),-3*sqrt(squeeze(P(i,i,:)))' + X(i,:)-X(i,1),'.');
        title(labels{i})
        ylabel(ylabels{i})
    end
    legend('estimate','+3 \sigma ','-3 \sigma')
    linkaxes(ax,'x');
end

% predict
function [X_hat,P_hat] = ekf_predict(X,P,u,dt)
    X_hat = BicycleModel_v2.f_dyn(X,u,dt);
    F = BicycleModel_v2.F_jac(X,u,dt);
    Q = BicycleModel_v2.Q_dyn(X,u,dt);
    P_hat = F*P*F'+Q;
end
% update
function [X,P,innov] = ekf_update(X_hat,P_hat,meas,dt,t)
    if strcmp(meas.type,'wheelspeed')
        [X,P,innov] = ekf_update_step(X_hat,P_hat,meas,dt,t);
        return
    else
        dtof = inf;
        tof_old = zeros(size(meas.sv_id));
        meas.tof = zeros(size(meas.sv_id));

        while (dtof) > 1e-8
            [X,P,innov] = ekf_update_step(X_hat,P_hat,meas,dt,t);
            pos = X(1:3);

            tof = vecnorm(meas.sv_pos-pos)/MeasurementModel.c;
            dtof = norm(tof-tof_old);
            tof_old = tof;
            meas.tof = tof;
        end
    end
end

function [X,P,innov] = ekf_update_step(X_hat,P_hat,meas,dt,t)
    z = meas.z;
    y_tilde = z-MeasurementModel.h_meas(X_hat,meas);
    H = MeasurementModel.H_jac(X_hat,meas);
    R = MeasurementModel.R_meas(X_hat,meas);
    innov_cov = R+H*P_hat*H';
    i = 1;
    i_kept = [];
    innov = {y_tilde,innov_cov,meas};
    n_prev = length(y_tilde);
    % exp: 
    % not filter if few measurements or within 1 minute of starting
    %while (i<=length(y_tilde)) && (length(y_tilde)>4) && t > 60
    while (i<=length(y_tilde)) &&  t > 60
        if abs(y_tilde(i))>15*sqrt(R(i,i))
        %if abs(y_tilde(i))>30*sqrt(R(i,i))
            i_kept(end+1) = 0;
            % exp: disable filter
            y_tilde(i) = [];
            R(i,:) = []; R(:,i) = [];
            H(i,:) = []; 
            %i = i+1;
        else
            i_kept(end+1) = 1;
            i = i+1;
        end
    end
    innov{end+1} = i_kept;
    if n_prev-length(y_tilde) > 5
        fprintf('removed %i meas\n',n_prev-length(y_tilde))
    end
    %innov = {y_tilde,innov_cov,meas};
    K = P_hat*H'*(R+H*P_hat*H')^-1;
    X = X_hat+K*y_tilde;

    %if (t>1300); keyboard; end

    %X(4) = wrapToPi(X(4));
    %        if X(5)<0
    %            X(5) = -X(5);
    %            X(4) = X(4)+pi;
    %        end
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
