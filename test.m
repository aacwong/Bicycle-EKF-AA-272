clear

load('initial_state.mat')
px = initial_state.pos_ecef(1);
py = initial_state.pos_ecef(2);
pz = initial_state.pos_ecef(3);
theta = 0;
ven = 1;
vu = 0;
b = initial_state.b_sec*MeasurementModel.c;
b_dot = 0;
psi = 0;
tau = 0;


X0 = [px,py,pz,theta,ven,vu,b,b_dot]';
hist = X0;
dt = .05;
tvec = 0:dt:10;
P = eye(8);
det_P_hist = det(P);

for t=tvec(2:end)
    [hist(:,end+1),P] = ekf_predict(hist(:,end),P,[0;0],dt);
    det_P_hist(end+1) = det(P);
end
hist = hist-hist(:,1);
plot(tvec, (hist(1:3,:))' )
legend

function [X_hat,P_hat] = ekf_predict(X,P,u,dt)
    X_hat = BicycleModel_v2.f_dyn(X,u,dt);
    F = BicycleModel_v2.F_jac(X,u,dt);
    Q = BicycleModel_v2.Q_dyn(X,u,dt);
    P_hat = F*P*F'+Q;
                        if any(isnan(P_hat),'all') || any(isinf(P_hat),'all'); keyboard; end;
end
