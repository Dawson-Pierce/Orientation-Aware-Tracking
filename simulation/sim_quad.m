function [t, p_out, R_out, q_out, X] = sim_quad(t_cmd, p_cmd, P)
%SIM_QUAD  High-fidelity quadrotor dynamics tracking a reference trajectory.
%
%   [t, p_out, R_out, q_out, X] = sim_quad(t_cmd, p_cmd)
%   [t, p_out, R_out, q_out, X] = sim_quad(t_cmd, p_cmd, P)
%
%   Inputs
%     t_cmd  -- N-by-1 time vector for reference positions
%     p_cmd  -- N-by-3 reference positions to track
%     P      -- (optional) quadrotor parameter struct; uses defaults if omitted
%
%   Outputs
%     t      -- M-by-1 internal time vector (fine dt)
%     p_out  -- M-by-3 actual positions
%     R_out  -- 3-by-3-by-M rotation matrices
%     q_out  -- M-by-4 quaternions [qw qx qy qz]
%     X      -- M-by-17 full state history

t_cmd = double(t_cmd(:));
p_cmd = double(p_cmd);
if size(p_cmd,1) ~= numel(t_cmd) || size(p_cmd,2) ~= 3
    error('p_cmd must be N-by-3 and t_cmd must be N-by-1')
end

%% Default parameters (matching sim_uav.m)
if nargin < 3 || isempty(P)
    P = struct();
end
if ~isfield(P,'m'),        P.m = 1.4;                          end
if ~isfield(P,'g'),        P.g = 9.80665;                      end
if ~isfield(P,'J'),        P.J = diag([0.03 0.03 0.05]);       end
if ~isfield(P,'arm'),      P.arm = 0.22;                       end

L = P.arm/sqrt(2);
if ~isfield(P,'r')
    P.r = [ L  L 0;
            L -L 0;
           -L -L 0;
           -L  L 0 ]';
end

if ~isfield(P,'kf'),       P.kf = 8.54858e-6;                  end
if ~isfield(P,'km'),       P.km = 1.6e-7;                      end
if ~isfield(P,'spin'),     P.spin = [1; -1; 1; -1];            end
if ~isfield(P,'cTau'),     P.cTau = P.km / P.kf;               end

if ~isfield(P,'tauOmega'), P.tauOmega = 0.04;                  end
if ~isfield(P,'OmegaMin'), P.OmegaMin = 0;                     end
if ~isfield(P,'OmegaMax'), P.OmegaMax = 1000;                  end

if ~isfield(P,'aMax'),     P.aMax = 6.0;                       end
if ~isfield(P,'tauMax'),   P.tauMax = [2.0; 2.0; 0.8];         end
if ~isfield(P,'u1Min'),    P.u1Min = 0;                        end
if ~isfield(P,'u1Max'),    P.u1Max = 2.2*P.m*P.g;              end

if ~isfield(P,'Kp'),       P.Kp = diag([2.0 2.0 6.0]);         end
if ~isfield(P,'Kd'),       P.Kd = diag([2.8 2.8 4.5]);         end
if ~isfield(P,'KR'),       P.KR = diag([8.0 8.0 1.8]);         end
if ~isfield(P,'Kw'),       P.Kw = diag([0.55 0.55 0.30]);      end

%% Build fine time grid and reference signals
dt = min(diff(t_cmd));
t  = (t_cmd(1):dt:t_cmd(end))';
p_ref = interp1(t_cmd, p_cmd, t, 'pchip');

v_ref = fdiff(p_ref, dt);
a_ref = fdiff(v_ref, dt);

%% Initial state  [p(3) v(3) q(4) w(3) Omega(4)]
%  Trim to hover at the initial reference: match position, velocity, and
%  the attitude required to produce the feed-forward acceleration.
x = zeros(17,1);
x(1:3) = p_ref(1,:)';
x(4:6) = v_ref(1,:)';

F0  = P.m*a_ref(1,:)' + [0;0;P.m*P.g];
R0  = R_from_Fyaw(F0, 0);
x(7:10) = R_to_quat(R0);

u1_0 = norm(F0);
x(14:17) = sqrt(max(u1_0/4, 0) / P.kf) * ones(4,1);

M = numel(t);
X = zeros(M,17);
X(1,:) = x';

%% RK4 integration
for k = 1:M-1
    tk = t(k);
    pr = p_ref(k,:)';
    vr = v_ref(k,:)';
    ar = a_ref(k,:)';

    k1 = quad_step(tk,           x,              pr, vr, ar, P);
    k2 = quad_step(tk+0.5*dt,    x+0.5*dt*k1,    pr, vr, ar, P);
    k3 = quad_step(tk+0.5*dt,    x+0.5*dt*k2,    pr, vr, ar, P);
    k4 = quad_step(tk+dt,        x+dt*k3,        pr, vr, ar, P);

    x = x + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
    x(7:10) = x(7:10) / norm(x(7:10));

    X(k+1,:) = x';
end

%% Pack outputs
p_out = X(:,1:3);
q_out = X(:,7:10);

R_out = zeros(3,3,M);
for k = 1:M
    R_out(:,:,k) = quat_to_R(q_out(k,:)');
end
end

%% ---- Local functions ----

function dx = quad_step(~, x, pr, vr, ar, P)
p  = x(1:3);
v  = x(4:6);
q  = x(7:10);
w  = x(11:13);
Om = x(14:17);

R  = quat_to_R(q);
b3 = R(:,3);

ep = pr - p;
ev = vr - v;

a_des = ar + P.Kp*ep + P.Kd*ev;
a_des = satvec(a_des, P.aMax);

F_des = P.m*a_des + [0;0;P.m*P.g];

yaw_des = 0;
R_des = R_from_Fyaw(F_des, yaw_des);

eR = 0.5*vee(R_des'*R - R'*R_des);
ew = w;

u1 = F_des.'*b3;
u1 = min(max(u1, P.u1Min), P.u1Max);

tau = -P.KR*eR - P.Kw*ew + cross(w, P.J*w);
tau = min(max(tau, -P.tauMax), P.tauMax);

u = [u1; tau];
f_cmd = mixer_from_geometry(u, P.r, P.cTau, P.spin);
f_cmd = max(f_cmd, 0);

Om_cmd = sqrt(f_cmd / P.kf);
Om_cmd = min(max(Om_cmd, P.OmegaMin), P.OmegaMax);
Om_dot = (Om_cmd - Om)/P.tauOmega;

f = P.kf*(Om.^2);
u1_act = sum(f);

tau_act = zeros(3,1);
for i = 1:4
    tau_act = tau_act + cross(P.r(:,i), [0;0;f(i)]);
end
tau_act(3) = tau_act(3) + sum(P.spin .* (P.cTau*f));

p_dot = v;
v_dot = (1/P.m) * (u1_act * b3) + [0;0;-P.g];
q_dot = 0.5 * quat_omega_mat(w) * q;
w_dot = P.J \ (tau_act - cross(w, P.J*w));

dx = [p_dot; v_dot; q_dot; w_dot; Om_dot];
end

function y = fdiff(x, dt)
y = zeros(size(x));
y(2:end-1,:) = (x(3:end,:) - x(1:end-2,:)) / (2*dt);
y(1,:)       = (x(2,:) - x(1,:)) / dt;
y(end,:)     = (x(end,:) - x(end-1,:)) / dt;
end

function v = satvec(v, vmax)
n = norm(v);
if n > vmax && n > 0
    v = (vmax/n)*v;
end
end

function R = R_from_Fyaw(F, yaw)
T = norm(F);
if T < 1e-9
    b3 = [0;0;1];
else
    b3 = F/T;
end
cpsi = cos(yaw);  spsi = sin(yaw);
b1c = [cpsi; spsi; 0];
b2 = cross(b3, b1c);
nb2 = norm(b2);
if nb2 < 1e-9
    b1c = [1;0;0];
    b2  = cross(b3, b1c);
    nb2 = norm(b2);
end
b2 = b2/nb2;
b1 = cross(b2, b3);
R  = [b1 b2 b3];
end

function f = mixer_from_geometry(u, r, cTau, spin)
B = zeros(4,4);
B(1,:) = 1;
for i = 1:4
    mi = cross(r(:,i), [0;0;1]);
    B(2,i) = mi(1);
    B(3,i) = mi(2);
    B(4,i) = mi(3) + spin(i)*cTau;
end
f = B \ u;
end

function R = quat_to_R(q)
q  = q / norm(q);
qw = q(1); qx = q(2); qy = q(3); qz = q(4);
R  = [1-2*(qy^2+qz^2)   2*(qx*qy-qw*qz)   2*(qx*qz+qw*qy);
      2*(qx*qy+qw*qz)   1-2*(qx^2+qz^2)   2*(qy*qz-qw*qx);
      2*(qx*qz-qw*qy)   2*(qy*qz+qw*qx)   1-2*(qx^2+qy^2)];
end

function Q = quat_omega_mat(w)
wx = w(1); wy = w(2); wz = w(3);
Q = [ 0   -wx -wy -wz;
      wx   0   wz -wy;
      wy  -wz  0   wx;
      wz   wy -wx  0  ];
end

function q = R_to_quat(R)
%R_TO_QUAT  Rotation matrix to unit quaternion [qw qx qy qz].
%  Shepperd's method — numerically robust for all orientations.
tr = trace(R);
d  = [R(1,1); R(2,2); R(3,3); tr];
[~, idx] = max(d);
switch idx
    case 4
        s  = 2*sqrt(1 + tr);
        qw = 0.25*s;
        qx = (R(3,2)-R(2,3))/s;
        qy = (R(1,3)-R(3,1))/s;
        qz = (R(2,1)-R(1,2))/s;
    case 1
        s  = 2*sqrt(1 + R(1,1) - R(2,2) - R(3,3));
        qw = (R(3,2)-R(2,3))/s;
        qx = 0.25*s;
        qy = (R(1,2)+R(2,1))/s;
        qz = (R(1,3)+R(3,1))/s;
    case 2
        s  = 2*sqrt(1 - R(1,1) + R(2,2) - R(3,3));
        qw = (R(1,3)-R(3,1))/s;
        qx = (R(1,2)+R(2,1))/s;
        qy = 0.25*s;
        qz = (R(2,3)+R(3,2))/s;
    case 3
        s  = 2*sqrt(1 - R(1,1) - R(2,2) + R(3,3));
        qw = (R(2,1)-R(1,2))/s;
        qx = (R(1,3)+R(3,1))/s;
        qy = (R(2,3)+R(3,2))/s;
        qz = 0.25*s;
end
q = [qw; qx; qy; qz];
if qw < 0, q = -q; end
end

function v = vee(M)
v = [M(3,2); M(1,3); M(2,1)];
end
