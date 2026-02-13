%% Multi-Sensor Point Cloud simulations
%
% Author: Dawson Pierce

clear; clc; close all

addpath('../../brew/')

pngDir = 'png_folder/';

if ~exist(pngDir,"dir")
    mkdir(pngDir);
end

%% Import STL for TR obj

TRsingle = stlread("drone_shape.stl");

num_drones = 10;

TR = repmat({TRsingle}, 1, num_drones);

%% Setup motion model / ETT PHD

target_motion = BREW.dynamics.Integrator_3D();

% generate measurements 

dt = 0.2;
tf = 10; 

t = 0:dt:tf;

measurements = {};

g = 9.81;
tilt_gain = 1.0;
min_speed = 0.2;

RzRyRx = @(yaw,pitch,roll) [ ...
    cos(yaw)*cos(pitch),  cos(yaw)*sin(pitch)*sin(roll)-sin(yaw)*cos(roll),  cos(yaw)*sin(pitch)*cos(roll)+sin(yaw)*sin(roll); ...
    sin(yaw)*cos(pitch),  sin(yaw)*sin(pitch)*sin(roll)+cos(yaw)*cos(roll),  sin(yaw)*sin(pitch)*cos(roll)-cos(yaw)*sin(roll); ...
    -sin(pitch),          cos(pitch)*sin(roll),                               cos(pitch)*cos(roll) ];

% Map world accel -> small-angle pitch/roll (quad tilts to accelerate)
% Pitch: nose-down positive ~ -atan2(a_x, g), Roll: right-wing-down ~ atan2(a_y, g)
tiltFromAccel = @(a) deal( ...
    atan2(norm(a(1:2)), g), ...    % total tilt magnitude (not used further)
    -tilt_gain * atan2(a(1), g), ... % pitch
     tilt_gain * atan2(a(2), g) );   % roll

paths = cell(1, num_drones);

paths{1}  = makeCirclePath([10,15, 2], 6, 0.6);
paths{2}  = makeCirclePath([20,15, 1], 5, 0.75);
paths{3}  = makeLemniscatePath([15,15, 2], 6, 0.5);
paths{4}  = makeLissajousPath([15,15, 1], [7,5], [1,2], 0.5);
paths{5}  = makeSineWeavePath([5, 5,  2], [20, 0, 0], [3, 0, 0], 0.8);
paths{6}  = makeSineWeavePath([5, 10, 1], [0, 20, 0], [0, 3, 0], 0.8);
paths{7}  = makeRoundedRectPath([8,8,  2], [22,22, 2], 3, 0.5);
paths{8}  = makeEllipsePath([15,14, 1.5], [8,4], 0.6);
paths{9}  = makeCircleWithZPath([12,18, 2], 5.5, 0.65, 1.0);
paths{10} = makeWanderPath([22, 8, 2], [0.5, 0.4, 0.2], 0.3);


active_paths = [1];

% Subset TR and paths to the active set
idx = active_paths(:)';        % row vector
TR_active = TR(idx);
paths_active = paths(idx);
num_active = numel(idx);

T = cell(num_active,1);

for kAct = 1:num_active
    [p0, v0, a0] = evalPVA(paths_active{kAct}, 0);
    vxy = [v0(1); v0(2)];
    if norm(vxy) < min_speed, vxy = [min_speed; 0]; end
    yaw0 = atan2(vxy(2), vxy(1));
    [~, pitch0, roll0] = tiltFromAccel(a0);
    R0 = RzRyRx(yaw0, pitch0, roll0);
    T{kAct} = [R0, p0(:); 0 0 0 1];
end


% initialize filter(s)
H = [eye(3) zeros(3, 3)];
R = 0.01 * eye(3); % Measurement noise
process_noise = 1 * eye(6); 
% inner_filter = BREW.filters.TrajectoryGGIWEKF('dyn_obj',target_motion, ...
%     'H',H,'process_noise',process_noise,'measurement_noise', R, 'L',20);

inner_filter = BREW.filters.GGIWEKF('dyn_obj',target_motion, ...
    'H',H,'process_noise',process_noise,'measurement_noise', R);



alpha = {10};
beta = {1}; 
mean = {[17.5; 17.5; 0; 0; 0; 0]}; 
covariance = {diag([15,15,10,2,2,2])}; 
IWdof = {10}; 
IWshape = {[1 0 0; 0 1 0; 0 0 1]}; 
weight = [1];

% birth_model = BREW.distributions.TrajectoryGGIWMixture( ...
%     'idx', {1}, ...
%     'alphas', alpha, ...
%     'betas', beta, ...
%     'means', mean, ...
%     'covariances', covariance, ...
%     'IWdofs', IWdof, ...
%     'IWshapes', IWshape, ...
%     'weights', weight);

birth_model = BREW.distributions.GGIWMixture( ...
    'alphas', alpha, ...
    'betas', beta, ...
    'means', mean, ...
    'covariances', covariance, ...
    'IWdofs', IWdof, ...
    'IWshapes', IWshape, ...
    'weights', weight);

% phd = BREW.multi_target.PHD('filter',inner_filter, 'birth_model', birth_model,...
%     'prob_detection', 0.8, 'prob_survive', 0.9, 'max_terms',50, ...
%     'cluster_obj',BREW.clustering.DBSCAN_obj(3,5),'extract_threshold',0.5);

%% Setup sensor_station objs

sensors{1} = sensor_station('location',[0 20 0],'view_vector',[1 0 0],'grid_slope',0.03);
% sensors{2} = sensor_station('location',[20 0 0],'view_vector',[0 1 0],'grid_slope',0.03);

% Setting up internal trackers this way so we can create copies since the
% classes for tracking inherit the handle class
for k = 1:length(sensors)
    sensors{k}.internal_tracker = drone_PHD( ...
        'filter',copy(inner_filter), ...
        'birth_model', copy(birth_model),...
    'prob_detection', 0.8, ...
    'prob_survive', 0.9, ...
    'max_terms',50, ...
    'cluster_obj',BREW.clustering.DBSCAN_obj(3,5), ...
    'extract_threshold',0.3);
end

points = cell(length(sensors),length(t)); 

%% Simulation

est_mix = cell(length(sensors),1);
est_colors = {'m','c'};

gifFilename = 'PointCloudDroneTrackingPOC.gif';

plot = 1;

if plot == 0 
   for k = 1:length(t) 
        for kAct = 1:num_active
            [p, v, a] = evalPVA(paths_active{kAct}, t(k));

            vxy = [v(1); v(2)];
            if norm(vxy) < min_speed, vxy = [min_speed; 0]; end
            yaw = atan2(vxy(2), vxy(1));

            [~, pitch, roll] = tiltFromAccel(a);
            R = RzRyRx(yaw, pitch, roll);

            T{kAct} = [R, p(:); 0 0 0 1];
        end
    end
elseif plot == 1
    f = figure; 
    subplot(1,2,1)
    ax = axes;
    axis equal; 
    view(3); hold on
    
    subplot(1,2,2)
    ax2 = axes;
    axis equal; 
    view(2); hold on

    for k = 1:length(t) 
        for kAct = 1:num_active
            [p, v, a] = evalPVA(paths_active{kAct}, t(k));

            vxy = [v(1); v(2)];
            if norm(vxy) < min_speed, vxy = [min_speed; 0]; end
            yaw = atan2(vxy(2), vxy(1));

            [~, pitch, roll] = tiltFromAccel(a);
            R = RzRyRx(yaw, pitch, roll);

            T{kAct} = [R, p(:); 0 0 0 1];
        end
        subplot(1,2,1)
        cla 
        trimmed_meas = cell(num_active);
        for kk = 1:length(sensors)
            points{kk,k} = sensors{kk}.get_points(TR_active,T);
            scatter3(points{kk,k}(:,1),points{kk,k}(:,2),points{kk,k}(:,3),6,'*'); hold on
            % sensors{kk}.plot()
            
            sensors{kk}.internal_tracker.predict(dt);
            [sensors{kk}.internal_tracker, trimmed_meas{kk}] = sensors{kk}.internal_tracker.correct(dt,points{kk,k}','basis_tracking',true);
            
            est_mix{kk} = sensors{kk}.internal_tracker.cleanup();

            for ii = 1:length(est_mix{kk})
                % temp_shape_mat = est_mix{kk}.distributions{ii}.IWshapes(:,:,end);
                % [B,D,Bp] = eig(temp_shape_mat);
                % p = est_mix{kk}.distributions{ii}.mean_history(1:3,end);
                % 
                % B(:,diag(D)<0)=-B(:,diag(D)<0);
                % D = abs(D);
                % 
                % [val,idx] = sort(diag(D));

                % B = B(:,idx);

                p = est_mix{kk}.distributions{ii}.mean(1:3);

                B = est_mix{kk}.distributions{ii}.basis;

                scale = 5;

                plot3([p(1) p(1)+scale*B(1,1)], [p(2) p(2)+scale*B(2,1)], [p(3) p(3)+scale*B(3,1)], 'r', 'LineWidth', 2);
                plot3([p(1) p(1)+scale*B(1,2)], [p(2) p(2)+scale*B(2,2)], [p(3) p(3)+scale*B(3,2)], 'g', 'LineWidth', 2);
                plot3([p(1) p(1)+scale*B(1,3)], [p(2) p(2)+scale*B(2,3)], [p(3) p(3)+scale*B(3,3)], 'b', 'LineWidth', 2);

            end
            
            est_mix{kk}.plot_distributions([1,2,3],'Color','r','LineWidth',2); hold off
        end
    
        scale_for_plots = 0.9;

        xlim(scale_for_plots*[0 30])
        ylim(scale_for_plots*[0 30])
        zlim(scale_for_plots*[-15 15]) 

        subplot(1,2,2)
        cla 
        for kk = 1:length(sensors)
            % points{kk,k} = sensors{kk}.get_points(TR_active,T);
            curr_points = trimmed_meas{kk};
            for ii = 1:length(curr_points)
                scatter3(curr_points{ii}(1,:),curr_points{ii}(2,:),curr_points{ii}(3,:),6,'*'); hold on
            end
            % sensors{kk}.plot()
        end
    
        xlim(scale_for_plots*[0 30])
        ylim(scale_for_plots*[0 30])
        zlim(scale_for_plots*[-15 15]) 

        drawnow;
        
        frame = getframe(f);
        img   = frame2im(frame);
        
        % convert to an indexed image with a fixed 256‐color map
        [A, map] = rgb2ind(img, 256);
        
        % write to GIF: first frame creates the file, subsequent frames append
        if k == 1
            imwrite(A, map, gifFilename, 'gif', ...
                    'LoopCount', Inf, ...      % make it loop forever
                    'DelayTime', 0.1);         % seconds between frames
        else
            imwrite(A, map, gifFilename, 'gif', ...
                    'WriteMode', 'append', ...
                    'DelayTime', 0.1);
        end

        saveas(gcf,[pngDir, sprintf('frame_%i',k),'.png'])
    
    end
end

%%

% figure;
% for kk = 1:length(sensors) 
%     scatter3(points{kk,1}(:,1),points{kk,1}(:,2),points{kk,1}(:,3),6,'*'); hold on
%     scatter3(points{kk,50}(:,1),points{kk,50}(:,2),points{kk,50}(:,3),6,'*'); hold on
%     sensors{kk}.plot()
% end

%% Functions

function path = makeCirclePath(c, r, w)
    path.p = @(t) [c(1)+r*cos(w*t); c(2)+r*sin(w*t); c(3)];
    path.v = @(t) [-r*w*sin(w*t);   r*w*cos(w*t);    0];
    path.a = @(t) [-r*w^2*cos(w*t); -r*w^2*sin(w*t); 0];
end

function path = makeCircleWithZPath(c, r, w, wz)
    path.p = @(t) [c(1)+r*cos(w*t); c(2)+r*sin(w*t); c(3)+0.8*sin(wz*t)];
    path.v = @(t) [-r*w*sin(w*t);   r*w*cos(w*t);    0.8*wz*cos(wz*t)];
    path.a = @(t) [-r*w^2*cos(w*t); -r*w^2*sin(w*t); -0.8*wz^2*sin(wz*t)];
end

function path = makeEllipsePath(c, ab, w)
    ax = ab(1); by = ab(2);
    path.p = @(t) [c(1)+ax*cos(w*t); c(2)+by*sin(w*t); c(3)];
    path.v = @(t) [-ax*w*sin(w*t);   by*w*cos(w*t);    0];
    path.a = @(t) [-ax*w^2*cos(w*t); -by*w^2*sin(w*t); 0];
end

function path = makeLemniscatePath(c, a, w)
    % Bernoulli lemniscate (figure-eight)
    path.p = @(t) c(:) + [a*cos(w*t)./(1+sin(w*t).^2); ...
                          a*(sin(w*t).*cos(w*t))./(1+sin(w*t).^2); 0];
    epsd = 1e-3;
    path.v = @(t) (path.p(t+epsd)-path.p(t-epsd))/(2*epsd);
    path.a = @(t) (path.v(t+epsd)-path.v(t-epsd))/(2*epsd);
end

function path = makeLissajousPath(c, ab, mn, w)
    m = mn(1); n = mn(2);
    ax = ab(1); by = ab(2);
    path.p = @(t) [c(1)+ax*sin(m*w*t+pi/6); c(2)+by*sin(n*w*t); c(3)];
    path.v = @(t) [ ax*m*w*cos(m*w*t+pi/6);  by*n*w*cos(n*w*t); 0];
    path.a = @(t) [-ax*(m*w)^2*sin(m*w*t+pi/6); -by*(n*w)^2*sin(n*w*t); 0];
end

function path = makeSineWeavePath(c, drift, amp, w)
    phi = pi/5; tf_guess = 10; k = 1/tf_guess;
    path.p = @(t) [c(1)+drift(1)*k*t + amp(1)*sin(w*t); ...
                   c(2)+drift(2)*k*t + amp(2)*sin(0.7*w*t+phi); ...
                   c(3)+drift(3)*k*t];
    path.v = @(t) [drift(1)*k + amp(1)*w*cos(w*t); ...
                   drift(2)*k + amp(2)*0.7*w*cos(0.7*w*t+phi); ...
                   drift(3)*k];
    path.a = @(t) [-amp(1)*w^2*sin(w*t); ...
                   -amp(2)*(0.7*w)^2*sin(0.7*w*t+phi); ...
                    0];
end

function path = makeRoundedRectPath(minCorner, maxCorner, r, w)
    % Rounded rectangle in XY plane, Z fixed.
    cx = [minCorner(1)+r, maxCorner(1)-r, maxCorner(1)-r, minCorner(1)+r];
    cy = [maxCorner(2)-r, maxCorner(2)-r, minCorner(2)+r, minCorner(2)+r];
    z0 = minCorner(3);
    per = 2*((maxCorner(1)-minCorner(1)-2*r) + (maxCorner(2)-minCorner(2)-2*r)) + 2*pi*r;
    wrap = @(u) mod(u, per);

    function P = xy(u)
        Ls = [(maxCorner(1)-minCorner(1)-2*r), pi*r/2, ...
              (maxCorner(2)-minCorner(2)-2*r), pi*r/2, ...
              (maxCorner(1)-minCorner(1)-2*r), pi*r/2, ...
              (maxCorner(2)-minCorner(2)-2*r), pi*r/2];
        cum = [0, cumsum(Ls)];
        uu = wrap(u);
        seg = find(uu >= cum(1:end-1) & uu < cum(2:end), 1, 'first');
        s = uu - cum(seg);
        switch seg
            case 1, P = [minCorner(1)+r + s, maxCorner(2)-r];
            case 2, th = s/r; P = [cx(2)+r*cos(-pi/2+th), cy(2)+r*sin(-pi/2+th)];
            case 3, P = [maxCorner(1)-r, maxCorner(2)-r - s];
            case 4, th = s/r; P = [cx(3)+r*cos(0+th),    cy(3)+r*sin(0+th)];
            case 5, P = [maxCorner(1)-r - s, minCorner(2)+r];
            case 6, th = s/r; P = [cx(4)+r*cos(pi/2+th), cy(4)+r*sin(pi/2+th)];
            case 7, P = [minCorner(1)+r, minCorner(2)+r + s];
            case 8, th = s/r; P = [cx(1)+r*cos(pi+th),   cy(1)+r*sin(pi+th)];
        end
    end

    pXY = @(t) xy(w*t);
    path.p = @(t) local_p(pXY, t, z0);
    epsd = 1e-3;
    path.v = @(t) (path.p(t+epsd)-path.p(t-epsd))/(2*epsd);
    path.a = @(t) (path.v(t+epsd)-path.v(t-epsd))/(2*epsd);

    function P3 = local_p(pXY, t, z)
        xyv = pXY(t);
        P3 = [xyv(1); xyv(2); z];
    end
end

function path = makeWanderPath(c, amp, w)
    path.p = @(t) [c(1)+amp(1)*sin(w*t); ...
                   c(2)+amp(2)*cos(0.9*w*t+0.7); ...
                   c(3)+amp(3)*sin(1.2*w*t)];
    path.v = @(t) [ amp(1)*w*cos(w*t); ...
                   -amp(2)*0.9*w*sin(0.9*w*t+0.7); ...
                    amp(3)*1.2*w*cos(1.2*w*t)];
    path.a = @(t) [-amp(1)*w^2*sin(w*t); ...
                   -amp(2)*(0.9*w)^2*cos(0.9*w*t+0.7); ...
                   -amp(3)*(1.2*w)^2*sin(1.2*w*t)];
end

function [p,v,a] = evalPVA(path, t)
    p = path.p(t);
    v = path.v(t);
    a = path.a(t);
end