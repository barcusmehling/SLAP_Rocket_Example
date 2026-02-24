% Run script to select rocket force and accelerometer locations from FEM
% Clicking location on plots saves to accel_nodes and shaker_nodes (forces)
% variables
% accel_nodes - nodes x 2 (local node number, fem node number)
% shaker_nodes - nodes x 5 (local node number, fem node number,xdof,ydof,zdof)
%% ===================== LOAD DATA =====================
clc;close all;clear all;
disp('Loading Rocket Modes...')
load ..\ModeShapes\FullRocketModes;   % loads nodes, phi, etc.
clc;
% ===================== BUILD LOCATIONS for rocket =====================
ths = 0:30:360-30;

r = [1.73 3.05*ones(1,length(6.5:13:110.5)) 3.05-3.05/4:-3.05/4:0];
zs = [0.175 6.5:13:110.5 117:6.5:117+3*6.5];

accel_locs = [];
c = 1;
for ii = 1:length(r)
    for jj = 1:length(ths)
        accel_locs(c,2:4) = [r(ii)*cosd(ths(jj)) r(ii)*sind(ths(jj)) zs(ii)];
        c = c + 1;
    end
end

accel_locs(:,2:4) = accel_locs(:,2:4) / 39.37;
accel_locs(:,1) = 1:size(accel_locs,1);
accel_locs = accel_locs(1:157,:);

% fins
r = [3.5 5.5 7.5];
zs = [9 11.5];

ths = [0 90 180 270];

accel_locs2 = [];
c = 1;
for ii = 1:length(r)
    for jj = 1:length(ths)
        for kk = 1:length(zs)
            accel_locs2(c,2:4) = [r(ii)*cosd(ths(jj)) r(ii)*sind(ths(jj)) zs(kk)]/39.37;
            c = c + 1;
        end
    end
end

accel_locs = [accel_locs;accel_locs2];
accel_locs(:,1) = 1:size(accel_locs,1);

% ===================== MAP TO FE NODES =====================
dmin = zeros(size(accel_locs,1),1);
accel_nodes = dmin;

for k = 1:size(accel_locs,1)
    d = sum((nodes(:,2:4) - accel_locs(k,2:4)).^2,2);
    [dmin(k), accel_nodes(k)] = min(d);
end

disp('Maximum distance between accel location and FE node:')
disp([num2str(max(sqrt(dmin))) ' m'])

rocket_nodes = nodes(accel_nodes,:);

% %%%%%%%%%%%%%%%% DUT and baseplate nodes
% DUT nodes
clc;close all;
barc_bpnodes = 1:202365;
basedut_FEMnodes = nodes(barc_bpnodes,:); % base and DUT FEM nodes in Rocket

ths = 0:30:360-30; % set array of angles for baseplate circle
r = 2.8; % close to outer radius of baseplate
z = 90.1678; % height of top surface of baseplate

base_accel_locs = []; % initialize empty matrix for accel locs on baseplate
for jj = 1:length(ths)
    base_accel_locs(jj,2:4) = [r*cosd(ths(jj)) r*sind(ths(jj)) z];
end

% assemble a few more accels on baseplate
vec1 = [-1.2 0 1.2];
base_accel_locs2 = [];
c = 1;
for ii = 1:length(vec1)
    for jj = 1:length(vec1)
        base_accel_locs2(c,2:4) = [vec1(ii) vec1(jj) z];
        c = c + 1;
    end
end

base_accel_locs = [base_accel_locs;base_accel_locs2]; % combine into single matrix
base_accel_locs(:,1) = 1:size(base_accel_locs,1); % number accel locs starting from one

% BARC
xs = [-0.4 -0.2 0 0.2 0.4];
ys1 = [-2.5 2.5];
ys2 = linspace(-2.5,2.5,7);
zs = linspace(90.75,92.62,5);

barc_accel_locs = [];
c = 1;
for ii = 1:length(xs)
    for kk = 1:length(zs)
        if zs(kk) == 92.62
            for jj = 1:length(ys2)
                barc_accel_locs(c,:) = [c xs(ii) ys2(jj) zs(kk)];
                c = c + 1;
            end
        else
            for jj = 2
                barc_accel_locs(c,:) = [c xs(ii) ys1(jj) zs(kk)];
                c = c + 1;
            end
        end
        
    end
end

ys3 = linspace(-2.5,-1.5,7);
x = -0.5;

for ii = 1:length(ys3)
    for jj = 1:length(zs)
        barc_accel_locs(c,:) = [c x ys3(ii) zs(jj)];
        c = c + 1;
    end
end

barc_base_accel_locs = [barc_accel_locs;base_accel_locs]; % assemble into one matrix
barc_base_accel_locs(:,2:4) = barc_base_accel_locs(:,2:4) / 39.37;

% find FEM nodes
dmin = zeros(size(barc_base_accel_locs,1),1); % next few lines get FEM locations corresponding to those in base_dut_accel_locs
barc_base_FEM_nodes = dmin; % initialize FEM nodes of baseplate and DUT to zeros

for k = 1:size(barc_base_accel_locs,1) % inc through locs
    d = sum((basedut_FEMnodes(:,2:4) - barc_base_accel_locs(k,2:4)).^2,2);
    [dmin(k), barc_base_FEM_nodes(k)] = min(d);
end

disp('Maximum distance between accel location and FE node:')
disp(max(sqrt(dmin)))

rocket_nodes2 = basedut_FEMnodes(barc_base_FEM_nodes,:);

rocket_nodes3 = [rocket_nodes;rocket_nodes2];

% ===================== SHARED STATE =====================
shaker_nodes = [];   % [order nodeID ux uy uz]
accel_nodes  = [];   % [order nodeID]
shaker_count = 0;
accel_count  = 0;

% ===================== SHAKER SELECTION =====================
figure('Position',[500 200 350 550])
h1 = scatter3(rocket_nodes(:,2),rocket_nodes(:,3),rocket_nodes(:,4),25,'k','filled');
xlabel('x'); ylabel('y'); zlabel('z')
grid on; axis equal;
title('Select flight force nodes')
hold on
set(h1,'ButtonDownFcn',@shakerCallback)

% ===================== ACCEL SELECTION =====================
figure('Position',[900 200 350 550])
h2 = scatter3(rocket_nodes2(:,2),rocket_nodes2(:,3),rocket_nodes2(:,4),25,'k','filled');
xlabel('x'); ylabel('y'); zlabel('z')
grid on; axis equal
title('Select accel nodes')
hold on
set(h2,'ButtonDownFcn',@accelCallback)

% save
% save('FlightAccelLocs','accel_locs')
% save('FlightForceLocs','shaker_locs')

%% ===================== CALLBACKS =====================
function shakerCallback(~,event)

    rocket_nodes = evalin('base','rocket_nodes');
    shaker_nodes = evalin('base','shaker_nodes');
    shaker_count = evalin('base','shaker_count');

    Pclick = event.IntersectionPoint;
    [~,idx] = min(sum((rocket_nodes(:,2:4)-Pclick).^2,2));

    P = rocket_nodes(idx,2:4);
    nodeID = rocket_nodes(idx,1);

    answ = inputdlg({'Ux','Uy','Uz'},'Shaker Direction',1,{'0','0','1'});
    if isempty(answ); return; end

    u = str2double(answ{1});
    v = str2double(answ{2});
    w = str2double(answ{3});

    mag = norm([u v w]);
    if mag == 0; return; end

    u = u/mag;
    v = v/mag;
    w = w/mag;

    shaker_count = shaker_count + 1;

    quiver3(P(1),P(2),P(3),0.25*u,0.25*v,0.25*w,...
            'b','LineWidth',2,'MaxHeadSize',2);
    scatter3(P(1),P(2),P(3),60,'r','filled')
    text(P(1),P(2),P(3),sprintf('  %d',shaker_count), ...
         'Color','r','FontWeight','bold')

    shaker_nodes(end+1,:) = [shaker_count nodeID u v w];

    assignin('base','shaker_nodes',shaker_nodes);
    assignin('base','shaker_count',shaker_count);
end

function accelCallback(~,event)
    rocket_nodes2 = evalin('base','rocket_nodes2');
    accel_nodes   = evalin('base','accel_nodes');
    accel_count   = evalin('base','accel_count');

    Pclick = event.IntersectionPoint;
    [~,idx] = min(sum((rocket_nodes2(:,2:4)-Pclick).^2,2));

    P = rocket_nodes2(idx,2:4);
    nodeID = rocket_nodes2(idx,1);

    accel_count = accel_count + 1;

    scatter3(P(1),P(2),P(3),60,'g','filled')
    text(P(1),P(2),P(3),sprintf('  %d',accel_count),'Color','g','FontWeight','bold')

    accel_nodes(end+1,:) = [accel_count nodeID];
    assignin('base','accel_nodes',accel_nodes);
    assignin('base','accel_count',accel_count);
end









