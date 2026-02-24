% Same as GetRocketLocationsScript.m but for lab config
clc;close all;clear all;
disp('Loading Lab Modes...')
load ..\ModeShapes\BARC_BaseplateModes;   % loads nodes, phi, etc.
clc;
% DUT and baseplate nodes
% DUT nodes

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

figure;
scatter3(barc_base_accel_locs(:,2),barc_base_accel_locs(:,3),barc_base_accel_locs(:,4),'b','filled')

% find FEM nodes
dmin = zeros(size(barc_base_accel_locs,1),1); % next few lines get FEM locations corresponding to those in base_dut_accel_locs
barc_base_FEM_nodes = dmin; % initialize FEM nodes of baseplate and DUT to zeros

for k = 1:size(barc_base_accel_locs,1) % inc through locs
    d = sum((nodes(:,2:4) - barc_base_accel_locs(k,2:4)).^2,2);
    [dmin(k), barc_base_FEM_nodes(k)] = min(d);
end

disp('Maximum distance between accel location and FE node:')
disp(max(sqrt(dmin)))

FEMnodes = nodes(barc_base_FEM_nodes,:);

figure;
scatter3(FEMnodes(:,2),FEMnodes(:,3),FEMnodes(:,4),'r','filled') % plot dut and base accels

% ===================== SHARED STATE =====================
shaker_nodes = [];   % [order nodeID ux uy uz]
accel_nodes  = [];   % [order nodeID]
shaker_count = 0;
accel_count  = 0;

% ===================== SHAKER SELECTION =====================
figure('Position',[500 200 350 550])
h1 = scatter3(FEMnodes(:,2),FEMnodes(:,3),FEMnodes(:,4),25,'k','filled');
xlabel('x'); ylabel('y'); zlabel('z')
grid on; axis equal;
title('Select flight force nodes')
hold on
set(h1,'ButtonDownFcn',@shakerCallback)

% ===================== ACCEL SELECTION =====================
figure('Position',[900 200 350 550])
h2 = scatter3(FEMnodes(:,2),FEMnodes(:,3),FEMnodes(:,4),25,'k','filled');
xlabel('x'); ylabel('y'); zlabel('z')
grid on; axis equal
title('Select accel nodes')
hold on
set(h2,'ButtonDownFcn',@accelCallback)

% save
% save('LabAccelLocs','accel_locs')
% save('LabShakerLocs','shaker_locs')

%% ===================== CALLBACKS =====================
function shakerCallback(~,event)

    rocket_nodes = evalin('base','FEMnodes');
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

    u = 0.1*u/mag;
    v = 0.1*v/mag;
    w = 0.1*w/mag;

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
    rocket_nodes2 = evalin('base','FEMnodes');
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









