% 1/29/2026
clc;clear all;
load RocketModes_Full;
uimplot5({nodes,elestruct.eles,fn,phi})

dut_node_inds = 28554:39284;
base_node_inds = 1:22476;

dut_nodes = nodes(dut_node_inds,:);
baseplate_nodes = nodes(base_node_inds,:);

nodes_dutbase = [baseplate_nodes;dut_nodes];

dut_phi = phi(dut_node_inds(1)*6-5:dut_node_inds(end)*6,:);
base_phi = phi(base_node_inds(1)*6-5:base_node_inds(end)*6,:);

phi_dutbase = [base_phi;dut_phi];

uimplot5({nodes_dutbase,elestruct.eles(1:26945,:),fn,phi_dutbase})

% 50 rocket modes below 3,000 Hz

%% 2/9/2026
clc;close all;clear all;
load FullRocketModesv03;
uimplot5({nodes,elestruct.eles,fn,phi})
return

inner_nodes = 207294:208503;
inner_els = 165350:166407;

uimplot5({nodes(inner_nodes,:),elestruct.eles(inner_els,:),fn,phi(6*inner_nodes(1)-5:6*inner_nodes(end),:)})

barc_bpnodes = 1:202365;
barc_bpels = 1:160529; 

uimplot5({nodes(barc_bpnodes,:),elestruct.eles(barc_bpels,:),fn,phi(1:6*barc_bpnodes(end),:)})

barcnodes = 76918:202365;
barcels = 58740:160529;
uimplot5({nodes(barcnodes,:),elestruct.eles(barcels,:),fn,phi(6*barcnodes(1)-5:6*barcnodes(end),:)})

payload_nodes = 202366:202893;
drogue_nodes = 202894:204142;
bnodes = 206142:207293;

payload_els = 160530:161033;
drogue_els = 161034:162257;
b_els = 164222:165349;

pd_nodes = [payload_nodes drogue_nodes];
pd_els = [payload_els drogue_els];

nodes(bnodes,:) = [];
phi(bnodes(1)*6-5:bnodes(end)*6,:) = [];

ind1 = find(elestruct.eles(:,1)==b_els(1));
ind2 = find(elestruct.eles(:,1)==b_els(end));

elestruct.eles(ind1:ind2,:) = [];

nodes(pd_nodes,:) = [];
phi(pd_nodes(1)*6-5:pd_nodes(end)*6,:) = [];

ind1 = find(elestruct.eles(:,1)==pd_els(1));
ind2 = find(elestruct.eles(:,1)==pd_els(end));

elestruct.eles(ind1:ind2,:) = [];

% % % nodes(:,1) = 1:length(nodes(:,1));
% % % elestruct.eles(:,1) = 1:length(elestruct.eles(:,1));

uimplot5({nodes,elestruct.eles,fn,phi})




%% DUT baseplate config
clc;close all;clear all;
load BARC_BaseplateModes;

barcnodes = 1:125448;
barcels = 1:101790;
uimplot5({nodes(barcnodes,:),elestruct.eles(barcels,:),fn,phi(1:barcnodes(end)*6,:)})

%% Plot accelerometer locations
clc;close all;clear all;
load BARC_BaseplateModes;
load LabAccelNodes;
load LabNodes;

ctrl_inds = 1:9;
dut_inds = 10:size(accel_nodes,1);

figure;
scatter3(FEMnodes(:,2),FEMnodes(:,3),FEMnodes(:,4),'k','filled') % plot dut and base accels

hold on;

scatter3(nodes(accel_nodes(ctrl_inds,2),2),nodes(accel_nodes(ctrl_inds,2),3),nodes(accel_nodes(ctrl_inds,2),4),'b','filled')
scatter3(nodes(accel_nodes(dut_inds,2),2),nodes(accel_nodes(dut_inds,2),3),nodes(accel_nodes(dut_inds,2),4),'r','filled')

xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)')

title('Control (Blue) and Modal Filter (Red) Accelerometer Locations','interpreter','tex')










