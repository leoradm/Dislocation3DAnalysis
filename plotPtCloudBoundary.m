function plotPtCloudBoundary(Boundary,solved_boundary,bNum,window_position,color)
% this code plots and formats the Boundary with its plane and normal vector

% define inputs
vector_length = 15; % length of normal vector in plot

% plot boundary
figure
pcshow(Boundary(bNum).ptCloud) % plot the pointcloud
hold on

% define & plot normal vector in center of cloud
refPT = mean(Boundary(bNum).XYZ_Points); 
normalVEC = solved_boundary * vector_length; 
quiver3(refPT(1),refPT(2),refPT(3),normalVEC(1),normalVEC(2),normalVEC(3),'LineWidth',3);
plotNormalPlane(refPT,normalVEC',color(bNum,:))
hold off

% format the plot nicely
set(gca,'XColor','k','YColor','k','ZColor','k','Color','w','GridColor','k','fontsize',26,'linewidth',2)
set(gcf,'Color','w','Position', window_position)
a=gca;
a.XLabel.Interpreter = 'latex';
a.YLabel.Interpreter = 'latex';
a.ZLabel.Interpreter = 'latex';
box on

end