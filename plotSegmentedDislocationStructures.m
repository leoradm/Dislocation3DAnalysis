% Leora Dresselhaus-Marais
% June 14, 2022
%
% This code will take the outputs saved from the Segmentation.m script, and
% plot the 3D dislcoation structures for the associated runs.

close all
clear all

%% USER INPUTS
% file specifiers
writePath = 'G:\My Drive\BoxMigration\Pin-Hua Files\Updated file\';
dataPath = writePath;
file = 'Pk200_WB_Segmented_3Dcube_FiberGradFilt_v1.mat';
tic
fprintf('Loading the file...')
load([writePath,file]); % load the files
fprintf('DONE!\n')
toc
clear X Y Z

% experiment parameters
M = 17;%x, (magnification)
viewingVector = [83.2297 43.7761]; % vector defining the view for 2D slices (final section)

% booleans to save/plot different things (for all Q's below, 1=Y; 0=N)
save3Dfig = 0; % Should we save a FIG file for the 3D visualization?
Spin3D = 1; % Should we save a video for the 3D spin movie?

% define diffraction peak:
PK = 'Pk200';

% Specify filename
VERSION = '1'; % what version number is this in the filenames?
FILENAME = [PK,'_WB_Segmented_'];


%% Plot the identified dislocation structure
% plot the figure using the PATCH_3Darray() function
fprintf('Plotting the arrays...')
figure('Position',[-744.3333 761 1536 1160])
clim = [0.1,1]; % define color limits
cmap = bone(64); % define colormap
PATCH_3Darray(CUBEn,xSc,ySc,zSc,clim,cmap,'col')
fprintf('DONE!\m')

Plot3D_CITATION = ['Adam A (2021). Plot a 3D array using patch '...
'(https://www.mathworks.com/matlabcentral/fileexchange/28497-plot-a-3d-'...
'array-using-patch), MATLAB Central File Exchange. Retrieved August 5, 2021.'];

% format plot
axes = gca;
fig = gcf;
set(axes,'FontSize',46,'FontName','Times New Roman','XGrid','on','YGrid','on','ZGrid','on','LineWidth',2,'xdir','normal','ydir','normal','zdir','normal')
set(fig,'Color',[1,1,1]);
xlabel('$\textit{x} (\mu m)$','interpreter','latex')
ylabel('$\textit{y} (\mu m)$','interpreter','latex')
zlabel('$\textit{z} (\mu m)$','interpreter','latex')

% Make lighting and transparency appropriate for all patches
for pNum=1:length(axes.Children)
    axes.Children(pNum).FaceAlpha = 0.2;
    axes.Children(pNum).EdgeAlpha = 0.2;
end
camlight left
lighting gouraud
box on

% plot an orange plane to show how the images slice through the structure
hold on
orange = [1, 0.4118, 0.1608];
RANGE = [0,max(xSc);0,max(ySc);0,max(zSc)]; % (range for images
xlim([-5,RANGE(1,2)+25]);     ylim([-25,RANGE(2,2)]);     zlim([0,RANGE(3,2)]);
zPos = 150; % position for the plane along the z-axis
vert = [RANGE(1,1),RANGE(2,1),zPos;...
        RANGE(1,2),RANGE(2,1),zPos;...
        RANGE(1,2),RANGE(2,2),zPos;...
        RANGE(1,1),RANGE(2,2),zPos]; %define vertices for patch
g = patch('Faces',[1,2,3,4],'Vertices',vert);
set(g,'FaceColor',orange,'FaceAlpha',0.5,'EdgeColor','none')

% draw crystallographic legend for the plot
CenterPos = [280,10,285]; % center of the axes
arrowLen = 15; % length of the arrow
fz = 36; % font size for crystallographic axis labels
ar1 = quiver3(CenterPos(1),CenterPos(2),CenterPos(3),...
    -arrowLen*sind(TwoTheta/2),0,arrowLen*cosd(TwoTheta/2),...
    'Color',[1,1,1]*0.5,'LineWidth',3,'MaxHeadSize',5);
text(CenterPos(1)-arrowLen*1.5,CenterPos(2),CenterPos(3)+arrowLen*1.5,'$\bf{[002]}$','FontSize',fz,'interpreter','latex','Color',[1,1,1]*0.5,'Rotation',TwoTheta/2);
ar2 = quiver3(CenterPos(1),CenterPos(2),CenterPos(3),...
    0,arrowLen,0,'Color','r','LineWidth',3,'MaxHeadSize',5);
text(CenterPos(1)-arrowLen,CenterPos(2)+arrowLen,CenterPos(3),'$\bf{[1\bar{1}0]}$','FontSize',fz,'interpreter','latex','Color','r','Rotation',TwoTheta/2);
ar3 = quiver3(CenterPos(1),CenterPos(2),CenterPos(3),...
    arrowLen*cosd(TwoTheta/2),0,arrowLen*sind(TwoTheta/2),...
    'Color',[0.4000 0.5098 0.2588],'LineWidth',3,'MaxHeadSize',5);
text(CenterPos(1)+arrowLen*3/4,CenterPos(2),CenterPos(3),'$\bf{[110]}$','FontSize',fz,'interpreter','latex','Color',[0.4000 0.5098 0.2588],'Rotation',TwoTheta/2);
daspect([1,1,1]) % set aspect ratio along each axis
axes.View = [-30,10]; % set the viewing vector

%% save the figure plot
if save3Dfig == 1
    fprintf('Now saving the full figure file... ')
    savefig([writePath,FILENAME,'Dislocations_',VERSION ,'.fig']);
    fprintf('DONE!\n')
end

%% Make the 3D spin movie of this map
if Spin3D==1
    tic
    fprintf('Now creating the movie...')
    rotationSteps = 300;
    fid = VideoWriter([writePath,FILENAME,'3DSpin_v',VERSION,'.mp4'],'MPEG-4');
    fid.FrameRate = 30;
    fid.Quality = 95;
    open(fid)
    
    % define initial viewing vector
    VIEW = axes.View;
    axes.View(2) = 11;
    
    % step through all possible rotations
    for rNum=1:rotationSteps
        VIEW(1) = VIEW(1) + 360/rotationSteps;
        axes.View = VIEW;
        writeVideo(fid,getframe(gcf));
        fprintf(['Saved frame ',num2str(rNum),' of ',num2str(rotationSteps),'.\n'])
    end
    close(fid)
    fprintf('DONE!\n')
    toc
end

