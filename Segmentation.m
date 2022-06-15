% Leora Dresselhaus-Marais
% Jan 27, 2020
%
% This code will read in scan data and compile a data cube from a z-scan,
% then project that scan onto another plane.

close all

%% USER INPUTS
% file specifiers
writePath = 'G:\My Drive\BoxMigration\Pin-Hua Files\Updated file\'; % define filepath for saved output files
dataPath = 'D:\blc12704\id06\S1_590A\S1_590A_Al_200_ff_saturday\'; % define the file path for your EDF data (include terminal slash)
% NOTE: the 200 scans were named incorrectly. This IS the correct scan, though!!!
scanName = 'Al_111_rockinglayer_';%_-75um_smy'; % specify the name of this run (should be a subfolder inside the specified 'dataPath' folder)
layers = dir([dataPath,scanName,'*']);
wbNum1 = '0014'; %which WB rocking layers are we packing into this array?
wbNum2 = '0004'; %which WB rocking layers are we packing into this array?

% experiment parameters
M = 17;%x, (magnification)
thresh = 110;% arbitrary thresholding value to find intensity corresponding to a dislocation (CHANGE ON EACH RUN)

% booleans to save/plot different things (for all Q's below, 1=Y; 0=N)
plotRawImON = 0; % Do you want to plot the raw images from the scan?
videoRawImON = 0; % Do you want to save the Raw Images you plotted to a movie?
save3Dcube = 1; % Should we save the cube of 3D data from the thresholded image stack?
save3Dfig = 0; % Should we save a FIG file for the 3D visualization?

% define diffraction peak:
PK = 'Pk200';

% Specify filename
VERSION = '1'; % what version number is this in the filenames?
FILENAME = [PK,'_WB_Segmented_'];

%% START BY EXPLORING WHICH SCAN NUMBER TO USE FOR WB
%{
% Run this section of code when initially starting analysis for a new
% experiment. This section plots all images in an angular scan so that you
% can manually identify the step that is characteristic of the weak-beam
% contrast for subsequent analysis.

figure
scNames = dir([dataPath,scanName,'00\*.edf']);
for sNum=1:length(scNames)
    im = edf_read([dataPath,scanName,'00\',scNames(sNum).name]);
    imagesc(im);
    title(['WB Scan ',num2str(sNum-1)])
    colorbar
    colormap gray
    caxis([85,500])
    set(gca,'FontSize',36)
    pause
end
%}
%% START WITH RESOLUTION TARGET...
%{
% When we start a new run, we begin by calibrating the magnification and
% associate number of um at the sample that describe each imaged pixel at our
% detector. To do this, we load the images called "alignment_cdx" which saves
% a stack of images collected on the same object, translated by a fixed known
% amount.

ResPath = ['D:\blc12704\id06\alignment\alignment_cdx\snaps\'];
files = dir([ResPath,'*.edf']);

% read resolution files
for fNum=4%1:length(files)
    im = edf_read([ResPath,files(fNum).name]);
    imagesc(im)
    title(strrep(files(fNum).name(1:10),'_',' '))
end

xlim([900,1075])
ylim([1025,1275])
set(gca,'ydir','normal','xdir','reverse','YAxisLocation','right','FontSize',36)
xlabel('x')
ylabel('y')
imdistline
%}
% so it is 134 pixels per 200-nm
UMperPX = 0.200/4.5; 
M = 6.500/UMperPX;
buffer = 250; % number of pixels to clip on the edges of the image to avoid edge artifacts


%% BEGIN RUNNING CODE:
% set slash direction based on operating system (will do this automatically):
if ispc==1
    slash = '\';
else
    slash = '/';
end

% load filenames, taking into account that WB frame changes over the scan
fNames = cell(1,length(layers));
for m=1:length(layers)
    % find number explaining layer
    index = strfind(layers(m).name,'_');
    NUM = layers(m).name(index(end)+1:end);
    num = str2num(NUM);
    
    % updated WB number manually to accommodate shifts over the scan
    if num <= 180
        wbNum_current = wbNum1;
    elseif num > 180 && num <= 230
        wbNum_current = '0014';
    elseif num > 230 && num <= 301
        wbNum_current = '0015';
    end
    wb{num+1} = wbNum_current;
    
    % save appropriate WB filename for this layer
    fNames{num+1} = [dataPath,layers(m).name,slash,'Al_111_rockinglayer_',NUM,'_',wbNum_current,'.edf'];
end

% load the information about scan from metadata
info1 = edf_info(fNames{1}); % load metadata
Dz = info1.motor.ffz; % position of the detector
Dx = info1.motor.mainx*-1; % position of the detector
TwoTheta = atand(Dz/Dx); % relevant 2theta for this story
y = [1:1:info1.dim_2]*UMperPX; % compute x-axis in image
x = [1:1:info1.dim_1]*UMperPX/sind(TwoTheta); % compute y-axis in image
z = [1:1:length(fNames)]*cosd(abs(info1.motor.diffry)); % vertical scanning axis

% remove buffered area (edge artifacts from detector)
x = x(buffer:end-buffer+1);
y = y(buffer:end-buffer+1);

% instantiate plotting parameters to load and plot files
if plotRawImON == 1
    figure('Position', [-744.3333 761 1536 1160]);%[793 761 1.5347e+03 1160]);
end
scale = 1;%/4; % scaling of the image size to reduce data load (i.e. binning step)
ind = 1; % counter for dislocation points
SE = strel('disk',50); % structuring element to identify dislocations
x = imresize(x,scale); % bin axis
y = imresize(y,scale); % bin axis

% open image if saving
if videoRawImON==1
    % create movie of the raw frames
    fid = VideoWriter([writePath,FILENAME,'OverlayImages_v',VERSION,'.avi']);
    fid.FrameRate = 5;
    open(fid)
end

% define noise floor for all three scans:
noiseFloor = 110;
refPosition = 50*-1;

% begin to scan through the filenames to load images
fSt = 1; % what number to start the fNum loop?
z = [1:1:length(fNames)]*cosd(abs(info1.motor.diffry)); % vertical scanning axis
filteringMethod = 2; % which type of filtering method should we use?
for fNum=fSt:length(fNames)
    % load file
    im = edf_read(fNames{fNum}); % read image
    im = im(buffer:end-buffer+1,buffer:end-buffer+1); % remove extraneous pixels at edges
    im = transpose(im); % flip

    [szX,szY] = size(im);
    
    % remove the edge pixes from each image:
    im2 = im;
    im2(:,1:3) = ones(szX,3)*noiseFloor; % remove noise at edges of image
    im2(:,end-2:end) = ones(szX,3)*noiseFloor; % remove noise at edges of image
    im2(1:3,:) = ones(3,szY)*noiseFloor; % remove noise at edges of image
    im2(end-2:end,:) = ones(3,szY)*noiseFloor; % remove noise at edges of image
    
    % remove the noise below the noise-floor:
    im2(im2<=noiseFloor) = noiseFloor;
    
    % save the raw file
    if fNum==fSt
        [len,wid] = size(im);
        CUBE_filt = zeros(len,wid,length(z));
    end
	CUBE_filt(:,:,fNum) = im2;
        
    if plotRawImON == 1
        % plot image
        i1 = imagesc(y,x,log(im));
        title(['Obj 10X Images (frame ',num2str(fNum),'), in WB step ',wb{fNum}])
        caxis([log(noiseFloor*2), log(max(max(im)))*0.9])
        colorbar
        colormap jet
        xlabel('y ($\mu$m)','interpreter','')
        ylabel('y ($\mu$m)','interpreter','')
        set(gca,'ydir','normal','xdir','reverse','YAxisLocation','right','FontSize',25)
        hold off
        
        if videoRawImON==1
            writeVideo(fid,getframe(gcf));
        else
            pause(0.1)
        end
    end

    fprintf(['Image ',num2str(fNum),' of ',num2str(length(fNames)),'.\n'])
end
if videoRawImON==1
    close(fid)
end
clear fid im2 im wb fNames

%% find the raw edges
% rescale the image for the fiber filtering method to work
factor = 0.6; % scaling factor to prevent Matlab from running out of memory
CUBE_filtScaled = imresize(CUBE_filt,factor); % scale image cube
clear CUBE_filt

% define binned x,y,z grid
[sX,sY,sZ] = size(CUBE_filtScaled);
dX = (x(end)-x(1))/sX;    dY = (y(end)-y(1))/sY;    dZ = (z(end)-z(1))/sZ;
[Y,X,Z] = meshgrid([1:sY]*dY,[1:sX]*dX,[1:sZ]*dZ);
xSc = linspace(dX,max(x),sX);
ySc = linspace(dY,max(y),sY);
zSc = linspace(dZ,max(z),sZ);

% use the fiber enhancement method to connect dislocation points
fprintf('Finding Dislocations with fiber filter method... ')
SSens = 5*diff(getrangefromclass(CUBE_filtScaled)); % Compute structure senstivity
CUBE = fibermetric(CUBE_filtScaled,'StructureSensitivity',SSens); % enhance
fprintf('DONE!\n')

%% define color limits for this system:
clim = [0.1,1]; % define color limits

% now strip the values below clim(1):
fprintf('Thresholding fiber-textured images...')
XYZ_Points = zeros(3,1);
CUBEn = CUBE;
counter = 1; % 
for m=1:sX
    for n=1:sY
        for o=1:sZ
            if CUBE(m,n,o) < clim(1)
                CUBEn(m,n,o) = NaN;
            else
                XYZ_Points(:,counter) = [xSc(m),ySc(n),zSc(o)];
                counter = counter+1;
            end
        end
    end
end
clear CUBE
fprintf('DONE!\n')

% save point list as input for "BoundaryAnalysis.m"
fprintf('Saving points for pointCloud operations...')
save([writePath,'DislocationXYZ_Points.mat'],'XYZ_Points','-v7.3') 
fprintf('DONE!\n')

% save the cube of data (if necessary)
if save3Dcube ==1
    fprintf('Now saving the full 3D Cube... ')
    save([writePath,FILENAME,'3Dcube_FiberGradFilt_v',VERSION,'.mat'],'CUBEn','xSc','ySc','zSc','clim','TwoTheta','-v7.3')
    fprintf('DONE!\n')
end

