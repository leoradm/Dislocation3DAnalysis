% Leora Dresselhaus-Marais
% Feb 19, 2022
%
% This code will create a pointcloud and use PCA to reduce the
% dimensionality of the dislocation boudnary plane into a vector.

format long

% define the files
fprintf('Loading files...')
load('DislocationXYZ_Points.mat');
fprintf('DONE!\n')

% load the bounds of the different boundaries
Boundary(1).xRange = [100,125];
Boundary(1).yRange = [20,100];
Boundary(1).zRange = [100, 149];

Boundary(2).xRange = [75,150];
Boundary(2).yRange = [0,100];
Boundary(2).zRange = [250,300];

Boundary(3).xRange = [5,28];
Boundary(3).yRange = [20,100];
Boundary(3).zRange = [180,300];

Boundary(4).xRange = [165,225];
Boundary(4).yRange = [0, 60];
Boundary(4).zRange = [0, 60];

Boundary(5).xRange = [140, 150];
Boundary(5).yRange = [4, 35];
Boundary(5).zRange = [0, 50];

% define plotting colors for the boundaries of interest
color = [1,0,0; 0,1,0; 0,1,1; 0,0,1; 0.5176,0.7490,0.9412; 1,0,1];

%% separate full point array into separate ones for each boundary

% initialize array for each boundary
B1_XYZ_Points = zeros(size(XYZ_Points))*NaN;
B2_XYZ_Points = zeros(size(XYZ_Points))*NaN;
B3_XYZ_Points = zeros(size(XYZ_Points))*NaN;
B4_XYZ_Points = zeros(size(XYZ_Points))*NaN;
B5_XYZ_Points = zeros(size(XYZ_Points))*NaN;

% save all points bounded in each region
fprintf('Sorting all points into boundary regions...')
tic
parfor ptNum = 1:length(XYZ_Points)
	% now sort points based on boundary locations
    if withinRange(Boundary,1,XYZ_Points(ptNum,:))==1
        B1_XYZ_Points(ptNum,:) = XYZ_Points(ptNum,:);
    elseif withinRange(Boundary,2,XYZ_Points(ptNum,:))==1
        B2_XYZ_Points(ptNum,:) = XYZ_Points(ptNum,:);
    elseif withinRange(Boundary,3,XYZ_Points(ptNum,:))==1
        B3_XYZ_Points(ptNum,:) = XYZ_Points(ptNum,:);
    elseif withinRange(Boundary,4,XYZ_Points(ptNum,:))==1
        B4_XYZ_Points(ptNum,:) = XYZ_Points(ptNum,:);
    elseif withinRange(Boundary,5,XYZ_Points(ptNum,:))==1
        B5_XYZ_Points(ptNum,:) = XYZ_Points(ptNum,:);
    end
end
fprintf('DONE!\n')
toc

% remove NaN points from the cells in the array
fprintf('Remove unnecessary cells in the array...')
B1_XYZ_Points = B1_XYZ_Points(~isnan(B1_XYZ_Points));
Boundary(1).XYZ_Points = reshape(B1_XYZ_Points,[length(B1_XYZ_Points)/3,3]);
B2_XYZ_Points = B2_XYZ_Points(~isnan(B2_XYZ_Points));
Boundary(2).XYZ_Points = reshape(B2_XYZ_Points,[length(B2_XYZ_Points)/3,3]);
B3_XYZ_Points = B3_XYZ_Points(~isnan(B3_XYZ_Points));
Boundary(3).XYZ_Points = reshape(B3_XYZ_Points,[length(B3_XYZ_Points)/3,3]);
B4_XYZ_Points = B4_XYZ_Points(~isnan(B4_XYZ_Points));
Boundary(4).XYZ_Points = reshape(B4_XYZ_Points,[length(B4_XYZ_Points)/3,3]);
B5_XYZ_Points = B5_XYZ_Points(~isnan(B5_XYZ_Points));
Boundary(5).XYZ_Points = reshape(B5_XYZ_Points,[length(B5_XYZ_Points)/3,3]);
fprintf('DONE!\n')

%% Define the coordinate transforms between the crystal and lab systems
TransformBasis = [1/sqrt(2), -1/sqrt(2),  0;
                  1/sqrt(2),  1/sqrt(2),  0;
                  0,          0,          1]; % transform for cut crystal orientation
ang = 10.38; % theta rotation for DFXM in the Bragg condition
RotationMatrix = [cosd(ang),  0, sind(ang);...
                  0,          1, 0         ;...
                 -sind(ang), 0, cosd(ang)];  % accounts for Bragg condition
TransformBasis = TransformBasis * RotationMatrix; % full transform from crystal to lab systems

%% Define pointClouds for each boundary to define its normal

for bNum=1:length(Boundary)
    % define this as a pointCloud object
    Boundary(bNum).ptCloud = pointCloud(Boundary(bNum).XYZ_Points); % format list into pointCloud object
    Boundary(bNum).normals = pcnormals(Boundary(bNum).ptCloud,10); % compute normals
    Boundary(bNum).ptCloud = pointCloud(Boundary(bNum).XYZ_Points,'Normal',Boundary(bNum).normals); % integrate into pointCloud

    % use the pointCloud to model the boundary
    maxDistance = 7; % distance offset (um) allowed to define plane/outliers
    [model,~,~,meanError] = pcfitplane(Boundary(bNum).ptCloud,maxDistance);
    Boundary(bNum).PlaneModel = model;
    Boundary(bNum).meanError = meanError;

    % define the crystallographic vector for this plane
    Boundary(bNum).ZoneAxis_Cart = Boundary(bNum).PlaneModel.Normal';
    Boundary(bNum).ZoneAxis_Cryst = TransformBasis*Boundary(bNum).ZoneAxis_Cart;
    Boundary(bNum).ZoneAxis_hkl = 1 ./ Boundary(bNum).ZoneAxis_Cryst;
    hkl = Boundary(bNum).ZoneAxis_hkl;
    hklTEXT = [num2str(hkl(1)),',',num2str(hkl(2)),',',num2str(hkl(3))];
    fprintf(['B',num2str(bNum),' has crystallographic vector ',hklTEXT,...
        ';     mean error = ',num2str(Boundary(bNum).meanError),'.\n'])

    % plot the results
    window_position = [1 41 1280 607.3333];
    solved_boundary = Boundary(bNum).PlaneModel.Normal;
    plotPtCloudBoundary(Boundary,solved_boundary,bNum,window_position,color)
    title({['Boundary ',num2str(bNum)],['$hkl = $',hklTEXT,', Error = ',...
        num2str(Boundary(bNum).meanError)]},'interpreter','latex','Color',[0,0,0],'fontsize',26)
end
fprintf(' \n')

%% round the solved hkl values to integers from MSAC output

% save the best answers we got above and associated errors
Boundary(1).Solved_hkl = [1;1;0];        Boundary(1).Solved_Error = 2.3006;
Boundary(2).Solved_hkl = [-2;0;1];       Boundary(2).Solved_Error = 1.9808;
Boundary(3).Solved_hkl = [1;1;-2];       Boundary(3).Solved_Error = 1.2041;
Boundary(4).Solved_hkl = [6;1;2];        Boundary(4).Solved_Error = 1.7356;
Boundary(5).Solved_hkl = [1;2;0];        Boundary(5).Solved_Error = 1.7089;

% redo assignments and hkl labels with the rounded values
for bNum=1:length(Boundary)
    % finish the assignments and coordinate transforms
    Boundary(bNum).Solved_uvw = 1./Boundary(bNum).Solved_hkl; 
    Boundary(bNum).Solved_uvw(Boundary(bNum).Solved_uvw==Inf) = 0;
    Boundary(bNum).Solved_xyz = transpose(TransformBasis) * Boundary(bNum).Solved_uvw;

    % print the results
    hklTEXT = [num2str(Boundary(bNum).Solved_hkl(1)),',',...
        num2str(Boundary(bNum).Solved_hkl(2)),',',...
        num2str(Boundary(bNum).Solved_hkl(3))];
    fprintf(['Testing if B',num2str(bNum),' has crystallographic vector ',hklTEXT,...
        ';     (which would imply a mean error = ',num2str(Boundary(bNum).Solved_Error),'.\n'])

    % plot the results
    window_position = [1 41 1280 607.3333];
    solved_boundary = Boundary(bNum).Solved_xyz;
    plotPtCloudBoundary(Boundary,solved_boundary,bNum,window_position,color)
    title({['Boundary ',num2str(bNum)],['$hkl = $',hklTEXT,', Error = ',...
        num2str(Boundary(bNum).Solved_Error)]},'interpreter','latex','Color',[0,0,0],'fontsize',26)
end

%% Manually reduce the 2D plane into 1D points for line vectors

for bNum=1:5
    % first define reference vectors for the dimensional reduction
    refPT = mean(Boundary(bNum).XYZ_Points); %Center point for the boundaries
    refVector = Boundary(bNum).Solved_xyz'; % reference normal
    temp = null(refVector); % create coordinate system
    v1=temp(:,1)'; % save in-plane X
    v2=temp(:,2)'; % save in-plane Y
    v3 = refVector; % save ref normal vector
    Boundary(bNum).View = [v1, refPT(2);v2,refPT(3);v3,refPT(1);0,0,0,1]; 
    
    % Reduce dimensionality of plane by projecting onto normal
    p2dArray = zeros(length(Boundary(bNum).XYZ_Points),4);
    for m=1:length(Boundary(bNum).XYZ_Points)
        p2dArray(m,:) = Boundary(bNum).View * [Boundary(bNum).XYZ_Points(m,:)';0];
    end
    Boundary(bNum).XY = p2dArray(:,1:2);
%    scatter(Boundary(bNum).XY(:,1),Boundary(bNum).XY(:,2))
    
    % find rotation for vertical lines
    angles = 0:1:180;%54;%
    meanBoundaryPoint = [mean(Boundary(bNum).XY(:,1)),mean(Boundary(bNum).XY(:,2))];
    inputPTS_x = Boundary(bNum).XY(:,1)-meanBoundaryPoint(1); % keep positions centered about boundary
	inputPTS_y = Boundary(bNum).XY(:,2)-meanBoundaryPoint(2); % keep positions centered about boundary
    for m = 1:length(angles)
        % let's see what rotation angle will give the most clearly defined frequencies
        ang = angles(m);
        tform = affine2d([cosd(ang),-sind(ang),0; sind(ang),cosd(ang),0; 0,0,1]); % define rotation matrix for points
        [trX,trY] = transformPointsForward(tform,inputPTS_x,inputPTS_y); % transform
%        Boundary(bNum).Transformed = [trX,trY]; % save
        
        % now let's represent and confirm the data from this rotation
        [nums,vals] = hist(trX,150); % reduce to 1D
        if m==1
            im_HIST=nums';
            im_FFT=fftshift(fft(nums))';
        else
            im_HIST(:,end+1)=nums';
            im_FFT(:,end+1) = fftshift(fft(nums))';
        end
%        plot(vals,nums+ang*1e3); hold on % plot 1D traces
%        plot(linspace(1/(vals(end)-vals(1)), 1/(vals(2)-vals(1)),length(vals)),fftshift(fft(nums))+ang*1e5); hold on % plot FFT for frequencies
%    	 scatter(trX,trY,2); xlim([-1,1]); ylim([-1,1]); title(['rotated ',num2str(ang)]); pause(0.1) % plot movie of the rotation
    end
    
    % define best angle
    if bNum==1
        Boundary(bNum).LineAngle = -(90-52); %degrees
    elseif bNum==2
        Boundary(bNum).LineAngle = -60; %degrees
    elseif bNum==3
        Boundary(bNum).LineAngle = 86; %degrees
    elseif bNum==4
        Boundary(bNum).LineAngle = -20; %degrees
    elseif bNum==5
        Boundary(bNum).LineAngle = -34; %degrees
    else
        Boundary(bNum).LineAngle = 0;
    end
    
    % convert vector back to XYZ
    DefiningPoints = [meanBoundaryPoint; tand(Boundary(bNum).LineAngle)+meanBoundaryPoint(1),1+meanBoundaryPoint(2)]; % define as a vector
    if max(max(DefiningPoints))==Inf % in case the angle is 90
        DefiningPoints = [meanBoundaryPoint; 1+meanBoundaryPoint(1),1/tand(Boundary(bNum).LineAngle)+meanBoundaryPoint(2)]; % define as a vector
    end    
    line_vector_bp = DefiningPoints(2,:) - DefiningPoints(1,:);
    line_vector_bp(line_vector_bp==Inf) = 1;
    line_vector_bp = line_vector_bp * 20;
    line_vector_xyz = inv(Boundary(bNum).View) * [line_vector_bp, 0,0]'; % transform to 2D boundary plane system
    line_vector_xyz = line_vector_xyz(1:3);
    %
    % finish the scatter plot
    figure
    ang = Boundary(bNum).LineAngle;
    tform = affine2d([cosd(ang),-sind(ang),0; sind(ang),cosd(ang),0; 0,0,1]); % define rotation matrix for points
	[trX,trY] = transformPointsForward(tform,inputPTS_x,inputPTS_y); % transform
    scatter(inputPTS_x,inputPTS_y,2); title(['rotated ',num2str(ang)]); pause(0.1) % plot movie of the rotation
    hold on
    quiver(0,0,line_vector_bp(1),line_vector_bp(2),'LineWidth',3,'Color','r')
    hold off
    xlabel('x_{boundary plane}')
    ylabel('y_{boundary plane}')
    title(['2D boundary plane'])
    legend('dislocation point', 'line vector')
    
%}
    % plot all together
    %
    figure
    x_ref = sort(unique(trX));
    x_axis_bp=linspace(min(trX),max(trX),100); % what is the x-axis in the boundary plane system?
    fft_x_axis_bp=linspace(-1/(x_ref(2)-x_ref(1)),1/(x_ref(2)-x_ref(1)),100); % what is the x-axis in the boundary plane system?
%    imagesc(x_axis_bp,[1:1:180],im_HIST); colorbar
    imagesc(fft_x_axis_bp,[1:1:180],abs(transpose(im_FFT))); xlim([0,1/(x_ref(2)-x_ref(1))]); colorbar
    caxis([0,1e4])
    xlabel('frequency probability')
%    xlabel('$x$_{bp}','interpreter','latex')
    ylabel('Angle (^o)')
%    title('1D histogram from all angles of 2D rotations')
    title('Fourier transform result of 2D rotation and 1D histogram')
    set(gca,'ydir','normal','fontsize',18,'FontName','Times New Roman')
    set(gcf,'Color','w')
    %}
    %
    figure
    pcshow(Boundary(bNum).ptCloud); hold on
    refPT = mean(Boundary(bNum).XYZ_Points); 
    normalVEC = Boundary(bNum).Solved_xyz * 15; 
    quiver3(refPT(1),refPT(2),refPT(3),normalVEC(1),normalVEC(2),normalVEC(3),'LineWidth',3);
    plotNormalPlane(refPT,normalVEC',color(bNum,:))
    title('Line vector in 3D point cloud')
    % plot the line vector in 3D point cloud
    quiver3(refPT(1),refPT(2),refPT(3),line_vector_xyz(1),line_vector_xyz(2),line_vector_xyz(3),'LineWidth',3,'Color','g')
    view(normalVEC)
    set(gcf,'color','w')
    set(gca,'ZColor','k','YColor','k','XColor','k','Color','w','FontName','Times New Roman',...
        'FontSize',18)
    box on
    
    %}
	% define line vector and convert coordinate system
    Boundary(bNum).DislocationLineVec_xyz = line_vector_xyz;
    Boundary(bNum).DislocationLineVec_uvw = TransformBasis * Boundary(bNum).DislocationLineVec_xyz;
    Boundary(bNum).DislocationLineVec_hkl = 1./Boundary(bNum).DislocationLineVec_uvw; 
    Boundary(bNum).DislocationLineVec_hkl(Boundary(bNum).DislocationLineVec_hkl==Inf) = 0;

    % print the result
    hklTEXT = [num2str(Boundary(bNum).DislocationLineVec_hkl(1)),',',...
        num2str(Boundary(bNum).DislocationLineVec_hkl(2)),',',...
        num2str(Boundary(bNum).DislocationLineVec_hkl(3))];
    fprintf(['We predict line vector of B',num2str(bNum),' has hkl ',hklTEXT,'.\n'])

    view(Boundary(bNum).DislocationLineVec_xyz)

end



