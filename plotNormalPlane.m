function plotNormalPlane(PT,vector,color)
% plot the plane normal to the current pointCloud, based on inputs:
%     PT      the reference point at the start of the vector
%     vector  the normal vector to reference against
%     color   the color of these objects

    % define the location of the reference point locally
   x1=PT(1)*1.0;  y1=PT(2)*1.0;  z1=PT(3)*1.0;

   % find the two orthonormal vectors that form a set with "vector"
   sz = size(vector);
   if sz(1) > 1
       vector = vector';
   end
   w = null(vector); 
   
   % set the grid for the plane and plot it
   [P,Q] = meshgrid(-100:100); % Provide a gridwork (you choose the size)
   X = x1+w(1,1)*P+w(1,2)*Q; % Compute the corresponding cartesian coordinates
   Y = y1+w(2,1)*P+w(2,2)*Q; %   using the two vectors in w
   Z = z1+w(3,1)*P+w(3,2)*Q;
   surf(X,Y,Z,'FaceAlpha',0.3,'FaceColor',color,'EdgeColor','none')
end