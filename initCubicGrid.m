function [coords, L] = initCubicGrid(nPart,density)

% Initialize with zeroes
coords = zeros(3,nPart);

% Get the corresponding box size
L = 40;


% Start positioning - use a 3D index for counting the spots
index = [0,0,0]';

% Assign particle positions
%need to change this to take into account birth position + distance
parameter=1;
		
for part=1:nPart
			
parameter=parameter+1;

% Set coordinate

	coords(:,part) = (index+[L*abs(randn),L*abs(randn),L*abs(randn)]');
                               
end
                                                        
end
