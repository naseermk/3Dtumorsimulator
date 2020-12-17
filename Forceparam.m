function [forces,gamma3,pressure] = Forcepara(coords,rad,poisson,modulus,...
                                             nPart,receptor,ligand)


etab=0.0001;

% Initialize all forces to 0
forces = zeros(size(coords));
gamma3 = zeros(size(rad));
pressure = zeros(size(rad));

% Get the number of particles
%nPart = size(coords,2);

f = 0.0001; %adhesion f^{ad} parameter

for part=1:nPart %change to parfor to call parallel looping in Matlab
%instead of looping over all pairs in the force calculation, find the distance between all
%pairs

dlist= sqrt(sum(bsxfun(@minus, coords(:,part), coords).^2, 1));

%index each set of nearest distances
[d, ind] = sort(rad(part,1)+rad'-dlist);
                
                %begin_index=find(d==2*rad(part,1));
                begin_index=find(d>0.0);
                
                %ind_closest = ind(2:(min(nPart,26))); %find the n nearest neighbors
                ind_closest=ind((begin_index):(end));
                
                coords_closest = coords(:,ind_closest);
                %rad_closest=rad(ind_closest,1);
                %poission_closest=poisson(partA,1)
                for partA=1:(size(ind_closest,2))
                
                dr =  coords(:,part) - coords_closest(:,partA);
                
                if norm(dr) > 0.0
                
                Rij = norm(dr);
                
                RijHat = dr/Rij;
                
                hij = (rad(ind_closest(partA),1) + rad(part,1) - Rij);
                
                Eij = ((1 - poisson(ind_closest(partA),1)^2)/(modulus(ind_closest(partA),1)) ...
                       + (1 - poisson(part,1)^2)/(modulus(part,1)));
                
                Rijf = (1/rad(part,1) + 1/rad(ind_closest(partA),1));
                
                areaint = pi*(1/Rijf)*hij;

                invDr2 = (hij)^(3/2); % 1/r^2
                
                forceFact = (invDr2/(0.75*(Eij)*sqrt(Rijf)))-areaint*f...
                *0.5*(receptor(ind_closest(partA),1)*ligand(part,1)+ ...
                      receptor(part,1)*ligand(ind_closest(partA),1));
                
                
                pressure(part,1) = pressure(part,1)+ ...
                abs(forceFact/areaint);
                
                
                forces(:,part) = forces(:,part) + (forceFact*RijHat);
                
                end
                
                end
                
                
                for partA=1:(size(ind_closest,2))
                
                dr =  coords(:,part) - coords_closest(:,partA);
                
                if norm(dr)>0.0
                
                Rij = norm(dr);
                
                RijHat = dr/Rij;
                
                hij = (rad(ind_closest(partA),1) + rad(part,1) - Rij);

                
                Rijf = (1/rad(part,1) + 1/rad(ind_closest(partA),1));
                
                mag=norm(forces(:,part));
                
                areaint = pi*(1/Rijf)*hij;
                
            gamma2 = etab*areaint*0.5*(1+(sum(bsxfun(@times, forces(:,part),RijHat)))/mag)...
            *0.5*(receptor(ind_closest(partA),1)*ligand(part,1)+ ...
            receptor(part,1)*ligand(ind_closest(partA),1));
 
                gamma3(part,1) = gamma3(part,1) + gamma2;
                
                end
                
                
                end
                
                
                
                
                end
                
                
end
