clear all;
close all;
 
rng('shuffle');

nSteps = 65000;

Polar=zeros(nSteps,1);

scat_fun=zeros(nSteps,1);  
 
cyclenum =1;
 
ds_msd=zeros(nSteps,cyclenum); 
 
radG=zeros(nSteps,cyclenum); 
 
den=zeros(nSteps,cyclenum);
 
df=zeros(nSteps,cyclenum);
 
ds_msdav=zeros(nSteps,1); 
 
radG_av=zeros(nSteps,1); 
 
den_av=zeros(nSteps,1);
 
df_av=zeros(nSteps,1);

liferows=1;
visdata=zeros(liferows,15);
visdata_row=zeros(liferows,15);
 
for cyclemsd=1:cyclenum
    
   numin=100;
   
    % ===================
    %     Initialize
    % ===================
    
    % Set configuration parameters
    nPart = numin;        % Number of particles
    density = 0.001;     % Density of particles
    mass = 1;           % Particles' mass
    nDim = 3;           % The dimensionality of the system (3D in our case)
    
    % Set simulation parameters
    dt = 10.0;         % Integration time
    dt2 = dt*dt;        % Integration time, squared
    
    %nSteps = 3800;       % Total simulation time (in integration steps)
    sampleFreq = 100;    % Sampling frequency
    sampleCounter = 1;  % Sampling counter
    printFreq = 1000;     % Printing frequency
    plotFreq = 100;  %plot frequency
    %volrate = 1.4; %rate of vol growth dv/dt
    radmitosis = 5.0; %mitotic radius aka critical size
    samplecounter = 1;
    eta=0.005; %ECM Viscosity
    gfac=1.0;
    taumin = gfac*54000; % second - minimum cell cycle time
    
    track(nDim,nSteps+1,numin)=zeros;
      
    % Set initial configuration
    [coords L] = initCubicGrid(nPart,density);
    
    % Set initial velocities with random numbers
    vels = zeros(nDim,nPart);
    
    
    %randomly assign a radius to each particle
    rad = zeros(nPart,1);
    modulus = zeros(nPart,1);
    poisson = zeros(nPart,1);
    receptor = zeros(nPart,1);
    ligand = zeros(nPart,1);
    countdr=0; 

    % The neutral mutations will be stored in a vector of vectors: 
    % each cell will inherit a list of mutations from its parent, 
    % and one new unique one. A "mutation" is represented as a number, 
    % this number is just an index used to track the drift of mutations 
    % through the population.

    UpperBound = 10^5;
    genotypes = cell(1,UpperBound);

    % Key: 0 = not a mutation (blank)
    %      1 = mutation present in all cells
    %      2 = first mutation to occur after division
    %      (etc)
    % for example, before division, a cell may have the mutations
    % [1,2,5,7]
    % and afterwards, both it and its daughter have
    % [1,2,5,7,23]
    % Following B Waclaw, the index of the latest new mutation (in the example, 23) 
    % should always be a unique integer: we will set it to the largest
    % index so far +1.

    latest_mutation = 1;

    % We will also need to update both the parents and the daughter genomes
    % to be something like
    % new_genome = [old_genome,[latest_mutation+1]]
    
    for part = 1:nPart
        rad(part,1) = randgaussrad(4.5,0.5); 
        modulus(part,1) = randgaussrad(10^-3,10^-4);
        poisson(part,1) = randgaussrad(0.5,0.02);
        receptor(part,1) = randgaussrad(0.9,0.02);
        ligand(part,1) = randgaussrad(0.9,0.02);
        lifetime(part,1) = 0;
        label(part,1)=part;

        genotypes{1,part} = [1]; % the first few cells are assumed to be clonal
        % and all contain the same mutation
    end
   
    
    % ===================
    % Molecular Dynamics
    % ===================
    
    time = 0; % Following simulation time
    
    volrate = (2*pi*(radmitosis)^3)/(3*taumin);
    for avini=1:numin
    initial(:,avini)=coords(:,avini);
    end
    
    for avini=1:numin
    track(:,1,avini)=coords(:,avini);
    end
      

    no_of_cells_ever_involved =  nPart;  
    
    for step = 1:nSteps
      
        centerM = zeros(3,1);
        % === Calculate new forces ===
        [forces,gamma3,pressure] = Forcepara(coords,rad,poisson,modulus,...
            nPart,receptor,ligand);
        
        % === Second integration step ===
        
        gammat = 0;
        % Implement the Andersen thermostat
        
	    for part =1:nPart
            gammat = 6*pi*eta*rad(part,1) + gamma3(part,1);
            %dropping the inertial term update the coords
            coords(:,part) = coords(:,part) + (dt*forces(:,part))/(gammat);
            % Update velocities - All velocities are updated at once
            vels(:,part) = forces(:,part)/(gammat);
            lifetime(part,1) = lifetime(part,1)+ dt; 
        end
         
        % === First integration step ===
        
        %update the size of particles - all radii updated at one time
        death = 10^(-6);
        deadpart=0;
	      deadnumin=0;
        
        for part=1:nPart
            
            if rand <= death*dt
                
                deadpart=deadpart+1;
                deadindex(deadpart,1) = part;
                
            elseif (rad(part,1) < radmitosis) && (pressure(part,1) < 0.0001)
                
                grate = (volrate/(4*pi*rad(part,1)*rad(part,1)));
                rad(part,1) = rad(part,1) + dt*randgaussrad(grate,10^-5);
                
                
            elseif (rad(part,1) >= radmitosis) && (pressure(part,1) < 0.0001)
                
                rad(end+1,1)=(2^(-1/3))*rad(part,1);
                rad(part,1) = (2^(-1/3))*rad(part,1); %new radius after division
                modulus(end+1,1) = randgaussrad(10^-3,10^-4);
                poisson(end+1,1) = randgaussrad(0.5,0.02);
                 
                
                receptor(end+1,1) = randgaussrad(0.9,0.02);
                ligand(end+1,1) = randgaussrad(0.9,0.02);
                

                lifetime(end+1,1) = 0; 
                no_of_cells_ever_involved= no_of_cells_ever_involved +1;
                label(end+1,1)= no_of_cells_ever_involved; 
                
                %must now add elements to vels, forces as well
                %due to coords=coords+dt*vels + 0.5*dt2*forces equation
                vels(:,end+1)= ([0,0,0]');
                
                %generating random numbers between 0,1
                a=0;
                b=1;
                r3=(b-a).*rand(1) + a;
                r4=pi*(b-a).*rand(1) + pi*a;
                r5=2*pi*(b-a).*rand(1) + 2*pi*a;
                psi=size(coords,2)+1;
                
                coords(1,psi) = coords(1,part)+radmitosis*(1-2^(-1/3))*sin(r4)*cos(r5);
                coords(2,psi) = coords(2,part)+radmitosis*(1-2^(-1/3))*sin(r4)*sin(r5);
                coords(3,psi) = coords(3,part)+radmitosis*(1-2^(-1/3))*cos(r4);
                
                coords(1,part)=coords(1,part)-radmitosis*(1-2^(-1/3))*sin(r4)*cos(r5);
                coords(2,part)=coords(2,part)-radmitosis*(1-2^(-1/3))*sin(r4)*sin(r5);
                coords(3,part) = coords(3,part)-radmitosis*(1-2^(-1/3))*cos(r4);

                % add neutral mutation to parent and daughter:
                % TODO test and finish
                old_genome = genotypes{1,part};
                new_genome = [old_genome,[latest_mutation+1]];
                genotypes{1,part} = new_genome;
                genotypes{1,end+1} = new_genome;
                latest_mutation = latest_mutation + 1;
                
            end
            
            
        end

        if deadpart > 0 
            
            coords(:,deadindex(1:deadpart))=[];
            rad(deadindex(1:deadpart))=[];
            modulus(deadindex(1:deadpart)) = [];
            poisson(deadindex(1:deadpart)) = [];
            receptor(deadindex(1:deadpart)) =[];
            ligand(deadindex(1:deadpart)) = [];
            vels(:,deadindex(1:deadpart))= [];
                lifetime(deadindex(1:deadpart))=[];
            label(deadindex(1:deadpart))= [];        
	
	        for part=1:deadpart
		   
                if deadindex(part,1)<= numin
                    deadnumin=deadnumin+1;
                    deadindexn(deadnumin,1)=deadindex(part,1);
                end
            
            end
        
        end
        
        if deadnumin > 0
           track(:,:,deadindexn(1:deadnumin,1))= [];
           initial(:,deadindexn(1:deadnumin,1))=[];            
        end

        numin=size(initial,2);

        
        nPart=size(coords,2);
        %this is the radius of gyration squared
        
	      centerM(1,1) = mean(coords(1,:));
        centerM(2,1) = mean(coords(2,:));
        centerM(3,1) = mean(coords(3,:));

	      for part=1:nPart
            radG(step,cyclemsd) = radG(step,cyclemsd) + ...
                ((norm(coords(:,part) - centerM(:,1)))^2)/nPart;  
        end
        
        den(step,cyclemsd) =nPart/((1/3)*4*pi*radG(step,cyclemsd)^(3/2));  
        
        radG_av(step) = radG_av(step) + (radG(step,cyclemsd)/cyclenum);
        
        den_av(step) = den_av(step) + (den(step,cyclemsd)/cyclenum);
        
        % === Move time forward ===
        time = time + dt;
        
        if mod(step,printFreq) == 0
            step % Print the step
            cyclemsd
            latest_mutation % print latest mutation
            time
        end
        
        if mod(step,plotFreq) == 0
            
            numcells(samplecounter,cyclemsd) = nPart;
            samplecounter = samplecounter + 1;
          
        end
        
        
        for avini=1:numin
            ds_msd(step,cyclemsd)=(norm(coords(:,avini)-initial(:,avini)))^2;  
            ds_msdav(step) = ds_msdav(step) + (ds_msd(step,cyclemsd))/(cyclenum*numin);
            
            df(step,cyclemsd)=(norm(coords(:,avini)-initial(:,avini)))^4;
            df_av(step)=df_av(step)+(df(step,cyclemsd))/(cyclenum*numin);
        end
        
        for avini=1:numin
            track(:,step+1,avini)=coords(:,avini);  
        end
        
        save('fnumin.txt','numin','-ascii','-append');
        
	    if mod(step,plotFreq) == 0
            countdr = countdr+1;
            deltart(countdr,cyclemsd) = rtumor(coords,centerM);
                for part=1:nPart
                     visdata(1,1:3)=coords(1:3,part);
                     visdata(1,4)= label(part,1);
                     visdata(1,5)= lifetime(part,1);
                     visdata(1,6)= step*dt;
                     visdata(1,7)= receptor(part,1);
                     visdata(1,8)= ligand(part,1);
                     visdata(1,9)= modulus(part,1);
                     visdata(1,10)= poisson(part,1);
                     visdata(1,11)= rad(part,1);
                     visdata(1,12:14)=vels(1:3,part);
                     visdata(1,15)=cyclemsd;
                     visdata_row= visdata(1,:);
                     save('lifetime1.txt','visdata_row', '-ascii','-append');
                end
        end	 

    end
 
    % compute and write out mutation frequency distribution
    freqs = zeros(latest_mutation,2);
    for gene = 1:latest_mutation
      count = 0;
      for part = 1:UpperBound
          this_genome = genotypes{1,part};
          for elem = 1:length(this_genome)
              if elem == gene
                  count = count + 1;
              end
          end
      end
      freqs(gene,:) = [gene count];
    end   

    freqs
    
    save('gene_freqs.txt', 'freqs', '-ascii');
    %writecell(genotypes, 'genotypes.txt', 'Delimiter', 'tab');
    
    % Simulation results
    % ===================
    radG_inst(1,1) = sqrt(radG(step,cyclemsd));
    rho = nPart/((4/3)*pi*radG_inst^3);
    
    L=2*radG_inst(1,1);
    dL = 1.0;
    
    save('initial.txt','initial','-ascii','-append');    
end
 
% Other simulation results: 
 

  save('numcellf.txt','numcells','-ascii');
  save('fds_msd.txt','ds_msd','-ascii');
  save('fds_msdav.txt','ds_msdav','-ascii');
  save('fradG.txt','radG','-ascii');
  save('fradG_av.txt','radG_av','-ascii');
  save('fden.txt','den','-ascii');
  save('fden_av.txt','den_av','-ascii');
  save('radius.txt','rad','-ascii');
  save('deltart.txt','deltart','-ascii'); 

