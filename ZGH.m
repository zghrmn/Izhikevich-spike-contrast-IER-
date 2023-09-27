% Created by Eugene M. Izhikevich, February 25, 2003
% Excitatory neurons   Inhibitory neurons

%Sep_K_edits for Z_GH

%The number of excitatory neurons in the network.  The mammalian cortex has
%about 4 times as many excitatory nerons as inhibitory ones.


%This code has a constant ratio of 4 to 1 _ Exc to Inh number of neurons

clear

clc


Emat = 0.1:0.5:1.1;
Imat = 0.1:0.5:1.1;




for EE = 1:length(Emat)
    Eweight = Emat(EE);
    for II = 1:length(Imat)
        
        firings =[];
        firings2 = [];
        Iweight = Imat(II);
        Ne=800;                Ni=200;
        re=rand(Ne,1);         ri=rand(Ni,1);
        %This will set the value of a for all excitatory neurons to 0.02 and the
        %value of a for inhibitory neurons to a random number between 0.02 and 0.1
        a=[0.02*ones(Ne,1);    0.02+0.08*ri];
        %This will allow b to range from 0.2-0.25
        b=[0.2*ones(Ne,1);     0.25-0.05*ri];
        %This will allow the spike reset membrane potential to range between -65
        %and -50
        c=[-65+15*re.^2;      -65*ones(Ni,1)];
        %This will allow the recovery reset value to range between 2 and 8
        d=[8-6*re.^2;          2*ones(Ni,1)];
        S=[Eweight*rand(Ne+Ni,Ne),Iweight*(-rand(Ne+Ni,Ni))];
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %using this part odf the code you can set the percentage value of
        %neurons NOT connected at all. If you want all of them connected to
        %eachother set it to zero
        percent_off=0.2;
        connections=randperm((Ne+Ni)^2);
        connections=connections(1:(floor(percent_off*length(connections))));
        for i=1:length(connections)
            S(connections(i))=0;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %The initial values for v and u
        v=-65*ones(Ne+Ni,1);  % Initial values of v
        u=b.*v;               % Initial values of u
        %Firings will be a two-column matrix.
        %The first column will indicate the time (1-1000)
        %that a neuron’s membrane potential crossed 30, and
        %the second column %will be a number between 1 and Ne+Ni
        %that identifies which neuron fired at that %time.
        %firings=[];
        firings=[];           % spike timings
        Simulationduration = 1000;
        for t=1:Simulationduration
            %Create some random input external to the network
            I=[5*randn(Ne,1);2*randn(Ni,1)]; % thalamic input
            %Determine which neurons crossed threshold at the
            %current time step t.
            fired=find(v>=30); % indices of spikes
            if ~isempty(fired)
                %Add the times of firing and the neuron number to firings.
                firings=[firings; t+0*fired, fired];
                %Reset the neurons that fired to the spike reset membrane potential and
                %recovery variable.
                v(fired)=c(fired);
                u(fired)=u(fired)+d(fired);
                %strengths of all other neurons that fired in the last time step connected to that
                %neuron.
                I=I+sum(S(:,fired),2);
            end
            %Move the simulation forward using Euler’s method.
            %             v=v+0.5*(0.04*v.^2+5*v+140-u+I);
            v=v+0.5*(0.04*v.^2+5*v+140-u+I);
            u=u+a.*(b.*v-u);
        end
        %Plot the raster plot of the network activity.
        
        
        
        
        
        firings2(:,1) = firings(:,2);
        firings2(:,2) = firings(:,1);
        
        ffs = zeros(Simulationduration,Ne+Ni);
        qq = unique(firings2(:,2));
        for xx = 1:numel(qq); idx = firings2(find(firings2(:,2)==qq(xx)),2);
            for j = 1:numel(idx)
                ffs(qq(xx), j) = idx(j);
            end
        end
        SS = SpikeContrast(ffs,Simulationduration);
        scsync(EE,II) = SS;
        A3 = SS;
        A1 = Eweight;
        A2 = Iweight;
        
        %Plot the raster plot of the network activity.
        figure()
        plot(firings(:,1),firings(:,2),'.');
        formatSpec = "Exc Weight: %d _ Inh Weight: %d _ SCSynchrony: %d ";
        
        title(sprintf(formatSpec,A1,A2,A3));
        
        
        
        
        
        
        
        
    end
end


imagesc(scsync)








