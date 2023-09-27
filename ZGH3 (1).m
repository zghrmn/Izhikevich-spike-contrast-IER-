% Created by Eugene M. Izhikevich, February 25, 2003
% Excitatory neurons   Inhibitory neurons

%Sep_K_edits for Z_GH




%This code has a different ratios of Exc to Inh ratio, synchrony evaluated 
clear

clc


for Eweight = 0.5:0.1:1.1

% Eweight = 0.5;
Iweight = 0.8;

E_to_I_prop = 0.1:0.01:0.9;

total_number_of_neurons = 4000;



for EIP = length(E_to_I_prop):length(E_to_I_prop)
    
        Ne =[];
        Ni = [];
    
        EIprop = E_to_I_prop(EIP);
        
        Ne = round(EIprop*total_number_of_neurons);
        Ni = total_number_of_neurons - Ne;
        
        firings =[];
        firings2 = [];
        
        
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
        %using this part of the code you can set the percentage value of
        %neurons NOT connected at all. If you want all of them connected to
        %eachother set it to zero
        percent_off=0.6; %arbitrary connection probability)
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
        Simulationduration = 4000;
        
        dt = 1;
        
        for t=1:Simulationduration/dt
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
            v=v+0.5*(0.04*v.^2+5*v+140-u+I).*dt;
            u=u+a.*(b.*v-u).*dt;
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
        scsync(EIP) = SS;
        A3 = SS;
        A1 = Eweight;
        A2 = Iweight;
        
        %Plot the raster plot of the network activity.
        figure()
        plot(firings(:,1),firings(:,2),'.');
        
end
end
        
        %plot the effect of exc to inh proportion
        
%         figure()
%         
%         plot(E_to_I_prop,scsync); xlabel('Excitatory to iInhibitory ratio'); ylabel('Contrast Synchrony')









