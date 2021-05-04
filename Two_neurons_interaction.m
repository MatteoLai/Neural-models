close all
clear
clc

%%
E0 = -65; % [mV] 
Vreset = -65; % [mV]
tau_m = 30; % [ms] - time constant of the membrane
tau_th = 5; % [ms] - time constant of the dynamic threshold (to model the refractory period)
r = 10; % [Mohm]
C = tau_m/r;
g = 1/r;
gsmax = 5/r;  %gs = gsmax*Ps
Ea = -90; % [mV]
Es1 = -70; % [mV]   - 1 : excitatory
Es2 = 0; % [mV]     - 2 : inhibitory 
% Note: There are two neurons in this model: one excitatory (1), and one
% inhinibitory (2). When 1 fires, it excites 2 (inhinibitory),then the 
% excitatory channels will open, which have Nernst potential Es2 = 0mV.
% Vice versa, the excitatory neuron will have Nernst potential Es1 = -70mV.

DPs1 = 0.3; % (0.03 - 0.6)
DPs2 = 0.2;
tau_s1 = 10; % [ms]
tau_s2 = 10; % [ms]

i1 = 4;
i2 = 0;


dt = 0.01;
t = 0:dt:100; % [ms]
L = length(t);

v1 = zeros(1,L);
v1(1) = E0;
v2 = zeros(1,L);
v2(1) = E0;
VthL = -55; % [mV]
VthH = 50;  % [mV]
Vth1 = zeros(1,L);
Vth1(1) = VthL;
Vth2 = zeros(1,L);
Vth2(1) = VthL;
spike1 = zeros(1,L); 
spike2 = zeros(1,L);
Pi1 = zeros(1,L);
Pe2 = zeros(1,L);

for k = 1:length(t)-1
    gi1 = gsmax*Pi1(k);
    ge2 = gsmax*Pe2(k);
    req1 = 1/(g+gi1);
    req2 = 1/(g+ge2);
    Eeq1 = (g*E0 + gi1*Es1)/(g+gi1); % Excitatory (1) neuron potential,
                 % which receives synaptic input from the inhibitory neuron
                 % (so it will be an inhibitory synapse in the model)
    Eeq2 = (g*E0 + ge2*Es2)/(g+ge2); % Inhibitory (2) neuron potential,
                 % which receives synaptic input from the excitatory neuron
    tau1 = req1*C;
    tau2 = req2*C;
    vinf1 = Eeq1 + req1*i1;
    vinf2 = Eeq2 + req2*i2;
    
    % Update the membrane potentials:
    v1(k+1) = (v1(k) - vinf1)*exp(-dt/tau1) + vinf1;
    v2(k+1) = (v2(k) - vinf2)*exp(-dt/tau2) + vinf2;
    
    % Update the threshold of the model:
    Vth1(k+1) = (Vth1(k) - VthL)*exp(-dt/tau_th) + VthL;
    Vth2(k+1) = (Vth2(k) - VthL)*exp(-dt/tau_th) + VthL;
    
    % Update of the probability of ion channels opening:
    Pi1(k+1) = Pi1(k)*exp(-dt/tau_s1);
    Pe2(k+1) = Pe2(k)*exp(-dt/tau_s2);
    
    if v1(k+1) > Vth1(k+1)
        % If the potential of neuron 1 (excitatory) exceeds the threshold, it
        % cause the spike (reset v1 and Vth1), and excite neuron 2 (by
        % modifying P2, opening probability of its channels)
        % Note: g2 = gmax*P2 --> if P2 increase, conductance also increases)
        v1(k+1) = Vreset;
        Vth1(k+1) = VthH;
        Pe2(k+1) = Pe2(k+1) + DPs2*(1-Pe2(k+1)); 
        spike1(k) = 1;
    end
    if v2(k+1) > Vth2(k+1) 
        % If the potential of neuron 2 (inhibitory) exceeds the threshold, it
        % cause the spike (reset v2 and Vth2), and inhibit neuron 1 (by
        % modifying P1, opening probability of its channels)
        % Note: g1 = gmax*P1 --> if P1 decrease, conductance also decrease)
        v2(k+1) = Vreset;
        Vth2(k+1) = VthH;
        Pi1(k+1) = Pi1(k+1) + DPs1*(1-Pi1(k+1));
        spike2(k) = 1;
    end
end

subplot(2,1,1)
plot(t,v1,t,Vth1,'r--')
xlabel('t [s]')
legend('Membrane potential','Threshold','Location','west')
title('Excitatory neuron')
subplot(2,1,2)
plot(t,v2,t,Vth2,'r--')
xlabel('t [s]')
legend('Membrane potential','Threshold','Location','west')
title('Inhibitory neuron')

figure
plot(t,v1,'g',t,v2,'r')
hold on
idx1 = find(spike1 == 1);
idx2 = find(spike2 == 1);
plot(idx1*dt,v1(idx1),'*g',idx2*dt,v2(idx2),'*r')
legend({'Excitatory (1)','Inhibitory (2)'},'Location','northwest')
xlabel('t [s]')

% Note: every time the excitatory neuron (1) fires - see asterisks-, it excites
%  2, so its membrane potential rises more rapidly towards the threshold.
% Vice versa, every time the inhibitory neuron (2) fires, it inhibits 1, 
% therefore it will take longer to reach the potential threshold.

figure
subplot(2,1,1)
plot(t,spike1)
xlabel('t [s]')
ylabel('spike')
title('Excitatory neuron')
subplot(2,1,2)
plot(t,spike2)
xlabel('t [s]')
ylabel('spike')
title('Inhibitory neuron')