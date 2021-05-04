close all
clear
clc

%%
choice = menu('Which rhythm do you want to simulate?','Theta','Alpha','Alpha-Beta','Else');
if choice == 1 
    Wep = 270; % theta
elseif choice == 2
    Wep = 135; % alpha
elseif choice == 3
    Wep = 68; % alpha-beta
else
    % Choose between: 108; 128; 675; 1350
    Wep = 108;  
end

Wpe = 0.8*Wep;
Wip = 0.25*Wep;
Wpi = 0.25*Wep; % weight of the inhibitory synapse
Ae = 3.25; % [mV]
Ai = 22; % [mV]              %22;  %17.6;
ae = 100; % = 1/tau_e        % 50 75 150 200
ai = 50; % = 1/tau_i
V0 = 6; % [mV]
k = 0.56; % [1/mV]
rmax = 5; % [1/s]
dt = 0.0001; % [ms]
t = 0:dt:200;
L = length(t);
yp = zeros(1,L);
zp = zeros(1,L);
ye = zeros(1,L);
ze = zeros(1,L);
yi = zeros(1,L);
zi = zeros(1,L);
Vp = zeros(1,L);
Vi = zeros(1,L);
Ve = zeros(1,L);

noise = 200*randn(1,L)+160;

for i = 1:L-1
    Vp(i) = Wpe*ye(i) - Wpi*yi(i);
    rp = rmax/(1+exp(-k*(Vp(i)-V0)));
    yp(i+1) = yp(i) + dt*zp(i);
    zp(i+1) = zp(i) + dt*(Ae*ae*rp-2*ae*zp(i)-ae^2*yp(i));
    
    Ve(i) = Wep*yp(i);
    re = rmax/(1+exp(-k*(Ve(i)-V0)));
    ye(i+1) = ye(i) + dt*ze(i);
    ze(i+1) = ze(i) + dt*(Ae*ae*(re+noise(i)/Wpe)-2*ae*ze(i)-ae^2*ye(i));
    
    Vi(i) = Wip*yp(i);
    ri = rmax/(1+exp(-k*(Vi(i)-V0)));
    yi(i+1) = yi(i) + dt*zi(i);
    zi(i+1) = zi(i) + dt*(Ai*ai*ri-2*ai*zi(i)-ai^2*yi(i)); 
end

start = 30/dt;
stop = 33/dt;
plot(t(start:stop),Vp(start:stop))
xlabel('t [s]')
ylabel('V_p [mV]')
if choice == 1 
   title('Theta rhythm')
elseif choice == 2
    title('Alpha rhythm')
elseif choice == 3
    title('Alpha-Beta rhythm')
else
    title('V_p(t)')
end