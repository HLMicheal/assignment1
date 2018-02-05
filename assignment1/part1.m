% By Huanyu Liu
% 100986552
% For elec4700 assignment1 part1

% Initialize the parameters
n=10; % number of particles
T=300; % temperture of the backgound
L=200e-9; % length of the frame (figure 1)
H=100e-9; % height of the frame
tao=0.2e-12; % the given mean time between collisions
m0=9.109e-31; % mass of a particle
mn=0.26*m0; % effective mass
kb=1.38e-23; % constant coeffient
vth=sqrt(2*kb*T/mn); % average speed of each particle
% Initialize the positions of each particle
Pox = L*rand(1,n);
Poy = H*rand(1,n);
% Initialize the speed of each particle and measure the initial temperature
for num=1:n
Vx(num) = randn()*vth/sqrt(2);
Vy(num) = randn()*vth/sqrt(2);
end
Tmeasured = sum(Vx.^2 + Vy.^2).*mn./(2*kb*n);
% draw the first locations of the particles and the blocks
figure(1)
plot(Pox,Poy,'.');
xlim([0 L]);
ylim([0 H]);
hold on
% more parameters that will be used in the loop
TStop = 1e-11; % max running time
t=0; % start time
dt=1e-14; % step time
path=0; % initialize the path length
 while t < TStop
     z=round(1+t/dt); % index, the z-th interval between collisions
     
         Tmeasured = sum(Vx.^2 + Vy.^2).*mn./(2*kb*n);
         Vact=sqrt(sum(Vx.^2+Vy.^2)/n); % the average speed of all the particles
         path=path+Vact*dt;
     tPx = Pox + Vx.*dt; % predict the position 
     tPy = Poy + Vy.*dt;
% when the particles go to the right and left border
     px1 = Pox >= L;
     Pox(px1) = Pox(px1) - L;
     px2 = Pox <= 0;
     Pox(px2) = Pox(px2) + L;
     
     py1 = tPy <= 0;
     Vy(py1) = Vy(py1) .* (-1);
     py2 = tPy >= H;
     Vy(py2) = Vy(py2) .* (-1);
         % now all velocity have been modified to the correct direction,
         % update the position
     Pox = Pox + Vx.*dt;
     Poy = Poy + Vy.*dt;
     
     figure(1)
     plot(Pox,Poy,'.');
     xlim([0 L]);
     ylim([0 H]);
     hold on

     figure(2)
     plot(t,Tmeasured,'.r');
     title('temperature plot');
     hold on
     
     fprintf('time: %g (%5.2g %%) temperature: %g \n', t/dt, t / TStop * 100, Tmeasured);
     pause(0.01)
     t=t+dt;
 end
     fprintf(' vth: %g\n',vth);
         MFP=path;
     fprintf(' MFP:%g m \n ', MFP);
  