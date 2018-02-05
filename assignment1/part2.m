% By Huanyu Liu
% 100986552
% For elec4700 assignment1 part2

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
intervals=round(TStop/dt); % number of steps 
Vz=zeros(1,intervals); % initial the size of all changing speed (will be used in hist)
ddt = 0; % time since last timestop
collisions=0; % number of timestops
time=0; % initialize the duration between collisions
path=zeros(1,n); % initialize the size of path length
 while t < TStop
     z=round(1+t/dt); % index, the z-th interval between collisions
     Pscat = 1-exp(-ddt/tao); % scattering posibility
     if Pscat > rand % if scatter
         time=time+ddt; % total time when scattering occur
         ddt=0; % reset the parameter for the possibility as required
         collisions=collisions+1; % one more collision occurs
         Vx = randn(1,n).*vth/sqrt(2);
         Vy = randn(1,n).*vth/sqrt(2); % velocity changes (in maxwell-boltzmann distribution)
         average_path_length(collisions)=sum(path)/n; % average path length for this interval
         path=zeros(1,n); % reset the path length
     else % nothing happens, same speed the next duration of time step
         path=path+sqrt(Vx.^2+Vy.^2).*dt; % add the next timestep's path length to the total path length
         ddt=ddt+dt; % add the timestep size to the parameter
     end
         Tmeasured = sum(Vx.^2 + Vy.^2).*mn./(2*kb*n);
         Vact=sqrt(sum(Vx.^2+Vy.^2)/n); % the average speed of all the particles
         Vz(z)=Vact; % will be used to get the distribution in hist
         
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
         figure(3)
         hist(Vz);
         title('velocity histogram');
         xlabel('velocity');
     fprintf(' vth: %g\n',vth);
         mean_time_between_collisions=time/collisions;
         MFP=sum(average_path_length)/collisions;
     fprintf(' MFP:%g m \n mean time between collisions: %g s\n', MFP, mean_time_between_collisions);
  