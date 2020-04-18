% Berken Utku Demirel - 2166221
clearvars
close all
%%
% Duration
duration = linspace(0,300,2001); 
durationStep = duration(2)-duration(1); 

% a)
Vmembrane = -90*ones(1,length(duration)); 
Vmembrane(duration>10) = 0;

% b)
% Vmembrane = -90*ones(1,length(duration)); 
% Vmembrane(duration>200) = -90;


% c)
% Vmembrane = -90*ones(1,length(duration)); 
% Vmembrane(duration>10) = -80;
% Vmembrane(duration>200) = -90;


%define constants
GNArest = 120; 
GKrest = 36; 
GL = 0.3; 
Cm = 1; 
Vrest = -90; 
VNA = 25; 
VK = -102;  
VL = -79.387;  

%definition of initials
m = 0.05*ones(1,length(duration));  
h = 0.54*ones(1,length(duration));  
n = 0.34*ones(1,length(duration));  
GNA = GNArest*ones(1,length(duration)); 
GK = GKrest*ones(1,length(duration)); 
alpha_m = 2.237*ones(1,length(duration));  
alpha_n = 0.0582*ones(1,length(duration)); 
alpha_h = 0.07*ones(1,length(duration)); 
beta_m = 4*ones(1,length(duration));  
beta_n = 0.125*ones(1,length(duration));
beta_h = 0.0474*ones(1,length(duration));
%% Loop for simulation

for i = 2:length(duration)

alpha_m(i) = 0.1*(25-Vmembrane(i-1)+Vrest)/(exp((25-Vmembrane(i-1)+Vrest)/10)-1);
alpha_h(i) = 0.07/exp((Vmembrane(i-1)-Vrest)/20);
alpha_n(i) = 0.01*(10-Vmembrane(i-1)+Vrest)/(exp((10-Vmembrane(i-1)+Vrest)/10)-1);
beta_m(i) = 4/exp((Vmembrane(i-1)-Vrest)/18);
beta_h(i) = 1/(exp((30-Vmembrane(i-1)+Vrest)/10)+1);
beta_n(i) = 0.125/exp((Vmembrane(i-1)-Vrest)/80);

%ionic conductances
m(i) = m(i-1) + durationStep*(alpha_m(i-1)*(1-m(i-1))-beta_m(i-1)*m(i-1));
h(i) = h(i-1) + durationStep*(alpha_h(i-1)*(1-h(i-1))-beta_h(i-1)*h(i-1));
n(i) = n(i-1) + durationStep*(alpha_n(i-1)*(1-n(i-1))-beta_n(i-1)*n(i-1));
GNA(i) = GNArest*m(i-1)^3*h(i-1);
GK(i) = GKrest*n(i-1)^4;

end
%% Plots

figure, plot(duration,alpha_n); 
title('alpha n');
hold on
plot(duration,beta_n);
title('alpha n and beta n');
xlabel('duration  (ms)');
legend('alpha_n','beta_n')

figure,
plot(duration,alpha_m);
hold on
plot(duration,beta_m);
title('alpha m and beta m');
xlabel('duration  (ms)');
legend('alpha_m','beta_m')

figure, 
plot(duration,alpha_h);
hold on
plot(duration,beta_h); 
title('alpha h and beta h');
xlabel('duration  (ms)');
legend('alpha_h','beta_h')

figure, 
plot(duration,n); 
hold on
plot(duration,m); 
title('n and m');
xlabel('duration  (ms)'); 
legend('n','m')

figure, 
plot(duration,h); 
title('h');
xlabel('duration  (ms)');
