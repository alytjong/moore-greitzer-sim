%%  MOORE-GREITZER EQUATION SIMULATION 2020/02/13
% NEW VERSION, NO FINITE DIFFERENCE, FROM STATE-SPACE FORM
% USE FOURIER SERIES ON PDE TO CONVERT IT TO ODE

%%  Parameters
%Params from 2008 MingQing Xiao

params6.l_c = 8; 
params6.B = 0.5; 
params6.m = 1.75;
params6.a = 1/3.5;
params6.nu = 0.1;
params6.H = 0.18;
params6.W = 0.25;
params6.gamma = 0.66;

params7 = params6;
params7.B = 2;
params7.gamma = 0.6;

params8 = params6;
params8.B = 0.5;
params8.gamma = 0.572;

params9 = params6;
params9.B = 0.72061;
params9.gamma = 0.572;

params = params9;

%%  Theta and time domain
L = 2*pi;
n = 512;  %domain, frequency resolution

tht2 = linspace(-L/2, L/2, n+1);
tht = tht2(1:n);

t = 0:1:1000;

k = (2*pi/L)*[0:n/2-1 -n/2:-1]';    %wave number

%%  Initial conditions
g0 = 0.005*sin(tht); %IC of PDE
gt0 = fft(g0)';    %fft of IC, make column vector

Phi0 = 0.51;
Psi0 = 0.66;

Y0 = [gt0; Phi0; Psi0];

%%  Solve, postprocessing
[t, Y] = ode45(@(t,y)MGEOM(t,y,k,params), t, Y0);
%or load data;
%%
load Yfig9.mat

gtsol = zeros(numel(t), numel(tht));
for j = 1:length(t)
    gtsol(j,:) = real(ifft(Y(j,1:end-2)));
end
Phi = Y(:,end-1);
Psi = Y(:,end);

%%  Plots

%Animation to check g data
f = figure('color', 'w');
line1 = plot(tht,gtsol(1,:));
title('Dynamics of $$g(\xi,\theta)$$ in stable configuration', ...
    'interpreter', 'latex', 'fontsize', 18);
xlabel('$$\theta$$' ,'interpreter', 'latex', 'fontsize', 15);
ylabel('$$g(\xi,\theta)$$', 'interpreter', 'latex', 'fontsize', 15);
xticks([-pi 0 pi]);
xticklabels({'-\pi','0','\pi'});
for iii = 1:40%1:size(gtsol,1)
    set(line1, 'ydata', gtsol(iii,:));
    axis([-pi pi -0.1 0.1]);
    drawnow;
    pause(0.01);
    
    %filename = ['ganimfig6_' num2str(iii) '.png'];
    %saveas(f, filename);
end

%%
figure; mesh(tht,t,gtsol);

%%
%Phi and psi
figure('color', 'w');
subplot(2,1,1); plot(t,Phi); ylabel('flow coefficient \Phi');
subplot(2,1,2); plot(t,Psi); ylabel('pressure rise \Psi');
%%
%Phase portrait
figure('color', 'w'); plot(Phi,Psi); hold on;
%unpack params
l_c = params.l_c; 
B = params.B; 
m = params.m;
a = params.a;
nu = params.nu;
H = params.H;
W = params.W;
gamma = params.gamma;
psi_c0 = 1.67*H;
psi_c = @(x) psi_c0 + H*(1 + 1.5.*(x./W-1)- 0.5.*(x./W-1).^3);

%calculate equilibrium
%r = fsolve(@(x) psi_c(x)- (x.^2)./(gamma^2), 0.4);
x = psi_c0*(W^3/H);
y = -(2/3)*(W^3/H)*(3*H/2/W^2 - 1/gamma^2);
sqroot = sqrt(x*(x-2*y^3));
Phie_a = nthroot(x-y^3 + sqroot,3) + nthroot(x-y^3 - sqroot,3) - y;

Phie = Phie_a;
Psie = psi_c(Phie);

Delta = Psie/Phie - a/4/B^2/nu;

plot(Phi,psi_c(Phi));
plot(Phi, (1/gamma^2)*Phi.^2);
plot(Phie,Psie,'rx');
plot(Phi(1),Psi(1),'bo');
legend('Phase portrait', ...
    'Compressor characteristic $$\psi_c$$', ...
    'Throttle characteristic $$F_T$$', ...
    'Equilibrium',...
    'Initial condition',...
    'location', 'southeast',...
    'interpreter', 'latex',...
    'fontsize', 12);
axis equal;
%title(['\Delta = ' num2str(Delta)]);

title('Dynamics of $$\Phi$$ and $$\Psi$$ in surge and stall combination', ...
    'interpreter', 'latex', 'fontsize', 18);
ylabel('$$\Psi$$', 'interpreter', 'latex', 'fontsize', 15);
xlabel('$$\Phi$$', 'interpreter', 'latex', 'fontsize', 15);

