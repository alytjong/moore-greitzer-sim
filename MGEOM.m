function dY = MGEOM(t, Y, k, params)

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
psi_c2 = @(g,x) H*(1.5*(g./params.W) - 0.5*(g./params.W).^3  ...
    - 1.5*((x./params.W-1).^2).*(g./params.W) ...
    - 1.5*(x./params.W-1).*(g./params.W).^2);

%% calculate rhs of eom

ghat = Y(1:end-2); %fft of g
g = ifft(ghat);    %g
Phi = Y(end-1);
Psi = Y(end);

Kinv = 1./(1 + m*a./k);

f = psi_c2(g,Phi);

dghat = -0.5* Kinv.*(nu*(k.^2).*ghat + 1i*k.*ghat) + a*Kinv.*(fft(f));
dPhi = (psi_c(Phi) - Psi)/l_c;
dPsi = (Phi - gamma*sqrt(Psi))/(4*(B^2)*l_c);

dY = [dghat; dPhi; dPsi];

