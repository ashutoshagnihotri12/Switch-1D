clc
close all
clear all
% 1D Danijela's approach
% This program uses adiabatic elimination approximation.
limiter = @newKorenLimiter;
[mu_e, D_e, n_0, l_0, t_0, E_0] = units(1.0);
L = 200;
m = 200;
h = L/m;
x = h*([1:m] - 1/2);
mu = 0.35;
gamma = 0.07;
phiLB = 17;  phiRB = 0 ;

% Initial condition(s) for electron and ion-density.

% Initial condition 1
% x0 = 0.9*L;
% sigma = 10;
% np = 1e-6*exp(-(x - x0).^2/sigma^2);
% NE = 0*np;
% x0 = 0.5*L;
% sigma = 10;
% np = 3e-3*exp(-((x-x0).^2)/sigma^2);
% NE = 0*np;

% Initial condition 2
% np = (0.5*L < x & x < 0.7*L)*1e-6;
% NE = 0*np;


% Initial condition 3
np = (0.5*L < x & x < 0.7*L)*1e-6;
NE = 0*np;

% Initial condition 4
% np = csvread('np.dat');
% NE = csvread('ne.dat');

counter = 0;
T = 0;
jE_int = 0*x;
WRT = 0;
% while(counter < 10*131)
while((T/mu)*t_0 <= 10000e-9)
%-------------------------------------------------------------------------%
np_t = np;
F_np = 0*x;
RK = 1;
while(RK <= 2)
    rho = np;
    rho(1) = rho(1) + (2*phiLB)/(h^2);                % corr. in source at LB
    rho(end) = rho(end) + (2*phiRB)/(h^2);            % corr. in source at RB
    A = spdiags(ones(m,1)*[1 -2 1],-1:1,m,m);
    A(1,1) = -3;    A(end,end) = -3;                  % Dirichlet b.c.
    A1D = A/(h^2);
    phi = -A1D\(rho');
    phi = phi';
    
    phiW = 2*phiLB - phi(1);                            % Left Bndy ghost cell
    phiE = 2*phiRB - phi(end);                          % Rght Bndy ghost cell
    phiL = [phiW phi(1:end-1)];                         % phi (j-1)
    phiR = [phi(2:end) phiE];                           % phi (j+1)
    
    EL = (phiL-phi)/h;
    ER = (phi-phiR)/h;
    E = (EL+ER)/2;
    %---------------------------------------------------------------------%
    dt = 0.5*(h/(eps+max(E)));
    
    npW = -np(1);   npWW = -np(2);
    npE = np(end);  npEE = np(end-1);
    
    npL = [npW np(1:end-1)];
    npLL = [npWW npW np(1:end-2)];
    npR = [np(2:end) npE];
    npRR = [np(3:end) npE npEE];
    
    [AL, AR, BL, BR] = StateValueEuler(np, npL, npLL, npR, npRR, limiter);
    GR = max(ER,0).*BL + min(ER,0).*BR;
    GL = max(EL,0).*AL + min(EL,0).*AR;
    
    if (min(E) < 0.0)
        break;
    end
    NE(end) = gamma*np(end);
    for j = m-1:-1:1
        NE(j) = (NE(j+1)*ER(j))/(EL(j) - h*abs(E(j))*exp(-1/abs(E(j))));
    end
    %     if (max(ne) <= 10^-30 || counter == 5000)
    %         csvwrite('ne_adia.dat',NE);
    %         csvwrite('np_adia.dat',np);
    %         break
    %     end
    Q = trapz(x,np-NE);
    F_np = - ((GR - GL)/h) + 1.0*NE.*abs(E).*exp(-1./abs(E)) + F_np;
    if(RK == 1)
        np = np + dt*F_np;
    end
    RK = RK+1;
end
% if(WRT == 0)
%   csvwrite('ne_init_adia.dat',mu*NE);
%   csvwrite('np_init_adia.dat',np_t);
%   WRT=WRT+1;
% end
np = np_t + 0.5*dt*F_np;
jE = mu_e*n_0*(E_0^2)*(mu*NE+mu*np).*(E.^2);           % ??
jE_int = jE*dt*(t_0/mu) + jE_int;

counter = counter + 1;

% F(counter) = getframe(gcf);
T = T + dt;
(T/mu)*t_0

subplot(2,2,1)
plot(x*l_0,NE*mu*n_0,'b-','LineWidth',2)
axis tight
grid on
title('n_{e}[cm^{-3}]')

subplot(2,2,2)
plot(x*l_0,np*n_0,'b-','LineWidth',2)
axis tight
grid on
title('n_{p}[cm^{-3}]')

% subplot(2,2,3)
% plot(x*l_0,E*E_0,'b-','LineWidth',2)
% axis tight
% grid on
% title('E[kV cm^{-1}]')

% subplot(2,2,3)
% plot(x*l_0,F_np/2,'b-','LineWidth',2)
% axis tight
% grid on

subplot(2,2,3)
semilogx((T/mu)*t_0 ,F_np(end),'b.','LineWidth',2)
axis tight
grid on
hold on

% subplot(2,2,4)
% plot(x*l_0,jE_int,'b-','LineWidth',2)
% axis tight
% grid on
% title('jE_{int}')

subplot(2,2,4)
plot(x*l_0,sign(np-mu*NE),'b.','LineWidth',2)
axis tight
grid on

% ax = axes;
% T1 = title(num2str((T/mu)*t_0));
% ax.Visible = 'off';
% T1.Visible = 'on';
% 
drawnow
end
% ff = figure;
% ff.Position = [1375 1072 863 187];
% subplot(2,2,1)
% semilogy(x*l_0,NE*mu*n_0+eps,'r-.','LineWidth',2)
% axis tight
% grid on
% title('n_{e}[cm^{-3}]')
% 
% subplot(2,2,2)
% semilogy(x*l_0,np*n_0+eps,'r-.','LineWidth',2)
% % plot(x*l_0,np*n_0+eps,'r-.','LineWidth',2)
% axis tight
% grid on
% title('n_{p}[cm^{-3}]')

% subplot(2,2,4)
% plot(x*l_0,jE_int,'r-.','LineWidth',2)
% axis tight
% grid on
% title('jE_{int}')

clc;
% ax = axes;
% T1 = title(num2str((T/mu)*t_0));
% ax.Visible = 'off';
% T1.Visible = 'on';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% video = VideoWriter('movieAdia_chk.avi','Motion JPEG AVI');
% open(video);
% writeVideo(video,F);
% close(video);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% subplot(2,2,1)
% plot(x*l_0,NE*mu*n_0,'b-','LineWidth',2)
% axis tight
% grid on
% title('n_{e}[cm^{-3}]')
% 
% subplot(2,2,2)
% plot(x*l_0,np*n_0,'b-','LineWidth',2)
% axis tight
% grid on
% title('n_{p}[cm^{-3}]')
% 
% subplot(2,2,3)
% plot(x*l_0,E*E_0,'b-','LineWidth',2)
% axis tight
% grid on
% title('E[kV cm^{-1}]')
% 
% subplot(2,2,4)
% plot(x*l_0,jE_int,'b-','LineWidth',2)
% axis tight
% grid on
% title('jE_{int}')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% np = np - ((GR - GL)/h)*dt + 1.0*NE.*abs(E).*exp(-1./abs(E))*dt;
% subplot(2,2,4)
% plot((T/mu)*t_0, trapz(x,np-NE*mu),'b.',(T/mu)*t_0, 0.0, 'r.')
% axis tight
% hold on
% title('(np-ne) vs t')
% np(70:130) = 1e-6;
% subplot(2,2,4)
% plot(x*l_0,jE,'b-','LineWidth',2)
% axis tight
% title('jE vs x')
 