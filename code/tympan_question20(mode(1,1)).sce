clf;
clear;

a = 1;
N_theta = 100;
theta = [1:N_theta+1]
d_theta = 2*%pi/(N_theta-1);
N_eta = 20;
eta = [0:N_eta-1];
d_eta = a/((N_eta-1));
CFL = 0.5;
c = 1;
lambda = 8.65373;

d_tau = CFL*d_eta*d_theta;
N = 20000;
T_simu = N;
Vect_Tmps = linspace(0, 1, N+1)*d_tau;
erreur = zeros(1, N+1);

//Position initiale
omega0 = zeros(N_eta, N_theta + 1)

pos = besselj(1, 3.83171*eta*d_eta);

for i = 1:(N_theta+1)
    omega0(:, i) = 0.5*cos(i*d_theta)*pos';
end
omega0(N_eta, :) = zeros(1, N_theta+1); //condition limite au bord de la membrane

omega_exact=omega0;

//Calcul des positions au cours du temps
x = eta'*d_eta*a*cos(theta*d_theta);
y = eta'*d_eta*a*sin(theta*d_theta);
[xf, yf, zf] = nf3d(x', y', omega0');
//affichage de la courbe


plot3d1(xf, yf, zf);
plot3d1(xf, yf, zf);
f=gcf();
f.color_map=jetcolormap(64);

//Tracé des positions successives
T_ligne = zeros(1, N_theta+1);
T_ligne(1) = -2;
T_ligne(2) = 1;
T = toeplitz(T_ligne);
for i = 1:N_theta+1
    T(i, :) = (1/i^2)*T(i,:);
end

T3_ligne = zeros(1, N_eta);
T3_ligne(1) = -2;
T3_ligne(2) = 1;
T3 = toeplitz(T3_ligne);

T2_ligne = ones(1, N_eta-1);
T2_sup = diag(T2_ligne, 1);
T2_inf = diag(T2_ligne*(-1), -1);
T2 = T2_sup + T2_inf;

for i = 1:N_eta
    T2(i, :) = (1/i)*T2(i,:);
end

omega1 = omega0;
t = 1;
im = 2;

while t <= T_simu
    tau=t*d_tau
    omega2=2*omega1- omega0 + d_tau^2 *((1/d_eta^2)*T3*omega1 + (1/(d_eta*d_theta)^2)*(T*omega1')' + (1/(d_theta*d_eta))*T2*omega1);

    //Calcul de omega2(1,1)
    S=0;
    for j = 1:N_theta
        S=S+omega1(2,j)-omega1(1,1);
    end
    omega2(1,1)=2*omega1(1,1)-omega0(1,1)+((4*d_tau^2)/(d_eta^2*N_theta))*S;
    //calcul de la première ligne
    for j = 2:N_theta
        omega2(1,j)=omega2(1,1);
    end

    //calcul première colonne
    for i = 2:N_eta-1
        omega2(i,1)=2*omega1(i,1)-omega0(i,1)+d_tau^2*((1/d_eta^2)*(omega1(i+1,1)-2*omega1(i,1)+omega1(i-1,1))+(1/(i*d_eta*d_theta))^2*(omega1(i,2)-2*omega1(i,1)+omega1(i,N_theta))+(1/(2*i*d_eta*d_theta))*(omega1(i+1,1)-omega1(i-1,1)));
    end
    omega2(:,N_theta+1)=omega2(:,1);
    omega2(N_eta,:)=zeros(1,N_theta+1);
    omega0 = omega1;
    omega1 = omega2;

    drawlater;
    clf();
    [xf, yf, zf] = nf3d(x', y', omega2');
    plot3d1(xf, yf, zf);
    a=gca();
    a.data_bounds=[-1,-1,-0.5;1,1,0.5];
    t = t+1;
    drawnow;
end
