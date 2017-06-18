clf;
clear;

a = 1;
N_theta = 80;
theta = [1:N_theta+1]
d_theta = 2*%pi/(N_theta-1);
N_eta = 40;
eta = [0:N_eta-1];
d_eta = a/((N_eta-1));
CFL = 0.5;
c = 1;
lambda0 = 2.40483;
lambda1 = 3.83171;

d_tau = CFL*d_eta*d_theta;
N = 2000;
T_simu = N;
Vect_Tmps = linspace(0, 1, N+1)*d_tau;
erreur = zeros(1, N+1);

//Position initiale
omega0 = zeros(N_eta, N_theta + 1)

J0 = besselj(0, lambda0*eta*d_eta);
J1 = besselj(1, lambda1*eta*d_eta);

for i = 1:(N_theta+1)
    omega0(:, i) = J0' + J1'*cos(i*d_theta);
end
omega0(N_eta, :) = zeros(1, N_theta+1); //condition limite au bord de la membrane

omega_exact=omega0;

//Calcul des positions au cours du temps
x = eta'*d_eta*a*cos(theta*d_theta);
y = eta'*d_eta*a*sin(theta*d_theta);
[xf, yf, zf] = nf3d(x', y', omega0');
//affichage de la courbe

subplot(1,2,1);
plot3d1(xf, yf, zf);
subplot(1,2,2);
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

erreur(1) = 0;

while t <= 50000000*d_tau
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
    subplot(1,2,1);
    [xf, yf, zf] = nf3d(x', y', omega2');
    plot3d1(xf, yf, zf);
    a=gca();
    a.data_bounds=[-1,-1,-2;1,1,2];
  
    for i = 1:(N_theta+1)
        omega_exact(:, i) = cos(lambda0*tau)*J0' + cos(i*d_theta)*cos(lambda1*tau)*J1';
    end
    subplot(1,2,2);

    [xf, yf, zf] = nf3d(x', y', omega_exact');
    plot3d1(xf, yf, zf);
    f=gcf();
    f.color_map=jetcolormap(64);
    a=gca();
    a.data_bounds=[-1,-1,-2;1,1,2];

//    subplot(2,1,2);
//    erreur=(cos(lambda*t)*omega_exact-omega2);
//    plot(eta,erreur(N_eta,N_theta+1));
//    a=gca();
//    a.data_bounds=[-1,-1,-1;1,1,1];
//
    drawnow;
//
    //Calcul de l'erreur en eta=0
    erreur(t+1) = 100*(omega2(1,1)-omega_exact(1,1))/omega_exact(1,1);
//    nom_image='image_'+string(im)+'.gif';
//    winnum=winsid();
//    xs2gif(winnum($),nom_image);
//    im = im+1;
    t = t+1;
end

//clf;
//plot(Vect_Tmps, erreur);
