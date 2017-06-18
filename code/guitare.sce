//Driver=("Rec");

clf();
clear;
//Definition des paramètres initiaux
SR = 44100;
B = 0.001;
f = 110;
TF = 1;
xo = 0.1;
co = 1;
rp = [0.3, 0.7];
loss = [100, 10; 1000, 8];

//pas de dicretisation
k=1/SR;

//déclaration des paramètres
gama=2*f;
kappa=(2*f*sqrt(B))/%pi;
N=50; // Il faut h >= sqrt(((gama^2)*(k^2)+sqrt((gama^4)*(k^4)+16*(kappa^2)*(k^2)))/2); 
h=1/N;

ksi1=((-gama^2)+sqrt((gama^4)+4*(kappa^2)*(loss(1,1)*2*%pi)^2))/(2*kappa^2);
ksi2=((-gama^2)+sqrt((gama^4)+4*(kappa^2)*(loss(2,1)*2*%pi)^2))/(2*kappa^2);

sigma1=(6*log(10)/(ksi2-ksi1))*((-1/loss(1,2))+1/loss(2,2));
sigma0=(6*log(10)/(ksi2-ksi1))*((ksi2/loss(1,2)-ksi1/loss(2,2)));

//création des matrices A, B et C
Dxxligne=zeros(1,N-1);
Dxxligne(1)=-2;
Dxxligne(2)=1;
Dxx=toeplitz(Dxxligne)*1/(h^2);

Dxxxx=Dxx * Dxx;
Id=eye(N-1,N-1);

A=(1+sigma0*k)*Id-sigma1*k*Dxx;
B=-2*Id-(gama^2)*(k^2)*Dxx+(kappa^2)*(k^2)*Dxxxx;
C=(1-sigma0*k)*Id+sigma1*k*Dxx;

//calcul de la suite U0
u0=zeros(N-1,1);
i=1;
while i*h <= xo
    u0(i)= i*h*co/xo
    i=i+1;
end

while i<=N-1
    u0(i)=(co/(xo-1))*i*h+co/(1-xo)
    i=i+1;
end
//Déclaration du vecteur stockant la position de la corde au niveau des micros
out=zeros(TF/k,2);
temps=zeros(TF/k,1);

//Calcul du premier itéré (U1)
inva=inv(A);
v=linspace(0,1,N+1)';
u1=-inva*(B*u0+C*u0);

//Remplissage vecteur de position
out(TF/k,1)=u0(floor(rp(1)/h));
out(TF/k,2)=u0(floor(rp(2)/h));
out(TF/k,1)=u1(floor(rp(1)/h));
out(2,2)=u1(floor(rp(2)/h));

//Mise à jour du vecteur temps
temps(1)=0;
temps(2)=k;

plot(v,[0;u0;0]);
xtitle("corde de guitare")
clf();
plot(v,[0;u1;0]);
xtitle("corde de guitare")

t=2*k;
c=0;   // compteur
im=0; // variable numérotant les images

//Calcul de la suit Un
while t<=TF
    temps(c+2)=t;
    u2=-inva*(B*u1+C*u0);
    u0=u1;
    u1=u2;
    t=t+k;
    
//Tracé toute les 30 images
    if modulo(c,30)==0
        drawlater;
        clf();
        subplot(2,1,1);
        plot(v,[0;u2;0]);
        xtitle("corde de guitare")
        a=gca();
        a.data_bounds=[0,-2 ; 1,2];
        subplot(2,2,3)
        plot(temps,out(:,1));
        xtitle("position corde micro 1")
        subplot(2,2,4)
        plot(temps,out(:,2));
        xtitle("position corde micro 2")
        drawnow;
//Enregistrement des images.gif
//        nom_image='image_'+string(im)+'.gif';
//        winnum=winsid();
//        xs2gif(winnum($),nom_image);
        im=im+1
 end
    //mise à jour du vecteur position
    out(t/k,1)=(u2(floor(rp(1)/h)) + u2(ceil(rp(1)/h)))/2;
    out(t/k,2)=(u2(floor(rp(2)/h)) + u2(ceil(rp(2)/h)))/2;
    c=c+1;

end
//enregistrement son.
playsnd(out(1 : size(out,1)), SR);
savewave("son.wav", out(1 : size(out, 1)), SR);


// Calcul de la transformée de Fourier
tf_s1=fft(out(:, 1));
plot((1:size(tf_s1, 1))/TF, tf_s1);



