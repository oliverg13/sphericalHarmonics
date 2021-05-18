%% Definition of coordinates
% a = polar angle \theta
% b = azimutal angle \phi
k=100; %resolution for harmonic
[b,a]=meshgrid(linspace(0,2*pi,k),linspace(0,pi,k));
da=pi/(k-1); db=2*pi/(k-1);

%% Definicion de funcion g(theta,phi)
g=sin(4*a+b)+exp(-( (a-pi/2) / (pi/6) ).^2   ).*exp(-( (b-pi/2) / (pi/6) ).^2     );
%g=ylm(2,-1,k,a,b)+ylm(3,1,k,a,b)+ylm(0,0,k,a,b);

%% Construccion de coeficientes glm
n=10;
glm=zeros(n,2*n-1); 
for l=0:(n-1)
    for m=(-l):l
        glm(l+1,m+n)=sum(sum(g.*conj(sphericalHarmonics(l,m,k)).*sin(a)))*da*db;
    end
end

%% Graficas para comparar
[X1,Y1,Z1] = sph2cart(a, b, abs(real(g)));
figure(1)
hold on
surf(X1,Y1,Z1)
title('$g(\theta,\phi)$','interpreter','latex')
hold off

figure(2)
hold on
[X2,Y2,Z2] = sph2cart(a, b, abs(real(gr)));
surf(X2,Y2,Z2)
title('$gr(\theta,\phi)$','interpreter','latex')
hold off

figure(3)
hold on
[X3,Y3,Z3] = sph2cart(a, b, abs(imag(gr)));
surf(X3,Y3,Z3)
title('$Im(gr(\theta,\phi))$','interpreter','latex')
hold off

%% Analisis de errores
t=abs((gr-g)./g);
ne=sum(sum(t<0.01));
pnot=ne/(k^2);

%% Representacion sobre una esfera
[X4,Y4,Z4] = sph2cart(a, b, ones(k,k));
figure(4)
hold on
surf(X4,Y4,Z4,real(g))
title('$g(\theta,\phi)$','interpreter','latex')
view(45,45)
hold off

[X5,Y5,Z5] = sph2cart(a, b, ones(k,k));
figure(5)
hold on
surf(X5,Y5,Z5,real(gr))
title('$gr(\theta,\phi)$','interpreter','latex')
view(45,45)
hold off
