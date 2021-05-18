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

%% Reconstruccion
gr=zeros(k,k); % prelocalizacion de funcion reconstruida
%qq=1; %contador
%figure(200)
hold on
for l=0:(n-1)
    for m=(-l):l
        gr=gr+glm(l+1,m+n)*sphericalHarmonics(l,m,k);
    end  
    %  subplot(2,3,qq)
    %    [Xq,Yq,Zq] = sph2cart(a, b, ones(k,k));
    %    surf(Xq,Yq,Zq,real(gr))
    %    subtitle(["l=" num2str(l) ", all m"])
    %    qq=qq+1;
end
hold off


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

%% sphericalHarmonics function
%% Input data
% L... Azimuthal quantum number
% M... Magnetic quantum number
% k... resolution of spacement for azimutal and polar angles
%%%%%%%%%%%%%%

%% Output
% Y... Spherical harmonic matrix representation k times k size
%%%%%%%%%%%%%%

function Y = sphericalHarmonics(L,M,k)
    %% Verification of data input
    if rem(L,1)~=0 && L<0
        error('L is not a postive integer number')
    elseif rem(M,1)~=0 
            error('M is not an integer')
    elseif abs(M)>=L+1
            error('abs(M) cannot be larger than L')
    elseif nargin==1 || nargin==2
        k=100;
    elseif nargin==3 && rem(k,1)~=0 && k<0
        error('l is not a postive interger number')
    end
    
    [b,a]=meshgrid(linspace(0,2*pi,k),linspace(0,pi,k)); % Definition of angular space
    
    m=1; %introduced local variable to work with negative M values
    if M<0
        m=M;
        M=-M;
    end
    n = sqrt((2*L+1)/( 4*pi*factorial(2*M)*nchoosek(L+M,L-M) ));
    plm=zeros(k,k); %prelocation of variable of legendre polinomial
    for c=1:k
        for cc=1:k
            pl=legendre(L,cos(a(c,cc)));
            plm(c,cc)=pl(M+1);
        end
    end
    Y=n*plm.*exp(1i*M*b);
    
    if m<0 
        Y=conj(Y)*(-1)^m;
    end
end

