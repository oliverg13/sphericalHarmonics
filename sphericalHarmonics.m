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
