%% 0.General
clear;clc

eo = 8.8541878128e-12;
mo = 4*pi*1e-7;

lamda = 1; % [m]
w = lamda;
a = 2.5*lamda;
d = lamda;

%% 1.Geometry
gd1 = [3 4 a+w a+w -(a+w) -(a+w) a+w -(a+w) -(a+w) a+w]'; % kentro
gd2 = [3 4 a+w a+w -(a+w) -(a+w) a+w+d a+w a+w a+w+d]'; % panw
gd3 = [3 4 a+w a+w -(a+w) -(a+w) -(a+w) -(a+w+d) -(a+w+d) -(a+w)]';
gd4 = [3 4 -(a+w) -(a+w) -(a+w+d) -(a+w+d) a+w -(a+w) -(a+w) a+w]';
gd5 = [3 4 a+w+d a+w+d a+w a+w a+w -(a+w) -(a+w) a+w]';
gd6 = [3 4 -(a+w) -(a+w) -(a+w+d) -(a+w+d) a+w+d a+w a+w a+w+d]';
gd7 = [3 4 a+w+d a+w+d a+w a+w a+w+d a+w a+w a+w+d]';
gd8 = [3 4 a+w+d a+w+d a+w a+w -(a+w) -(a+w+d) -(a+w+d) -(a+w)]';
gd9 = [3 4 -(a+w) -(a+w) -(a+w+d) -(a+w+d) -(a+w) -(a+w+d) -(a+w+d) -(a+w)]';
gd10 = [1 0 0 a 0 0 0 0 0 0]';

gd = [gd1 gd2 gd3 gd4 gd5 gd6 gd7 gd8 gd9 gd10];
sf = '(R1+R2+R3+R4+R5+R6+R7+R8+R9)-Rc'; % Rc -> R circle
ns = [82*ones(1,10) ; [49:57 99] ];

[d1,bt] = decsg(gd,sf,ns);
% pdegplot(d1,'SubdomainLabels','on');

% matlab region numbers
matre = zeros(1,size(bt,2)-1);
for k = 1:length(matre)
    matre(k) = find(bt(:,k)==1);
end

%% 2.Mesh
[p, e, t] = initmesh(d1);
Nref = 4; % number of "refinements"
for k=1:Nref
    [p, e, t] = refinemesh(d1,p,e,t);
end

% figure, pdeplot(p,e,t);
% xlabel('x[m]'), ylabel('y[m]')
% title('Plane Discretization')

Nn = size(p,2); % number of 'n'odes
Nd = size(e,2); % number of e'd'ges
Ne = size(t,2); % number of 'e'lements

% adjust region numeration to desired numbers
for ie = 1:Ne
    t(4,ie) = find(matre==t(4,ie));
end

R = 1e-6; theta = 0;
beta = -lamda*log(R)/(4*pi*d*cosd(theta));
alpha = 1-1i*beta;

Lx = [1/alpha; alpha; alpha];  Ly = [alpha; 1/alpha; alpha];  Lxy = [1; 1; alpha^2];

M = zeros(3,length(matre));
M(:,1) = ones(3,1);
% e_zz
M(1,2:3) = Ly(3)*ones(1,2);
M(1,4:5) = Lx(3)*ones(1,2);
M(1,6:end) = Lxy(3)*ones(1,length(matre)-5);
% m_xx, m_yy
M(2:3,2:3) = Ly(1:2)*ones(1,2);
M(2:3,4:5) = Lx(1:2)*ones(1,2);
M(2:3,6:end) = Lxy(1:2)*ones(1,length(matre)-5);

M(1,1:end) = eo*M(1,1:end);
M(2:3,1:end) = mo*M(2:3,1:end);

% incident wave
E0 = 1; % [V/m]
Ei = zeros(Nn,1);
for in = 1:Nn
    if ( p(2,in)<a && p(2,in)>-a && p(1,in)>0 )
        Ei(in) = 0;
    else
        Ei(in) = E0*exp(-1i*2*pi*p(1,in)/lamda);
    end
end

% pdeplot(p,e,t,'XYdata',abs(Ei))

% "1" for the unknown nodes, "0" for the known nodes
node_id = ones(Nn,1);
% initialize the final solution, X0 -> Xs (scattering)
X0 = zeros(Nn,1);

for id = 1:Nd
    % check if the id-th edge belongs to an interface
    if ( (e(6,id)==0) || (e(7,id) == 0) )
        % check which edges belong to the a-circle -> Escat = -Ei
        % sqrt(x^2 + y^2) < 3a/2
        if ( sqrt( p(1,e(1,id))^2 + p(2,e(1,id))^2 ) < 3*a/2 )
             % mark both two nodes of the edge as known
            node_id(e(1,id)) = 0;
            node_id(e(2,id)) = 0;
            X0(e(1,id)) = -Ei(e(1,id));
            X0(e(2,id)) = -Ei(e(2,id));
        end
    end
end

% "re-number" the unknown nodes
index = zeros(Nn,1);
sum = 0;
for in = 1:Nn % in-th node
    if(node_id(in)>0)
        sum = sum + 1;
        index(in) = sum;
    end
end

%% 4.Main

Nf = sum; % total number of the unknown nodes
Af = spalloc(Nf,Nf,7*Nf);
Bf = zeros(Nf,1);

for ie = 1:Ne % ie-th element
    n(1:3) = t(1:3,ie); % 3-nodes 
    rg = t(4,ie); % region of the ie-th triangle
    x(1:3) = p(1,n(1:3)); y(1:3) = p(2,n(1:3)); % nodes' coordinates (x,y)
    De = det([ones(3,1) x' y']);
    Ae = abs(De/2); % area of the triangle
    
    % z(i) = a(i) + b(i)*x + c(i)*y, i={1,2,3}
    b = zeros(3,1);  c = zeros(3,1);
    b(1) = (y(2)-y(3))/De; c(1) = (x(3)-x(2))/De;
    b(2) = (y(3)-y(1))/De; c(2) = (x(1)-x(3))/De;
    b(3) = (y(1)-y(2))/De; c(3) = (x(2)-x(1))/De;
     
    % Af & Bf
    Se = zeros(3);  Pe = zeros(3);
    Te = (Ae/12)*ones(3); 
    Te = Te - diag(diag(Te)) + diag((Ae/6)*ones(1,3));
    Te = M(1,rg)*Te;
    for i=1:3
        for j=1:3
            Se(i,j) = ( b(i)*(1/M(3,rg))*b(j) + c(i)*(1/M(2,rg))*c(j) )*Ae;
            Pe(i,j) = Se(i,j) - (2*pi*3e8/lamda)^2*Te(i,j);
            if (node_id(n(i))~=0) 
                if (node_id(n(j))~=0) % i,j:unknown
                    Af(index(n(i)),index(n(j))) = Af(index(n(i)),index(n(j))) + Pe(i,j);
                else % i:unknown  j:known
                    Bf(index(n(i))) = Bf(index(n(i))) - Pe(i,j)*X0(n(j));
                end
            end
        end
    end  
end

% direct solver
X = Af\Bf;

% complete the X0 matrix with the computed values of the unknown nodes
for i = 1:Nn
    if(index(i)>0)
        X0(i) = X(index(i));
    end
end

figure, pdeplot(p,e,t,'XYdata',abs(Ei+X0))
axis equal; axis tight;
colormap('jet');
xlim([-a-w a+w]), ylim([-a-w a+w])
xlabel('x[m]'), ylabel('y[m]'), legend('E_{tot}')
title(['E_{tot} = E_s + E_i, for  a = ', num2str(lamda/a), '\lambda', ', d = \lambda'])
