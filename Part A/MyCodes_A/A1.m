%% 0.General
clear;clc;close all;

V = 1; % V[Volt]
eo = 8.854187817e-12; % absolute dielectric permittivity of a vacuum
er = 1; % relative dielectric permittivity of the air
dc = er*eo; % 'd'ielectric 'c'onstant
Z0 = 50;

%% 1.Geometry
ra = 0.7605e-3;
rb = 1.75e-3;

gd = [1 1; 0 0; 0 0; ra rb]; % geometry description matrix
ns = [82 82; 49 50];
sf = 'R2-R1';

%% 2.Mesh
d1 = decsg(gd,sf,ns);
[p, e, t] = initmesh(d1);
Nref = 3; % number of "refinements"
for i=1:Nref
    [p, e, t] = refinemesh(d1,p,e,t);
end
     
figure, pdeplot(p,e,t);
xlabel('x[m]'), ylabel('y[m]')
title('Plane Discretization')

%% 3."Re-number" the unknown nodes

Nn = size(p,2); % number of 'n'odes
Nd = size(e,2); % number of e'd'ges
Ne = size(t,2); % number of 'e'lements

% "1" for the unknown nodes, "0" for the known nodes
node_id = ones(Nn,1);
% initialize the final solution 
X0 = zeros(Nn,1);

for id = 1:Nd
    % check if the id-th edge belongs to an interface
    if ( (e(6,id)==0) || (e(7,id) == 0) )
        % mark both two nodes of the edge as known
        node_id(e(1,id)) = 0;
        node_id(e(2,id)) = 0;
        % check which edges belong to the ra-circle -> +V(volt)
        % sqrt(x^2 + y^2) < (a+b)/2
        if (p(1,e(1,id))^2 + p(2,e(1,id))^2 < (ra+rb)^2/4)
            X0(e(1,id)) = V;
            X0(e(2,id)) = V;
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
S = spalloc(Nn,Nn,7*Nn);
Sff = spalloc(Nf,Nf,7*Nf);
B = zeros(Nf,1);

for ie = 1:Ne % ie-th element
    n(1:3) = t(1:3,ie); % 3-nodes 
    rg = t(4,ie); % region of the ie-th triangle
    x(1:3) = p(1,n(1:3)); y(1:3) = p(2,n(1:3)); % nodes' coordinates (x,y)
    De = det([ones(3,1) x' y']);
    Ae = abs(De/2); % area of the triangle
    
    % z(i) = a(i) + b(i)*x + c(i)*y, i={1,2,3}
    b = zeros(3,1);
    c = zeros(3,1);
    for k = 0:2
        xk = circshift(x,-k);
        yk = circshift(y,-k);
        b(k+1) = (yk(2)-yk(3))/De;
        c(k+1) = (xk(3)-xk(2))/De;
    end
    
    % S_ff & B
    Se = zeros(3);
    for i=1:3
        for j=1:3
            Se(i,j) = dc*(b(i)*b(j)+c(i)*c(j))*Ae;
            S(n(i),n(j)) = S(n(i),n(j)) + Se(i,j);
            if (node_id(n(i))~=0) 
                if (node_id(n(j))~=0) % i,j:unknown
                    Sff(index(n(i)),index(n(j))) = Sff(index(n(i)),index(n(j))) + Se(i,j);
                else % i:unknown  j:known
                    B(index(n(i))) = B(index(n(i))) - Se(i,j)*X0(n(j));
                end
            end
        end
    end  
end

% direct solver
tic
X = Sff\B;
toc

% iterative solver: "biconjugate gradient" -> Ax = b
% [x,flag] = bicg(A,b,tol,maxit)
% If tol is [], then bicg uses the default, tol=1e-6
% maxit -> maximum number of iterations
% flag=0 -> bicg converged to the desired tolerance tol within maxit iterations
% flag=1 -> bicg iterated maxit times but did not converge.
tic
[X1, flag] = bicg(Sff,B,[],170); % 170 iterations for Nref=3 in order to converge (tol=1e-6)
toc

% complete the X0 matrix with the computed values of the unknown nodes
for i = 1:Nn
    if(index(i)>0)
        X0(i) = X(index(i));
    end
end

figure, pdeplot(p,e,t,'XYdata',X0)
hold on;
[ux,uy] = pdegrad(p,t,-X0);
pdeplot(p,t,'FlowData',[ux; uy])
xlabel('x[m]'), ylabel('y[m]')
legend('\bf \phi','\bf E')

En = X0'*S*X0/2; % Energy
C = X0'*S*X0/V^2; % FEM capacity
Cth = 120*pi*eo*sqrt(er)/Z0; % theoretical capacity

%% phi -> theoretical approach
K = -V/log(rb/ra);
p = linspace(ra, rb, 200);
v = K*log(p./rb);

figure, plot(1000*p,v)
xlim(1000*[ra,rb])
xlabel('\rho [mm]'), ylabel('\phi [V]')
title('\phi(\rho) = -Vln(\rho/b)/ln(b/a)')

line(1000*[p(100),p(100)],[0,v(100)],'Color','red','LineStyle','--')
line(1000*[0,p(100)],[v(100),v(100)],'Color','red','LineStyle','--')
text(1000*p(100),v(100),' \leftarrow \phi=0.4, at \rho=(a+b)/2')