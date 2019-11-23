%% 0.General
clear; clc; 
V = 100; % V[Volt]
eo = 8.85418e-12; % absolute dielectric permittivity of a vacuum
er = 2.2; % relative permittivity

%% 1.Geometry
w = 0.03;
h = 2e-3;
d = 0.01;

gd1 = [w/2 w/2 -w/2 -w/2 d/2 d/2+h d/2+h d/2]';
gd2 = -gd1;
gd3 = 2.5*w*[1 1 -1 -1 -1 1 1 -1]';

gd = [3 3 3; 4 4 4; gd1 gd2 gd3]; % geometry description matrix
ns = [82 82 82; 49 50 51];
sf = 'R3-R1-R2';

%% 2.Mesh
dl = decsg(gd,sf,ns);

[p, e, t] = initmesh(dl);
Ne = size(t,2); % number of 'e'lements

Nref = 1; % number of "refinements"
for i=1:Nref
    [p, e, t] = refinemesh(dl,p,e,t);
end

% refine the triangles that belong to a (w/10)x(d)-rectangle
% center: (w/2,0) and (-w/2,0)
j = 1;
% re = zeros(1,400);
for ie = 1:Ne
    for en = 1:3
        if ( abs(p(1,t(en,ie)))>=w/2-w/20 && abs(p(1,t(en,ie)))<=w/2+w/20 && abs(p(2,t(en,ie)))<=d/2 )
            re(j) = ie;
            j = j + 1;
        end
    end
end
      
[p, e, t] = refinemesh(dl,p,e,t,re');

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
        
        % The nodes of the edge are known only if they belong to the
        % "internal" interface, (i.e. capacitor). Thus, only the nodes
        % that belong to an interface, inside a (WxD)-rectangle (W>D), are KNOWN.
        % *******SOS*******
        % We assume that the "external" nodes are UNKNOWN, because of 
        % the Neumann boundary conditions applied to that interface.
        
        if ( abs(p(1,e(1,id)))<w && abs(p(2,e(1,id)))<d )
            % mark both two nodes of the edge as known
            node_id(e(1,id)) = 0;
            node_id(e(2,id)) = 0;
            % check which conductor we refer to
            if ( p(2,e(1,id))>0 ) % conductor with y>0
                X0(e(1,id)) = V/2;
                X0(e(2,id)) = V/2;
            else
                X0(e(1,id)) = -V/2;
                X0(e(2,id)) = -V/2;
            end
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
Nf = sum; % number of unknown nodes
S = spalloc(Nn,Nn,7*Nn);
Sff = spalloc(Nf,Nf,7*Nf);
B = zeros(Nf,1);

for ie = 1:Ne % ie-th element
    
    n(1:3) = t(1:3,ie); % 3-nodes 
    rg = t(4,ie); % region of the ie-th triangle
    x(1:3) = p(1,n(1:3)); y(1:3) = p(2,n(1:3)); % nodes' coordinates (x,y)
    if ( max(abs(x))<=w/2 && max(abs(y))<=d/2 )
        dec = er*eo; % dielectric of the capacitor
    else
        dec = eo; % dialectric outside the capacitor
    end
%     dec = eo; 
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
            Se(i,j) = dec*(b(i)*b(j)+c(i)*c(j))*Ae;
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
X = Sff\B;

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
C = X0'*S*X0/V^2; % Capacity