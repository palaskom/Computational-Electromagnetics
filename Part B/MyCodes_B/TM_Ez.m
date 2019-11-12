% 1.Geometry
clear;clc
r = 0.01;
gd = [1 0 0 r]'; % geometry description matrix

% 2.Mesh
d1 = decsg(gd);
[p, e, t] = initmesh(d1);
Nref = 3; % number of "refinements"
for i=1:Nref
    [p, e, t] = refinemesh(d1,p,e,t);
end

% 3."Re-number" the unknown nodes

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

% 4.Main

Nf = sum; % total number of the unknown nodes
Sff = spalloc(Nf,Nf,7*Nf);
Tff = spalloc(Nf,Nf,7*Nf);

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
    Te = (Ae/12)*ones(3); 
    Te = Te - diag(diag(Te)) + diag((Ae/6)*ones(1,3));
    for i=1:3
        for j=1:3
            Se(i,j) = (b(i)*b(j)+c(i)*c(j))*Ae;
            if (node_id(n(i))~=0) 
                if (node_id(n(j))~=0) % i,j:unknown
                    Sff(index(n(i)),index(n(j))) = Sff(index(n(i)),index(n(j))) + Se(i,j);
                    Tff(index(n(i)),index(n(j))) = Tff(index(n(i)),index(n(j))) + Te(i,j);
                end
            end
        end
    end  
end

[V,D] = eigs(Sff,Tff,12,'smallestabs');
fc = 3e8*sqrt(diag(D))/2/pi;
pnm = sqrt(diag(D))*r;

X1 = zeros(Nn,12);
Nm = 0; % number of mode
for i = 1:length(fc)
    Nm = Nm+1;
    X = V(:,i);
    % complete the X0 matrix with the computed values of the unknown nodes
    for k = 1:Nn
        if(index(k)>0)
            X0(k) = X(index(k));
        end
    end  
    X1(:,i) = X0;
    subplot(4,3,Nm)
    pdeplot(p,e,t,'XYdata',X0)
    title(['f_c = ', num2str(round(fc(i)*1e-9,2)), 'GHz'])
    colormap('jet'); axis equal; axis tight;
              
end
suptitle('TM modes')

Ez_TM = [X1(:,1) X1(:,3) X1(:,4) X1(:,6) X1(:,7) X1(:,10)];
fc_TM = [fc(1) fc(3) fc(4) fc(6) fc(7) fc(10)]*1e-9; % GHz
save('Ez_TM.mat','Ez_TM')
save('fc_TM','fc_TM')
    