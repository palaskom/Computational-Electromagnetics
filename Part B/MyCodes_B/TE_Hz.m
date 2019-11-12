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

% 4.Main

S = spalloc(Nn,Nn,7*Nn);
T = spalloc(Nn,Nn,7*Nn);

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
    
    % S_ff=S & T
    Se = zeros(3);  
    Te = (Ae/12)*ones(3); 
    Te = Te - diag(diag(Te)) + diag((Ae/6)*ones(1,3));
    for i=1:3
        for j=1:3
            Se(i,j) = (b(i)*b(j)+c(i)*c(j))*Ae;
            S(n(i),n(j)) = S(n(i),n(j)) + Se(i,j); 
            T(n(i),n(j)) = T(n(i),n(j)) + Te(i,j);
        end
    end  
end

[V,D] = eigs(S,T,13,'smallestabs');
fc = 3e8*sqrt(diag(D))/2/pi;
pnm = sqrt(diag(D))*r;

Nm = 0;
for i = 2:length(fc)
    Nm = Nm+1;
    X0 = V(:,i);      
    subplot(4,3,Nm)
    pdeplot(p,e,t,'XYdata',X0)
    title(['f_c = ', num2str(round(fc(i)*1e-9,2)), 'GHz'])
    colormap('jet'); axis equal; axis tight;
end
suptitle('TE modes')

Hz_TE = [V(:,3) V(:,4) V(:,6) V(:,7) V(:,9) V(:,12)];
fc_TE = [fc(3) fc(4) fc(6) fc(7) fc(9) fc(12)]*1e-9; % GHz

save('Hz_TE.mat','Hz_TE')
save('fc_TE','fc_TE')    
