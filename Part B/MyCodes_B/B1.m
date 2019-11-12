% 1.Geometry
clear;clc
r = 0.01;
gd = [1 0 0 r]'; % geometry description matrix

% 2.Mesh
d1 = decsg(gd);
[p, e, t] = initmesh(d1);
Nref = 2; % number of "refinements"
for Nm=1:Nref
    [p, e, t] = refinemesh(d1,p,e,t);
end
     
load('Ez_TM.mat'), load('fc_TM')
load('Hz_TE.mat'), load('fc_TE')

% TE
for Nm = 1:length(fc_TE)
    subplot(2,3,Nm)
    pdeplot(p,e,t,'XYdata',Hz_TE(:,Nm))
    title(['f_c = ', num2str(round(fc_TE(Nm),2)), 'GHz'])
    colormap('jet'); axis equal; axis tight;
end
suptitle('TE modes')

% TM
figure
for Nm = 1:length(fc_TM)
    subplot(2,3,Nm)
    pdeplot(p,e,t,'XYdata',Ez_TM(:,Nm))
    title(['f_c = ', num2str(round(fc_TM(Nm),2)), 'GHz'])
    colormap('jet'); axis equal; axis tight;
end
suptitle('TM modes')

% First 9 modes & fc_error
str = ["TE_{11}", "TM_{01}", "TE_{21}", "TM_{11}", "TE_{01}",... 
       "TE_{31}", "TM_{21}", "TE_{41}", "TE_{12}"];

% TEnm -> p'nm -> J'n(p'nm)=0
% TMnm -> pnm -> Jn(pnm)=0

pth = [1.8412 2.4048 3.0542 3.8317 3.8317 4.2012 5.1356 5.3175  5.3314];
fct = 0.3*pth/(2*pi*r); % fc theoretical (GHz)

i = 1; % TM
j = 1; % TE
fcc = zeros(1,9);
figure
for Nm = 1:9
     if ( fc_TM(i)<fc_TE(j) )
        subplot(3,3,Nm)
        pdeplot(p,e,t,'XYdata',Ez_TM(:,i))
        colormap('jet'); axis equal; axis tight; 
        fcc(Nm) = fc_TM(i);
        i = i+1;
     else
        subplot(3,3,Nm)
        pdeplot(p,e,t,'XYdata',Hz_TE(:,j))
        colormap('jet'); axis equal; axis tight;
        fcc(Nm) = fc_TE(j);
        j = j+1;
     end
     title(str(Nm))
     
end

fc_error = round(100*abs(fcc-fct)./fct,6);
