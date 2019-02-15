% 306

function [npanels] = number_of_panels(Uinf,alpha,npanels,c,K,NACA,tol2,int2)

Cl = [0 1];
i = 2;
M = NACA(1);
P = NACA(2);
TT = NACA(3);

if exist('n_panels.mat','file') && exist('tol.mat','file') && exist('int.mat','file')
    load('n_panels.mat')
    load('tol.mat')
    load('int.mat')
else
    while abs(Cl(i)-Cl(i-1))>tol2
        npanels = npanels + int2;
        i = i+1;
        Cl(i) = const_vortex_code(Uinf,alpha,M,P,TT,c,npanels);
    end
    
end
if exist('n_panels.mat','file') && exist('tol.mat','file') && exist('int.mat','file')
    if tol ~= tol2 || int ~= int2
        npanels = 100;
        while abs(Cl(i)-Cl(i-1))>tol2
            npanels = npanels + int2;
            i = i+1;
            Cl(i) = const_vortex_code(Uinf,alpha,M,P,TT,c,npanels);
        end
    end
end
tol = tol2;
int = int2;
save('n_panels','npanels')
save('tol','tol')
save('int','int')
end