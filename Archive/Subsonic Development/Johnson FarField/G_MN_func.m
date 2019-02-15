
% Jake Davis
% 02/25/2018

function G_MN = G_MN_func(pnt1,pnt2,M,N,triGeom)

a = triGeom.a;
nu_xi = triGeom.nu_xi;
nu_eta = triGeom.nu_eta;

x1 = pnt1(1); y1 = pnt1(2);
x2 = pnt2(1); y2 = pnt2(2);

if abs(nu_eta) <= abs(nu_xi)
    if M == 1
        N = N + 1;
        E_MN0 = (x2^(M-1)*y2^(N-1)) - (x1^(M-1)*y1^(N-1));
        D_MN = E_MN0;
        G_MN = (1/(N*nu_xi)) * D_MN;
    else
        G_MN = -(nu_eta/nu_xi)*G_MN_func(pnt1,pnt2,M-1,N+1,triGeom)...
              + (a/nu_xi)*G_MN_func(pnt1,pnt2,M-1,N,triGeom);
    end
else
    if N == 1
        M = M + 1;
        E_MN0 = (x2^(M-1)*y2^(N-1)) - (x1^(M-1)*y1^(N-1));
        D_MN = E_MN0;
        G_MN = -(1/(M*nu_eta)) * D_MN;
    else
        G_MN = -(nu_xi/nu_eta)*G_MN_func(pnt1,pnt2,M+1,N-1,triGeom)...
              + (a/nu_eta)*G_MN_func(pnt1,pnt2,M,N-1,triGeom);
    end
end

end
