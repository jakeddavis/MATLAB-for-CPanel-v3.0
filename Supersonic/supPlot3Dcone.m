
% Jake Davis
% 01/20/2019

function supPlot3Dcone(min,max,P,M,varargin)

mu = asin(1/M);
len = max - min;
r = len * tan(mu);
n = 100;

[z, y, x] = cylinder([0 r],n);
if nargin > 4
    if varargin{1} == 1
        X = x(:,[1:26,76:end]);
        Y = y(:,[1:26,76:end]);
        Z = z(:,[1:26,76:end]);
    else
        X = x(:,26:76);
        Y = y(:,26:76);
        Z = z(:,26:76);
    end
else
    X = x;
    Y = y;
    Z = z;
end

X = -len*X + P(1);
Y = Y + P(2);
Z = Z + P(3);

surf(X,Y,Z,'FaceColor',[.5 .5 .5],'EdgeColor','none','FaceAlpha',0.5)

end