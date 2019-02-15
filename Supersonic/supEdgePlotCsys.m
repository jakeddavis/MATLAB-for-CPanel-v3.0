
% Jake Davis
% 01/01/2019

function supEdgePlotCsys(pnt1,pnt2,P)

% Redfine Stuff
x = P(1); y = P(2); z = P(3);
x1 = pnt1(1); y1 = pnt1(2);
x2 = pnt2(1); y2 = pnt2(2);
m  = (y2-y1) / (x2-x1);

vec = linspace(P(2),2,100);
for i = 1:length(vec)
    s1 = vec(i) - y1;
    s2 = vec(i) - y2;
    xm(i) = (x-x1) - (vec(i)-y1)/m;
    xm(i) = (x-x1) - (vec(i)-y1)/m;
    sm1  = xm(i) + s1/m;
    sm2  = xm(i) + s2/m;
    ym1(i)  = s1 - sm1/m;
    ym2(i)  = s2 - sm2/m;
    
    xmc(i)  = m*xm(i);
    ym1c(i) = m*ym1(i);
    ym2c(i) = m*ym2(i);
end

vec1 = linspace(x1,2,100);
for i = 1:length(vec1)
    yy1(i) = (vec1(i)-x1)/m + y1;
end
vec2 = linspace(x2,2,100);
for i = 1:length(vec2)
    yy2(i) = (vec2(i)-x2)/m + y2;
end
vec3 = linspace(y1,-2,100);
for i = 1:length(vec3)
    xx3(i) = m*(vec3(i)-y1) + x1;
end
vec4 = linspace(y2,-2,100);
for i = 1:length(vec4)
    xx4(i) = m*(vec4(i)-y2) + x2;
end

if abs(m) > 1
    figure
    hold on
    title('sup edge')
    plot(vec,xm,'k')
    plot(vec,ym1,'b')
    plot(vec,ym2,'r')
    axis equal
    grid on
    
    figure(2)
    hold on
    plot(vec1,yy1,'b')
    plot(vec2,yy2,'r')
    axis equal
    grid on
else
    figure
    hold on
    title('sub edge')
    plot(vec,xmc,'k')
    plot(vec,ym1c,'b')
    plot(vec,ym2c,'r')
    axis equal
    grid on
    
    figure(2)
    hold on
    plot(xx3,vec3,'b')
    plot(xx4,vec4,'r')    
    axis equal
    grid on
    
    figure(2)
    hold on
    plot(vec1,yy1,'x')
    plot(vec2,yy2,'o')
    axis equal
    grid on
end

end