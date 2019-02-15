function hrmplot(hrm)

hold on;

ncmpn=length(hrm);
for ii=1:ncmpn
  
  c=hrm(ii);
  cmpn=c.pts;
  
  surf(cmpn(:,:,1),cmpn(:,:,2),cmpn(:,:,3),ii*ones([size(cmpn,1) size(cmpn,2)]));
  %   surf(cmpn(:,:,1),cmpn(:,:,2),cmpn(:,:,3));
end
view(3)
hold off
axis equal
axis off

end
