function hrmwrite(hrm,fname)


ncmpn=length(hrm);

%create and begin inputing information in the hermite file
out=fopen(fname,'wt');

fprintf(out,' HERMITE INPUT FILE\n\n');
fprintf(out,' NUMBER OF COMPONENTS = %d\n',ncmpn);

for ii = 1:ncmpn
  
  c=hrm(ii);
  cmpname=c.name;
  group=c.group;
  typ=c.type;
  ncross=c.ncross;
  npt=c.npt;
  cmpn=c.pts;
  
  fprintf(out,'\n');
  fprintf(out,'%s \n',cmpname);
  fprintf(out,' GROUP NUMBER      = %d \n',group);
  fprintf(out,' TYPE              = %d\n',typ);
  fprintf(out,' CROSS SECTIONS    = %d \n',ncross);
  fprintf(out,' PTS/CROSS SECTION = %d \n',npt);
  
  
  for jj=1:ncross
    for kk=1:npt
      fprintf(out,'%10.5f  %10.5f  %10.5f\n',cmpn(jj,kk,1),cmpn(jj,kk,2),cmpn(jj,kk,3));
    end
  end
end

fclose(out);
