function [hrm, ncmpn] = hrmread(fname)

fp=fopen(fname,'r');
fposchar(fp,'=');

ncmpn=fscanf(fp,'%i',[1 1]);

hrm(ncmpn) = struct('name',[],'group',[],'type',[],'ncross',[],'npt',[],'pts',[]);

for ii=1:ncmpn
  
  c = struct('name',[],'group',[],'type',[],'ncross',[],'npt',[],'pts',[]);
  
  fgetl(fp); % Skip past end of last line read.
  fgetl(fp); % Skip blank line.
  c.name=strtrim(fgetl(fp));
  
  fposchar(fp,'=');
  c.group=fscanf(fp,'%i',[1 1]);
  
  fposchar(fp,'=');
  c.type=fscanf(fp,'%i',[1 1]);
  
  fposchar(fp,'=');
  ncross=fscanf(fp,'%i',[1 1]);
  c.ncross=ncross;
  
  fposchar(fp,'=');
  npt=fscanf(fp,'%i',[1 1]);
  c.npt=npt;
  
  cmpn=zeros([ncross npt 3]);
  
  for jj=1:ncross
    for kk=1:npt
      for ll=1:3
        cmpn(jj,kk,ll)=fscanf(fp,'%f',[1 1]);
      end
    end
  end
  c.pts=cmpn;
  
  hrm(ii)=c;
end
fclose(fp);


end

function fposchar(fp, ch)
jnk=0;
while(jnk~=ch)
  jnk=fscanf(fp,'%1c',[1 1]);
end
end