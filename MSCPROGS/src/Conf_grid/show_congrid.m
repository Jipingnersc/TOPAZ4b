
idm=1400;  jdm=1320;

mlon=loada('regional.grid.a',1,idm,jdm);
mlat=loada('regional.grid.a',2,idm,jdm);
mdep=loada('regional.depth.a',1,idm,jdm);
mask_0=find(mdep<=1|mdep>20000);
var0=mdep;  var0(mask_0)=nan;

m_proj('Satellite','lon',-15,'lat',50,'radius',50)

figure('color',[1 1 1])
  m_coast('patch',[.6 .6 .6]); m_grid; hold on;
  P0=m_pcolor(mlon,mlat,var0);  set(P0,'linestyle','none');
 

%m_proj('Lambert','lon',[-40 98],'lat',[58 86])
m_proj('Lambert','lon',[-20 70],'lat',[60 85.5],'rec','on')
%m_proj('Mercator','lon',[-25 98],'lat',[58 83])

figure('color',[1 1 1])
  m_coast('patch',[.6 .6 .6]); m_grid; hold on;
  P0=m_pcolor(mlon,mlat,var0);  set(P0,'linestyle','none');
 


