function RMSbulletin(regions,vartype,innovtype1,innovtype2,filename,interval)
%
%
nc=netcdf(filename);

varname1=['RMS_' innovtype1 '_' vartype];
varname2=['RMS_' innovtype2 '_' vartype];

% For several regions; divide into 2D array
nrow=ceil(sqrt(prod(size(regions))));

figure(1);
clf;

% For each region
for i=1:prod(size(regions))


   vardims1=prod(size(dim(nc{varname1})));
   vardims2=prod(size(dim(nc{varname2})));
   if (vardims1~=vardims2) 
      disp(['Dimensions do not match for variables']);
      return ;
   elseif (vardims1==2)
      var1=nc{varname1}(:,regions(i));
      var2=nc{varname2}(:,regions(i));
      I=find(var1<-900); var1(I)=nan;
      I=find(var2<-900); var2(I)=nan;
   elseif (vardims1==3)
      var1=nc{varname1}(:,regions(i),interval);
      var2=nc{varname2}(:,regions(i),interval);
      I=find(var1<-900); var1(I)=nan;
      I=find(var2<-900); var2(I)=nan;
   elseif (vardims1==1)
      var1=nc{varname1}(:);
      var2=nc{varname2}(:);
      I=find(var1<-900); var1(I)=nan;
      I=find(var2<-900); var2(I)=nan;
      plot(var1,'bs','LineWidth',3,'MarkerSize',10);
      hold on;
      plot(var2,'ro','LineWidth',3,'MarkerSize',10);
      T=title([vartype ' RMS stats ']);
      xlabel('Regions')
      ylabel(vartype)
      legend(innovtype1,innovtype2)
      return
   end

   regname=nc{'RegionName'}(regions(i),:);

   subplot(nrow,nrow,i)
   % Only works for surface vars
   plot(var1,'b','LineWidth',2);
   %axis([0 10 0 4])
   %set(gca,'FontSize',14);
   %set(gca,'FontWeight','b');
   hold on;
   plot(var2,'r','LineWidth',2);
   axis([0 10 0 4])
   if (vardims1==2)
      T=title([vartype ' RMS stats for ' regname]);
   elseif (vardims1==3)
      %'jepp'
      %T=title([vartype ' RMS stats for ' regname ' interval ' ...
      %         num2str(interval)]);
      T=title([vartype ' RMS stats for  interval ' ...
               num2str(interval) ' ' regname]);
   end
   set(T,'Interpreter','none');
   grid on;

   % Only first box
   if (i==1)
      legend(innovtype1,innovtype2)
   end

   %Only Last Row
   if (floor((i-1)/nrow) == (nrow-1))
      xlabel('Forecast Length [days]')
   end 

   % Only first column
   if (mod(i,nrow)==1)
      ylabel(vartype)
   end
end
close(nc);



