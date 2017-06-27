close all
Calcium=zscore(GoodCalcium,1,2);
xx = pdist(Calcium,'correlation'); %Calcium should be Neurons/time
ND=size(Calcium,1);
%To get the distance between the ith and jth observations (i < j), use the formula xx((i-1)*(m-i/2)+j-i). m is size of dataset

sortDistRow = sort(xx);
percent=2.0;
position=round(length(xx)*percent/100)
dc = sortDistRow(position);
clear sortDistRow;
rho=zeros(1,ND);
maxd=max(max(xx));

for i=1:ND-1
  for j=i+1:ND
     rho(i)=rho(i)+exp(-(xx((i-1)*(ND-i/2)+j-i)/dc)*(xx((i-1)*(ND-i/2)+j-i)/dc));
     rho(j)=rho(j)+exp(-(xx((i-1)*(ND-i/2)+j-i)/dc)*(xx((i-1)*(ND-i/2)+j-i)/dc));
  end
end

[rho_sorted,ordrho]=sort(rho,'descend');
delta(ordrho(1))=-1.;
nneigh(ordrho(1))=0;

%((ordrho(ii)-1)*(ND-ordrho(ii)/2)+ordrho(jj)-ordrho(ii))

for ii=2:ND
   delta(ordrho(ii))=maxd;
   for jj=1:ii-1
     if(xx(((ordrho(ii)-1)*(ND-ordrho(ii)/2)+ordrho(jj)-ordrho(ii)))<delta(ordrho(ii)))
        delta(ordrho(ii))=xx(((ordrho(ii)-1)*(ND-ordrho(ii)/2)+ordrho(jj)-ordrho(ii)));
        nneigh(ordrho(ii))=ordrho(jj);
     end
   end
end
delta(ordrho(1))=max(delta(:));
disp('Generated file:DECISION GRAPH')
disp('column 1:Density')
disp('column 2:Delta')

fid = fopen('DECISION_GRAPH', 'w');
for i=1:ND
   fprintf(fid, '%6.2f %6.2f\n', rho(i),delta(i));
end
clearvars i j

disp('Select a rectangle enclosing cluster centers')
scrsz = get(0,'ScreenSize');
figure('Position',[6 72 scrsz(3)/4. scrsz(4)/1.3]);
for i=1:ND
  ind(i)=i;
  gamma(i)=rho(i)*delta(i);
end
subplot(2,1,1)
tt=plot(rho(:),delta(:),'o','MarkerSize',5,'MarkerFaceColor','k','MarkerEdgeColor','k');
title ('Decision Graph','FontSize',15.0)
xlabel ('\rho')
ylabel ('\delta')
subplot(2,1,1)
rect = getrect(1);
rhomin=rect(1);
deltamin=rect(2);
NCLUST=0;
cl=-ones(1,ND);
for i=1:ND
  if ( (rho(i)>rhomin) && (delta(i)>deltamin))
     NCLUST=NCLUST+1;
     cl(i)=NCLUST;
     icl(NCLUST)=i;
  end
end
fprintf('NUMBER OF CLUSTERS: %i \n', NCLUST);
disp('Performing assignation')

%assignation
for i=1:ND
  if (cl(ordrho(i))==-1)
    cl(ordrho(i))=cl(nneigh(ordrho(i)));
  end
end


%halo
halo=cl;

if (NCLUST>1)
  bord_rho=zeros(1,NCLUST);
  for i=1:ND-1
    for j=i+1:ND
      if ((cl(i)~=cl(j))&& (xx((i-1)*(ND-i/2)+j-i)<=dc))
        rho_aver=(rho(i)+rho(j))/2.;
        if (rho_aver>bord_rho(cl(i))) 
          bord_rho(cl(i))=rho_aver;
        end
        if (rho_aver>bord_rho(cl(j))) 
          bord_rho(cl(j))=rho_aver;
        end
      end
    end
  end
  for i=1:ND
    if (rho(i)<bord_rho(cl(i)))
      halo(i)=0;
    end
  end
end
for i=1:NCLUST
  nc=0;
  nh=0;
  for j=1:ND
    if (cl(j)==i) 
      nc=nc+1;
    end
    if (halo(j)==i) 
      nh=nh+1;
    end
  end
  fprintf('CLUSTER: %i CENTER: %i ELEMENTS: %i CORE: %i HALO: %i \n', i,icl(i),nc,nh,nc-nh);
end

cmap=colormap;
for i=1:NCLUST
   ic=int8((i*64.)/(NCLUST*1.));
   subplot(2,1,1)
   hold on
   plot(rho(icl(i)),delta(icl(i)),'o','MarkerSize',8,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
end
subplot(2,1,2)
disp('Performing 2D nonclassical multidimensional scaling')
Y1 = mdscale(xx, 2, 'criterion','metricstress');
plot(Y1(:,1),Y1(:,2),'o','MarkerSize',2,'MarkerFaceColor','k','MarkerEdgeColor','k');
title ('2D Nonclassical multidimensional scaling','FontSize',15.0)
xlabel ('X')
ylabel ('Y')
for i=1:ND
 A(i,1)=0.;
 A(i,2)=0.;
end
for i=1:NCLUST
  nn=0;
  ic=int8((i*64.)/(NCLUST*1.));
  for j=1:ND
    if (halo(j)==i)
      nn=nn+1;
      A(nn,1)=Y1(j,1);
      A(nn,2)=Y1(j,2);
    end
  end
  hold on
  plot(A(1:nn,1),A(1:nn,2),'o','MarkerSize',2,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
end

%for i=1:ND
%   if (halo(i)>0)
%      ic=int8((halo(i)*64.)/(NCLUST*1.));
%      hold on
%      plot(Y1(i,1),Y1(i,2),'o','MarkerSize',2,'MarkerFaceColor',cmap(ic,:),'MarkerEdgeColor',cmap(ic,:));
%   end
%end
faa = fopen('CLUSTER_ASSIGNATION', 'w');
disp('Generated file:CLUSTER_ASSIGNATION')
disp('column 1:element id')
disp('column 2:cluster assignation without halo control')
disp('column 3:cluster assignation with halo control')
for i=1:ND
   fprintf(faa, '%i %i %i\n',i,cl(i),halo(i));
end

figure;
Mean_cl=zeros(NCLUST,size(GoodCalcium,2));
for i=1:NCLUST
    idx=find(cl==i);
    Mean_cl(i,:)=mean(Calcium(idx,:),1);
end
plot(Mean_cl');
