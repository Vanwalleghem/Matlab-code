GMModels = {};
options = statset('MaxIter',500);
for k = 81:130
    GMModels{k} = fitgmdist(temp,k,'Options',options,'CovarianceType','diagonal','Options',options, 'Regularize', 1e-5);
    BIC(k)= GMModels{k}.BIC;
end
[minBIC,numComponents] = min(BIC);
numComponents
BIC_smooth=smooth(BIC');
figure;plot(BIC_smooth);

GMModel=fitgmdist(ZS,78,'Options',options,'CovarianceType','diagonal','Options',options, 'Regularize', 1e-5);

x = linspace(0.2,size(ZS,2)/5,size(ZS,2));
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 900]);
y = linspace(1,size(ZS,1),size(ZS,1));
imagesc(x,y,ZS(randperm(size(ZS,1)),:),[0 4]);colormap hot;set(gca,'YTickLabel',[]);

GoodBetas=GoodBetas_ZS;
x = linspace(1,size(Cmap,2),size(Cmap,2));
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);
counter=1;xplot=floor(sqrt(length(GoodBetas)));yplot=ceil(length(GoodBetas)/xplot);
for i=GoodBetas   
    NumberOfCells=length(find(idxKmeans==i));
    %subplot(5,1,counter);plot(x,Cmap(i,:),x,Model_DF(i).Fitted);title(num2str(NumberOfCells))
    subplot(yplot,xplot,counter);plot(Cmap(i,:));title(num2str(NumberOfCells))
    hold on;plot((FinalFlow/1000)-0.01);
    xlim([0 size(Cmap,2)])
    counter=counter+1;
end

%GoodBetas_select=GoodBetas([1 9 4 10 6 5 7 8]);
GoodBetas_select=GoodBetas([1 9 4 10 2 5 7 8]);

[idx,nlogl,P,logpdf,M] = cluster(GMModels{numComponents},ZS);

for k = 1:length(GMModels)
    BIC(k)=GMModels{k}.BIC;
end
plot(BIC);
BIC_smooth=smooth(BIC');
plot(BIC_smooth);
[~,minBIC]=min(BIC_smooth);

numComponents=78;
GMModels{numComponents}=GMModel;
Cmap_GM=GMModels{numComponents}.mu;
idxKmeans_GM=cluster(GMModels{numComponents},ZS);
[Model_GM,GoodBetas_GM]=Test_Regress(Cmap_GM,NewFlow,idxKmeans_GM,0.3);

NewFlow=zeros(6,size(ZS,2));
GCaMP6=[-0.104392135015146,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
%GCaMP6=[0.000256990000000000;0.00850739000000000;0.0654158300000000;0.0784609000000000;0.0764130100000000;0.0665958600000000;0.0579028900000000;0.0467942900000000;0.0232079800000000;0.0144564400000000;0.00695772000000000;0.00526551000000000;0.00299500000000000;0.00198520000000000;0.00128512000000000;0.00134175000000000;0.000403170000000000;0];
back=    [56 256 556 1006 1106 1466 1826 2086]; %Infuse
back_off=[106 356 706 1056 1206 1516 1926 2176];
fwd=    [156 406 756 856 1256 1316 1526 1626 1986 2236]; %Withdraw
fwd_off=[206 506 806 956 1306 1366 1576 1776 2036 2286] ;
for i=1:length(back)
    NewFlow(1,back(i):back(i)+size(GCaMP6,1)-1)=GCaMP6';
    NewFlow(2,back_off(i):back_off(i)+size(GCaMP6,1)-1)=GCaMP6';
    NewFlow(5,back(i):back_off(i))=1;
end
for i=1:length(fwd)
    NewFlow(3,fwd(i):fwd(i)+size(GCaMP6,1)-1)=GCaMP6';
    NewFlow(4,fwd_off(i):fwd_off(i)+size(GCaMP6,1)-1)=GCaMP6';
    NewFlow(6,fwd(i):fwd_off(i))=1;
end
clearvars GCaMP6 back back_off fwd fwd_off;

Threshold=0.3;
[Model_ZS,GoodBetas_ZS]=Test_Regress(Cmap_ZS,NewFlow,idxKmeans_ZS,Threshold);

GoodBetas_select=GoodBetas_ZS([4 5 8 9 10]);
x = linspace(1,size(Cmap_ZS,2),size(Cmap_ZS,2));
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;xplot=floor(sqrt(length(GoodBetas_select)));yplot=ceil(length(GoodBetas_select)/xplot);
for i=GoodBetas_select    
    NumberOfCells=length(find(idxKmeans_ZS==i));
    %subplot(5,1,counter);plot(x,Cmap_ZS(i,:),x,Model_DF(i).Fitted);title(num2str(NumberOfCells))
    subplot(xplot,yplot,counter);plot(Cmap_ZS(i,:));title(num2str(NumberOfCells))
    %subplot(xplot,yplot,counter);imagesc(DF(find(idxKmeans_ZS==i),:),[0 0.3]);colormap hot;title(num2str(NumberOfCells))
    xlim([0 size(Cmap_ZS,2)])
    counter=counter+1;
end

ModelResultsSeg_ZS=[];
parfor i=1:length(ZS)
    %mdl=stepwiselm(Regressor',DataSet(i,:),'linear','Criterion','adjrsquared','Intercept',false,'Upper','interactions','Verbose',0);
    mdl=stepwiselm(NewFlow',ZS(i,:),'linear','Criterion','adjrsquared','Upper','linear','Verbose',0);
    ModelResultsSeg_ZS(i).coef=mdl.Coefficients;
    ModelResultsSeg_ZS(i).MSE=mdl.MSE;
    ModelResultsSeg_ZS(i).Fitted=mdl.Fitted;
    ModelResultsSeg_ZS(i).rsquared=mdl.Rsquared.Adjusted;
end
idx_rsq_ZS=find([ModelResultsSeg_ZS.rsquared]>0.1);
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);
imagesc(ZS(idx_rsq_ZS,:),[0 5]);colormap hot
temp=ZS(idx_rsq_ZS,:);

options = statset('UseParallel',1); [idxKmeans_ZS_filter Cmap_ZS_filter]=kmeans(temp,50,'Options',options,'Distance','correlation','Replicates',5,'MaxIter',1000,'Display','final');
[Model_ZS,GoodBetas_ZS]=Test_Regress(Cmap_ZS_filter,NewFlow,idxKmeans_ZS_filter,Threshold);

%Nice presentation of flow stimuli
NewFlow=zeros(2,500);
back=    [56 256 556 1006 1106 1466 1826 2086]; %Infuse
back_off=[106 356 706 1056 1206 1516 1926 2176];
fwd=    [156 406 756 856 1256 1316 1526 1626 1986 2236]; %Withdraw
fwd_off=[206 506 806 956 1306 1366 1576 1776 2036 2286] ;
%back=back/5;back_off=back_off/5;
%fwd=fwd/5;fwd_off=fwd_off/5;
slow_fast=[1 2 3 2 1 2 2 1];
slow_fast_fwd=[2 1 1 2 2 2 2 3 1 1];
for i=1:length(back)
    NewFlow(1,back(i):back_off(i))=slow_fast(i);
end
for i=1:length(fwd)    
    NewFlow(2,fwd(i):fwd_off(i))=slow_fast_fwd(i);
end
clearvars GCaMP6 back back_off fwd fwd_off;

find(NewFlow(1,:)==3);
temp=zeros(size(ans));
%temp(1:10)=1;temp(11:20)=2;temp(21:31)=1;
temp(1:50)=1;temp(51:100)=2;temp(101:151)=1;
NewFlow(1,ans)=temp;


find(NewFlow(2,:)==3);
temp=zeros(size(ans));
%temp(1:10)=1;temp(11:20)=2;temp(21:31)=1;
temp(1:50)=1;temp(51:100)=2;temp(101:151)=1;
NewFlow(2,ans)=temp;

FinalFlow=NewFlow(1,:)-NewFlow(2,:);
FinalFlow=FinalFlow*5;

parfor i=1:size(ZS,1)
    mdl=stepwiselm(NewFlow',ZS(i,:),'Upper','linear','Intercept',false,'Criterion','bic','verbose',0);
    model_ZS(i).coef=mdl.Coefficients;
    model_ZS(i).MSE=mdl.MSE;
    model_ZS(i).Fitted=mdl.Fitted;
    model_ZS(i).rsquared=mdl.Rsquared.Adjusted;
end

Threshold=0.1;
idx_rsq_ZS=find([model_ZS.rsquared]>Threshold);ZS_rsq=ZS(idx_rsq_ZS,:);
figure;imagesc(ZS(idx_rsq_ZS,:),[0 4]);colormap hot
options = statset('UseParallel',1); [idxKmeans_ZS Cmap_ZS]=kmeans(ZS(idx_rsq_ZS,:),10,'Options',options,'Distance','correlation','Replicates',10,'MaxIter',2000,'Display','final');
[Model_ZS,GoodBetas_ZS]=Test_Regress(Cmap_ZS,NewFlow,idxKmeans_ZS,0.2);

Flow_profile=[0,0,0,0,0,0,0,0,0,0,5,5,5,5,5,5,5,5,5,5,5,0,0,0,0,0,0,0,0,0,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,0,0,0,0,0,0,0,0,0,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,0,0,0,0,0,0,0,0,0,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,0,0,0,0,0,0,0,0,0,5,5,5,5,5,5,5,5,5,5,10,10,10,10,10,10,10,10,10,10,5,5,5,5,5,5,5,5,5,5,5,0,0,0,0,0,0,0,0,0,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,0,0,0,0,0,0,0,0,0,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,0,0,0,0,0,0,0,0,0,10,10,10,10,10,10,10,10,10,10,10,0,0,0,0,0,0,0,0,0,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,0,0,0,0,0,0,0,0,0,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,0,0,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,10,10,10,10,10,10,10,10,10,10,10,0,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,0,0,0,0,0,0,0,0,0,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-10,-10,-10,-10,-10,-10,-10,-10,-10,-10,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,0,0,0,0,0,0,0,0,0,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,0,0,0,0,0,0,0,0,0,0,0,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,0,0,0,0,0,0,0,0,0,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,0,0,0,0,0,0,0,0,0,0,0,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
Flow_profile2=interp(Flow_profile,5);Flow_profile2=int8(Flow_profile2);Flow_profile2(Flow_profile2==1)=0;Flow_profile2(Flow_profile2==-1)=0;Flow_profile2(Flow_profile2==6)=5;Flow_profile2(Flow_profile2==-6)=-5;Flow_profile2(Flow_profile2==11)=10;Flow_profile2(Flow_profile2==-11)=-10;Flow_profile2(Flow_profile2==-4)=-5;Flow_profile2(Flow_profile2==4)=5;Flow_profile2(Flow_profile2==2)=0;


Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;xplot=floor(sqrt(length(GoodBetas_ZS)));yplot=ceil(length(GoodBetas_ZS)/xplot);
for i=GoodBetas_ZS
    idx_temp=find(idxKmeans_ZS==i);
    subplot(xplot,yplot,counter);plot(x,mean(ZS_rsq(idx_temp,:),1));hold on;plot(x,(double(Flow_profile2)/10)-1)
    counter=counter+1;
end

for k = 61:80
GMModels{k} = fitgmdist(ZS_rsq,k,'Options',options,'CovarianceType','diagonal','Options',options, 'Regularize', 1e-5);
BIC(k)= GMModels{k}.BIC;
end
BIC_smooth=smooth(BIC');
figure;plot(BIC_smooth);

%To compute DSI, only look at "simple" stimuli
back= [56  256 1006 1106 1826 2086]; %Infuse
fwd=  [156 406 756  856  1986 2236]; %Withdraw
Direction_Selective_clusters=zeros(length(GoodBetas_ZS),2);
for i=1:length(Direction_Selective_clusters)
    idx_temp=find(idxKmeans_ZS==GoodBetas_ZS(i));
    temp=mean(ZS_rsq(idx_temp,:),1);
    back_str=zeros(length(back),1);
    for j=1:length(back)
        back_str(j)=max(temp(back(j):back(j)+20));
    end
    fwd_str=zeros(length(fwd),1);
    for j=1:length(fwd)
        fwd_str(j)=max(temp(fwd(j):fwd(j)+20));
    end
    Direction_Selective_clusters(i,1)=mean(back_str)/mean(fwd_str);
    Direction_Selective_clusters(i,2)=mean(fwd_str)/mean(back_str);
    %Direction_Selective_clusters(i,2)=(((std(back_str)/mean(back_str))^2)+(std(fwd_str)/mean(fwd_str))^2)^1/2;
end

Numbers=[1 [MatFiles.GoodNumber]];
counter=1;
idx_Plane=nan(length(ZS),1);
idx_Fish=nan(length(ZS),1);
name=strcat(MatFiles(1).name);
% %[Plane,~]=regexp(name,'\d+_(\d+)_','tokens','match');Plane=str2num(Plane{1}{1});
% [Plane,~]=regexp(name,'\d\D(\d+)um','tokens','match');Plane=str2num(Plane{1}{1});
% %[Fish,~]=regexp(name,'(\d+)_\d+_','tokens','match');Fish=str2num(Fish{1}{1});
% [Fish,~]=regexp(name,'(\d)\D\d+um','tokens','match');Fish=str2num(Fish{1}{1});
% idx_Plane(1:Numbers(2))=Plane;
% idx_Fish(1:Numbers(2))=Fish;
for i=1:length(MatFiles)
    %[Fish,~]=regexp(files{i},'(\d+)_','tokens','match');Fish=str2num(Fish{1}{1});
    name=strcat(MatFiles(i).name);    
    [Plane,~]=regexp(name,'_(\d+)um','tokens','match');Plane=str2num(Plane{1}{1});
    [Fish,~]=regexp(name,'Fish(\d)_','tokens','match');Fish=str2num(Fish{1}{1});
    if Plane>1000
        Plane=Plane-12376;
    end
    idx_Plane(Numbers(i):Numbers(i+1))=Plane;
    idx_Fish(Numbers(i):Numbers(i+1))=Fish;
end
clearvars i Fish Plane name counter


counter=1;counter2=1;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
GoodClusterData=[];
%name=['ELO';'CLO'];
for i=GoodBetas_ZS  
    idx=find(idxKmeans_ZS==i);
    GoodClusterData(counter).mean=mean(ZS_rsq(idx,:),1);
    GoodClusterData(counter).ZS=ZS_rsq(idx,:);
    GoodClusterData(counter).planes=idx_Plane(idx_rsq_ZS(idx));
    GoodClusterData(counter).fish=idx_Fish(idx_rsq_ZS(idx));
    subplot(length(GoodBetas_ZS),3,counter2);plot(GoodClusterData(counter).mean)%,'color',colors(counter,:));title(num2str(length(idx))),xlim([0 size(ZS,2)]),ylim([-1 4])    
    subplot(length(GoodBetas_ZS),3,counter2+1);imagesc(GoodClusterData(counter).ZS, [0 4]); colormap hot;title(num2str(length(idx)))    
    subplot(length(GoodBetas_ZS),3,counter2+2);histogram(GoodClusterData(counter).planes);
    %subplot(length(GoodBetas_ZS),4,counter2+3);histogram(GoodClusterData(counter).fish);
    %subplot(length(GoodBetas_select),4,counter2+3);bar([sum(GoodClusterData(counter).State==3) sum(GoodClusterData(counter).State==4)]);set(gca,'xticklabel',name);
    counter2=counter2+3;
    counter=counter+1;
end

coefficients={};
for idx=1:length(model_ZS)
    coef=[model_ZS(idx).coef];
    temp=coef.Properties.RowNames;temp=regexp(temp,'x(\d)','tokens');
    if ~isempty(temp)
        temp=[temp{:}];temp=[temp{:}];temp=[temp{:}];%temp=str2num(temp);
        for coef_idx=1:height(coef)
            if coef.pValue(coef_idx)<0.05
                coefficients{idx,str2num(temp(coef_idx))}=coef.Estimate(coef_idx);
            end
        end
    end
end
idxempty=cellfun('isempty',coefficients);
coefficients(idxempty)={0};
clearvars idxempty idx coef_idx coef
coefficients=cell2mat(coefficients);
coefficients_rsq=coefficients(idx_rsq_ZS,:);

eva = evalclusters(coefficients_rsq,'kmeans','CalinskiHarabasz','KList',[1:30]);
options = statset('UseParallel',1); [idxKmeans_coef Cmap_coef]=kmeans(coefficients_rsq,10,'Options',options,'Distance','cityblock','Replicates',10,'MaxIter',2000,'Display','final');

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;xplot=floor(sqrt(size(Cmap_coef,1)));yplot=ceil(size(Cmap_coef,1)/xplot);
for i=1:size(Cmap_coef,1)
    idx_temp=find(idxKmeans_coef==i);
    subplot(xplot,yplot,counter);plot(x,mean(ZS_rsq(idx_temp,:),1));hold on;plot(x,(double(Flow_profile2)/10)-1)
    counter=counter+1;
end


counter=1;counter2=1;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
%name=['ELO';'CLO'];
%for i=1:size(Cmap_coef,1)

for i=[7 4 6 2 5 1 9 3 8 10]
    idx_temp=find(idxKmeans_coef==i);  
    subplot(length(GoodBetas_ZS),3,counter2);plot(x,mean(ZS_rsq(idx_temp,:),1),'LineWidth',2);xlim([0 size(ZS,2)]);ylim([-1 3]);%,'color',colors(counter,:));title(num2str(length(idx))),xlim([0 size(ZS,2)]),ylim([-1 4])    
    subplot(length(GoodBetas_ZS),3,counter2+1);imagesc(ZS_rsq(idx_temp,:), [0 3]); colormap hot;title(num2str(length(idx_temp)))    
    subplot(length(GoodBetas_ZS),3,counter2+2);histogram(idx_Plane(idx_rsq_ZS(idx_temp)));
    %subplot(length(GoodBetas_ZS),4,counter2+3);histogram(GoodClusterData(counter).fish);
    %subplot(length(GoodBetas_select),4,counter2+3);bar([sum(GoodClusterData(counter).State==3) sum(GoodClusterData(counter).State==4)]);set(gca,'xticklabel',name);
    counter2=counter2+3;
    counter=counter+1;
end

%second stimulus complex (~150) = second basic stimulus (~155)+80
%eight complex stimulus (~1005) = first basic stimulus (~255)

Corr_BasicCLust=zeros(size(Basic_Clusters,1),size(ZS_rsq,1));

parfor i=1:size(ZS_rsq,1)
    corr_temp=[];
    ZS_rsq_temp=[ZS_rsq(i,start_complex(1):start_complex(1)+80) ZS_rsq(i,start_complex(2):start_complex(2)+80)];
    for j=1:size(Basic_Clusters,1)
        Basic_temp=[Basic_Clusters(j,start_basic(1):start_basic(1)+80) Basic_Clusters(j,start_basic(2):start_basic(2)+80)];
        temp=corrcoef(ZS_rsq_temp, Basic_temp);
        corr_temp(j)=temp(1,2);        
    end
    Corr_BasicCLust(:,i)=corr_temp;
end

[MaxCorr_BasiClust MaxCorr_BasiClust_ind]=max(abs(Corr_BasicCLust),[],1);
idx_BasicClust_corr=find(MaxCorr_BasiClust>0.7);
BasicClustData=[];
for basic=1:max(MaxCorr_BasiClust_ind)
    idx_BasicClust_temp=find(MaxCorr_BasiClust_ind(idx_BasicClust_corr)==basic);
    BasicClustData(basic).ZS=ZS_rsq(idx_BasicClust_corr(idx_BasicClust_temp),:);
    BasicClustData(basic).mean=mean(BasicClustData(basic).ZS,1);    
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;xplot=floor(sqrt(max(MaxCorr_BasiClust_ind)));yplot=ceil(max(MaxCorr_BasiClust_ind)/xplot);
for i=1:max(MaxCorr_BasiClust_ind)    
    subplot(xplot,yplot,counter);plot(x,BasicClustData(i).mean);hold on;plot(x,(double(Flow_profile2)/10)-1)
    counter=counter+1;
end

%Basic cluster coloring

All_ROIs=[];
ROIs_idx=[];
for i = 1:length(MatFiles)
    name=strcat(MatFiles(i).name);
    Rs=load(name, 'ROIs');
    Rs=Rs.ROIs;
    F=load(name, 'idx_components');
    F=F.idx_components+1;
    Rs=Rs(:,F);
    All_ROIs{i}=Rs;
    if i==1
        ROIs_idx(i)=length(F);
    else
        ROIs_idx(i)=ROIs_idx(i-1)+length(F);
    end
end
clearvars GC C S F N name i;

Numbers=[0 [ROIs_idx]];
temp=[];
for basic=1:max(MaxCorr_BasiClust_ind)
    idx_BasicClust_temp=find(MaxCorr_BasiClust_ind(idx_BasicClust_corr)==basic);
    temp{basic}=idx_rsq_ZS(idx_BasicClust_corr(idx_BasicClust_temp));
end
clearvars i basic j k l m n o p

Start=min(cellfun(@min, temp));Start=find(Numbers<Start,1,'last');
filename=MatFiles(Start).name;

Test=zeros(8,size(ZS,2));
for basic=1:max(MaxCorr_BasiClust_ind)
    Test(basic,:)=mean(ZS(temp{basic},:),1);    
end

%colors = [1,0,0;0,1,0;0,0,1;1,0.600000000000000,0;1,0,0.7];
%colors = distinguishable_colors(9,[1 1 1; 0 0 0]);
colors = [0         0    1.0000
         0    0.5000    1.0000
    1.0000         0         0
    1.0000    0.1034    0.7241
    1.0000    0.5000    0.3000
         0    0.7000    0.2000
    0.5000    0.5000         0
         0    0.5000    0.5000];
colors = colors*256;
for idx=Start:length(MatFiles)
    filename=MatFiles(idx).name;
    ROIsNb=[];ClusterNb=[];
    %for k = 1 : length(temp)
    for k = 1 : length(temp)
        tempROIsNb=find([temp{k}]<=Numbers(idx+1));
        if tempROIsNb            
            ROIsNb=[ROIsNb temp{k}(tempROIsNb)];
            temp{k}(tempROIsNb)=[];
            ClusterNb=[ClusterNb ; repmat(k,length(tempROIsNb),1)];
        end
    end
    if ROIsNb
        imagename=regexp(filename,'_output_analysis','split');
        %imagename=regexp(imagename,'_output_analysis_matlab2.mat','split');
        imagename=strcat('AVG_',imagename{1},'.tif');
        image=double(imread(imagename));image=image/max(max(image));image=image*128;
        image=uint8(image);
        image2=zeros(size(image(:,:,1)));
        image3=repmat(image,1,1,3);
        ROIs=All_ROIs{idx};       
        ROIsNb=ROIsNb-Numbers(idx);
        ROIs=ROIs(:,ROIsNb);
        for k = 1 : size(ROIs,2)
            image2=zeros(size(image(:,:,1)));
            ROI=full(reshape(ROIs(:,k),size(image,1),size(image,2)));            
            image2=image2+ROI;image2=(image2/max(max(image2)));image2=uint8(image2);
            for j=1:3
                image3(:,:,j)=image3(:,:,j)+image2*colors(ClusterNb(k),j);
            end
        end
        %image3(:,:,3)=image;
        name=strcat('BasicClust_',imagename(4:end));
    imwrite(image3,name,'tif');
    end
    %image3=uint8(image3);

end
clearvars idx i temp tempFileNb fileNb AVG_files filename image counter Numbers image2 image3 k ROI ROIs ROIsNb Start tempROIsNb name imagename tempidx Raster

ZS_diff_rsq = zeros(size(ZS_rsq));
parfor i=1:size(ZS_rsq,1)
    temp = TVRegDiff(ZS_rsq(i,:), 50, 1e-1,[],'large',1e-8,[],0,0);
    ZS_diff_rsq(i,:)=temp;
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);x = linspace(0.2,size(Cmap,2)/5,size(Cmap,2));
counter=1;counter2=1;xplot=floor(sqrt(max(MaxCorr_BasiClust_ind)));yplot=ceil(max(MaxCorr_BasiClust_ind)/xplot);
for i=1:max(MaxCorr_BasiClust_ind)
    if counter==3
        counter=counter+1;
    end
    subplot(3,3,counter);plot(x,BasicClustData(i).mean,'color',colors(counter2,:)/256);%axis([0 131 -1 4]);rectangle('FaceColor','r','Position',[11 -1 10 0.25]);rectangle('FaceColor','r','Position',[51 -1 10 0.25]);rectangle('FaceColor','r','Position',[91 -1 10 0.25]);rectangle('FaceColor','b','Position',[31 -1 10 0.25]);rectangle('FaceColor','b','Position',[71 -1 10 0.25]);rectangle('FaceColor','b','Position',[111 -1 10 0.25]);
    counter=counter+1;
    counter2=counter2+1;
end
