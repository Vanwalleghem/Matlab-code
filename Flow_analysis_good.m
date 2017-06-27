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

back=    [56 256 556 1006 1106 1466 1826 2086]; %Infuse
fwd=    [156 406 756 856 1256 1316 1526 1626 1986 2236]; %Withdraw

Direction_Selective_clusters=zeros(length(GoodBetas_ZS),1);
for i=1:length(Direction_Selective_clusters)
    temp=mean(ZS
end
