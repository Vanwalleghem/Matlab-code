MatFiles=dir('*matlab*.mat');
name=strcat(MatFiles(1).name);
Calcium=load(name, 'DenoisedTraces');
Calcium=Calcium.DenoisedTraces(:,1:2500);
MatFiles(1).number=size(Calcium,1);
Spikes=load(name, 'Spikes');
Spikes=Spikes.Spikes(:,1:2500);
Noise=load(name, 'Noise');
Noise=Noise.Noise(:,1:2500);
Fitness=load(name, 'idx_components');
Fitness=Fitness.idx_components+1;
GoodCalcium=Calcium;
GoodSpikes=Spikes;
MatFiles(1).GoodNumber=length(Fitness);
for i = 2:length(MatFiles)
name=strcat(MatFiles(i).name);
C=load(name, 'DenoisedTraces');
C=C.DenoisedTraces(:,1:2500);
%     if i==3
%         C=[C(:,1) C(:,1) C(:,1:58)];
%     end
S=load(name, 'Spikes');
S=S.Spikes(:,1:2500);
N=load(name, 'Noise');
N=N.Noise(:,1:2500);
F=load(name, 'idx_components');
F=F.idx_components+1;
GC=C(F,:);
GS=S(F,:);
Noise=vertcat(Noise,N);
Calcium=vertcat(Calcium,C);
Spikes=vertcat(Spikes,S);
Fitness=horzcat(Fitness,F);
GoodCalcium=vertcat(GoodCalcium,GC);
GoodSpikes=vertcat(GoodSpikes,GS);
MatFiles(i).number=size(Calcium,1);
MatFiles(i).GoodNumber=MatFiles(i-1).GoodNumber+length(F);
MatFiles(i).GC=GC;
end
clearvars GC C S F N name i GS;

ZS=zscore(GoodCalcium,1,2);
x = linspace(0.2,size(ZS,2)/5,size(ZS,2));y = linspace(1,size(ZS,1),size(ZS,1));
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 900]);
imagesc(x,y,ZS,[0 5]);colormap hot;set(gca,'YTickLabel',[]);

Numbers=[1 [MatFiles.GoodNumber]];
counter=1;
idx_Plane=nan(length(GoodCalcium),1);
idx_Fish=nan(length(GoodCalcium),1);
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
    if isempty(regexp(name,'\d\D(\d+)um','tokens'))
        [Plane,~]=regexp(name,'_(\d)um_','tokens','match');Plane=str2num(Plane{1}{1});
        Fish=1;
    else
        [Plane,~]=regexp(name,'\d\D(\d+)um','tokens','match');Plane=str2num(Plane{1}{1});
        [Fish,~]=regexp(name,'(\d)\D\d+um','tokens','match');Fish=str2num(Fish{1}{1});
    end
    %[Plane,~]=regexp(name,'\d+_(\d+)_','tokens','match');Plane=str2num(Plane{1}{1});
    %[Fish,~]=regexp(name,'(\d+)_\d+_','tokens','match');Fish=str2num(Fish{1}{1});
   
    idx_Plane(Numbers(i):Numbers(i+1))=Plane;
    idx_Fish(Numbers(i):Numbers(i+1))=Fish;
end
clearvars i Fish Plane name counter

%for Fish 16 : load('FlowAnalysis_Kmeans_100clust.mat', 'FluoTraces')
FluoTracesb=FluoTraces;
FluoTraces=[FluoTracesb; FluoTraces16];

DF=DeltaF2(double(FluoTraces),21,11);
max_DF=max(DF,[],2);
min_DF=min(DF,[],2);
idx_Threshold_5to200=find(max_DF>0.05 & max_DF<2 & min_DF>-0.1);
DF_5to200=DF(idx_Threshold_5to200,:);

flow=zeros(6,655);
GCaMP6=[-0.104392135015146,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
%GCaMP6=[0.000256990000000000;0.00850739000000000;0.0654158300000000;0.0784609000000000;0.0764130100000000;0.0665958600000000;0.0579028900000000;0.0467942900000000;0.0232079800000000;0.0144564400000000;0.00695772000000000;0.00526551000000000;0.00299500000000000;0.00198520000000000;0.00128512000000000;0.00134175000000000;0.000403170000000000;0];
back=[57 257 457];
back_off=[106 306 506];
fwd=[157 357 557];
fwd_off=[207 407 607];
flow(1,back(1):back(1)+size(GCaMP6,1)-1)=GCaMP6';
flow(1,back(2):back(2)+size(GCaMP6,1)-1)=GCaMP6';
flow(1,back(3):back(3)+size(GCaMP6,1)-1)=GCaMP6';
flow(2,back_off(1):back_off(1)+size(GCaMP6,1)-1)=GCaMP6';
flow(2,back_off(2):back_off(2)+size(GCaMP6,1)-1)=GCaMP6';
flow(2,back_off(3):back_off(3)+size(GCaMP6,1)-1)=GCaMP6';
flow(3,fwd(1):fwd(1)+size(GCaMP6,1)-1)=GCaMP6';
flow(3,fwd(2):fwd(2)+size(GCaMP6,1)-1)=GCaMP6';
flow(3,fwd(3):fwd(3)+size(GCaMP6,1)-1)=GCaMP6';
flow(4,fwd_off(1):fwd_off(1)+size(GCaMP6,1)-1)=GCaMP6';
flow(4,fwd_off(2):fwd_off(2)+size(GCaMP6,1)-1)=GCaMP6';
flow(4,fwd_off(3):fwd_off(3)+size(GCaMP6,1)-1)=GCaMP6';
flow(5,back(1):back(1)+43)=1;
flow(5,back(2):back(2)+43)=1;
flow(5,back(3):back(3)+43)=1;
flow(6,fwd(1):fwd(1)+43)=1;
flow(6,fwd(2):fwd(2)+43)=1;
flow(6,fwd(3):fwd(3)+43)=1;
clearvars GCaMP6 back back_off fwd fwd_off;

x = linspace(0.2,size(DF_5to200,2)/5,size(DF_5to200,2));
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 900]);
y = linspace(1,size(DF_5to200,1),size(DF_5to200,1));
figure;imagesc(x,y,DF_5to200(randperm(size(DF_5to200,1)),:)*100,[0 30]);colormap hot;set(gca,'YTickLabel',[]);

ZS=zscore(GoodCalcium,1,2);
% max_ZS=max(ZS,[],2);
% idx_ZS=find(max_ZS>2);
% ZS_select=ZS(idx_ZS,:);
random_5k=datasample(ZS,5000);
eva = evalclusters(random_5k,'kmeans','silhouette','Distance','correlation','KList',[1:100]);
options = statset('UseParallel',1); [idxKmeans Cmap]=kmeans(ZS,eva.OptimalK,'Options',options,'Distance','correlation','Replicates',5,'MaxIter',1000,'Display','final');

Threshold=0.3;
[Model_DF,GoodBetas]=Test_Regress(Cmap,flow,idxKmeans,Threshold);
[Sorted_rsq,order_rsq]=sort([Model_DF.rsquared]);
idx=find(Sorted_rsq>Threshold);
counter=1;xplot=floor(sqrt(length(idx)));yplot=ceil(length(idx)/xplot);
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 900]);
GoodCluster=[];
for i=idx
    %GoodCluster(counter,:)=mean(ZS(find(idxKmeans==i),:),1);
    GoodCluster(counter,:)=Cmap(order_rsq(i),:);
    subplot(xplot,yplot,counter);plot(Cmap(order_rsq(i),:));xlim([0 size(Cmap,2)])
    counter=counter+1;
end
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 900]);
imagesc(GoodCluster,[0 4]);colormap hot


x = linspace(0.2,size(DF_5to200,2)/5,size(DF_5to200,2));
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 900]);
y = linspace(1,size(DF_5to200,1),size(DF_5to200,1));
imagesc(Cmap*100,[0 30]);colormap hot;set(gca,'YTickLabel',[]);

[coeff,score,~,~,explained,~] = pca(ZS);

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1200, 900]);
imagesc(coeff(:,1:5)*100,[0 20]);colormap hot;set(gca,'YTickLabel',[]);

%GoodBetas_select=GoodBetas([2 3 5 6 8]);
GoodBetas_select=GoodBetas_GC(6);
x = linspace(1,size(Cmap,2),size(Cmap,2));
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);
counter=1;xplot=floor(sqrt(length(GoodBetas_select)));yplot=ceil(length(GoodBetas_select)/xplot);
for i=GoodBetas_select    
    NumberOfCells=length(find(idxKmeans==i));
    %subplot(5,1,counter);plot(x,Cmap(i,:),x,Model_DF(i).Fitted);title(num2str(NumberOfCells))
    subplot(xplot,yplot,counter);plot(Cmap(i,:));title(num2str(NumberOfCells))
    xlim([0 size(Cmap,2)])
    counter=counter+1;
end

idx_fwd=find(idxKmeans==16);
idx_back=find(idxKmeans==18);

idx_final_fwd=idx_Threshold_5to200(idx_fwd);
idx_final_back=idx_Threshold_5to200(idx_back);

counter=1;
idx_fish_CTRL=zeros(length(DF),1);
idx_Plane=nan(length(DF),1);
[Fish,~]=regexp(files{1},'(\d+)_','tokens','match');Fish=str2num(Fish{1}{1});
idx_fish_CTRL(1:index(1))=Fish;
[Plane,~]=regexp(files{1},'\d+_(\d+)','tokens','match');Plane=str2num(Plane{1}{1});
idx_Plane(1:index(1))=Plane;
for i=2:length(index)
	[Fish,~]=regexp(files{i},'(\d+)_','tokens','match');Fish=str2num(Fish{1}{1});
    [Plane,~]=regexp(files{i},'\d+_(\d+)','tokens','match');Plane=str2num(Plane{1}{1});
    idx_fish_CTRL(index(i-1):index(i))=Fish;
    idx_Plane(index(i-1):index(i))=Plane;
end
idx_fish_CTRL(idx_fish_CTRL==0)=16;

counter=1;
for i=GoodBetas_select
    idx_temp=find(idxKmeans==i);
    distribution_fish{counter}=idx_fish_CTRL(idx_Threshold_5to200(idx_temp));
    distribution_plane{counter}=idx_Plane(idx_Threshold_5to200(idx_temp));
    counter=counter+1;
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 300, 900]);
xplot=floor(sqrt(length(GoodBetas_select)));yplot=ceil(length(GoodBetas_select)/xplot);
for i=1:length(GoodBetas_select)
    subplot(5,1,i);histogram(distribution_fish{i})
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 300, 900]);
edges = [0:10:250];
xplot=floor(sqrt(length(GoodBetas_select)));yplot=ceil(length(GoodBetas_select)/xplot);
for i=1:length(GoodBetas_select)
    subplot(5,1,i);histogram(distribution_plane{i},edges)
end

parfor i=1:size(ZS,1)
    mdl=stepwiselm(flow',ZS(i,:),'Upper','linear','Intercept',false,'Criterion','bic','verbose',0);
    model_DF_Thr5(i).coef=mdl.Coefficients;
    model_DF_Thr5(i).MSE=mdl.MSE;
    model_DF_Thr5(i).Fitted=mdl.Fitted;
    model_DF_Thr5(i).rsquared=mdl.Rsquared.Adjusted;
end

rsq=[model_DF_Thr5.rsquared];
figure;histogram(rsq);
idx_rsq=find(rsq>0.3);
ALL_DF_02=ZS(idx_rsq,:);

Merged_DF=DF_5to200;%   ALL_DF_02=Dataset (neurons x time)
MergedIDX=cell(1,length(Merged_DF)); % Where the IDs of merged ROIs will be stored

counter=1;threshold=0.75;
progressbar % can be skipped if you haven't progressbar installed    
while counter<=size(Merged_DF,1) % can also be for counter=1:size(Dataset,1)    
    if isnan(Merged_DF(counter,1)) % if the timeserie has been merged, it's set to nan, so this prevents spending time on it
        max_corr=0;
        MergedIDX{counter}=[];
    else
        max_corr=1;    % hack to get the equivalent of a do... while loop (at least one iteration)
        if isempty(MergedIDX{counter})
            MergedIDX{counter}=counter;
        end
    end    
    while max_corr>threshold % standard threshold of Bianco et al
        corr_temp=zeros(size(Merged_DF,1),1);        
        parfor j=1:size(Merged_DF,1)    % parallelized, can replace by for loop if not needed
        	if j>counter && ~isnan(Merged_DF(j,1))
                temp=corrcoef(Merged_DF(counter,:), Merged_DF(j,:));
                corr_temp(j)=temp(1,2);                                
            end
        end        
        [max_corr idx_max]=nanmax(corr_temp);
        if max_corr>threshold
            if max_corr>(threshold+0.1) % shortcut to merge in one go all time series with >0.85 correlation, can be skipped or changed
                idx_max=find(corr_temp>(threshold+0.1));
                MergedIDX{counter}=[MergedIDX{counter} idx_max'];
            else
                MergedIDX{counter}=[MergedIDX{counter} idx_max]; % merge things one by one, so it's long
            end
            Merged_DF(counter,:)=nanmean(DF_5to200(MergedIDX{counter},:),1); %average from dataset, easier than adding one at a time with a factor                
            Merged_DF(idx_max,:)=nan;
            corr_temp(idx_max)=nan;            
        end
    end
    counter=counter+1;
    progressbar(counter/size(Merged_DF,1)); % can be skipped if you haven't progressbar installed    
end
clearvars i j idx_max counter max_corr temp

numMerged=zeros(numel(MergedIDX),1);
for i=1:numel(MergedIDX)
    numMerged(i)=numel(MergedIDX{i});
end
temp=find(numMerged>25);
Representative_clust=[];
nb_threshold=10;nb_fish_threshold=5;
for i=1:numel(temp)
    temp_idx=MergedIDX{temp(i)};
    temp_fish=idx_Fish(temp_idx);
    if numel(unique(temp_fish))>=nb_fish_threshold   
        nbcells=[];
        for j=1:max(idx_Fish)
            nbcells{j}=sum(temp_fish==j);
        end
        nbcells=cell2mat(nbcells);
        nbcells=sum(nbcells>=nb_threshold);
        if nbcells>=nb_fish_threshold
            Representative_clust=[Representative_clust temp(i)];
        end
    end
end
clearvars temp i temp_idx temp_fish nbcells

MeanClust=zeros(length(Representative_clust),size(GoodCalcium,2));
counter=1;
for idx=Representative_clust
MeanClust(counter,:)=mean(GoodCalcium(MergedIDX{idx},:),1);
counter=counter+1;
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);
edges = [0:10:250];
xplot=floor(sqrt(size(MeanClust,1)));yplot=ceil(size(MeanClust,1)/xplot);
for i=1:size(MeanClust,1)
    subplot(xplot,yplot,i);plot(MeanClust(i,:));title(num2str(length(MergedIDX{Representative_clust(i)})))
end

counter=1;
for i=1:length(GoodBetas_Bianco)
    subplot(2*length(GoodBetas_Bianco),2,counter);plot(MeanClust(GoodBetas_Bianco(i),:))
    subplot(2*length(GoodBetas_Bianco),2,counter+1);imagesc(zscore(GoodCalcium(MergedIDX{Representative_clust(GoodBetas_Bianco(i))},:),1,2));colormap hot    
    counter=counter+2;
end

coefficients={};
for idx=1:length(model_DF_Thr5)
    coef=[model_DF_Thr5(idx).coef];
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
[max_coef,idx_max]=max(coefficients,[],2);

eva = evalclusters(GoodCalcium,'kmeans','silhouette','Distance','correlation','KList',[1:100]);
if eva.OptimalK<100
    options = statset('UseParallel',1); [idxKmeans Cmap]=kmeans(GoodCalcium,eva.OptimalK,'Options',options,'Distance','correlation','Replicates',5,'MaxIter',1000,'Display','final');
else
    eva = evalclusters(GoodCalcium,'kmeans','silhouette','Distance','correlation','KList',[101:150]);
    if eva.OptimalK<150
        options = statset('UseParallel',1); [idxKmeans Cmap]=kmeans(GoodCalcium,eva.OptimalK,'Options',options,'Distance','correlation','Replicates',5,'MaxIter',1000,'Display','final');
    else
        eva = evalclusters(GoodCalcium,'kmeans','silhouette','Distance','correlation','KList',[151:200]);
        options = statset('UseParallel',1); [idxKmeans Cmap]=kmeans(GoodCalcium,eva.OptimalK,'Options',options,'Distance','correlation','Replicates',5,'MaxIter',1000,'Display','final');
    end
end
save('CNMF_Flow','-v7.3');
system('shutdown -s')

%GoodBetas_select=GoodBetas_GC([1 2 4 5 6 7 8 9 10]);
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);
counter=1;xplot=floor(sqrt(length(GoodBetas_select)));yplot=ceil(length(GoodBetas_select)/xplot);
for i=GoodBetas_select
idx=find(idxKmeans==i);
NumberOfCells=length(idx);
subplot(length(GoodBetas_select),1,counter);plot(Cmap(i,:));%,'color',colors(counter,:)/256);title(num2str(NumberOfCells))
counter=counter+1;
end

GoodBetas_select=GoodBetas([1 2 4 5 6 7 8 9 10]);
GoodBetas_select=GoodBetas_select([1 2 4 5 6 7 8 9 10]);
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);
counter=1;xplot=floor(sqrt(length(GoodBetas_select)));yplot=ceil(length(GoodBetas_select)/xplot);
GoodBetas_selectbackup=GoodBetas_select;
GoodBetas_select=GoodBetas;
x = linspace(0.2,size(Cmap,2)/5,size(Cmap,2));
for i=GoodBetas_select
idx=find(idxKmeans==i);
NumberOfCells=length(idx);
subplot(length(GoodBetas_select),3,counter);plot(x,Cmap(i,:));title(num2str(NumberOfCells));axis([0 131 -0.04 0.1]);rectangle('FaceColor','r','Position',[10 -0.04 10 0.005]);rectangle('FaceColor','r','Position',[50 -0.04 10 0.005]);rectangle('FaceColor','r','Position',[90 -0.04 10 0.005]);rectangle('FaceColor','b','Position',[30 -0.04 10 0.005]);rectangle('FaceColor','b','Position',[70 -0.04 10 0.005]);rectangle('FaceColor','b','Position',[110 -0.04 10 0.005]);
subplot(length(GoodBetas_select),3,counter+1);imagesc(zscore(GoodCalcium(idx,:),1,2),[0 4]);colormap hot
subplot(length(GoodBetas_select),3,counter+2);histogram(idx_Plane(idx),[0:20:300]);
counter=counter+3;
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);
counter=1;xplot=floor(sqrt(length(GoodBetas_select)));yplot=ceil(length(GoodBetas_select)/xplot);
for i=GoodBetas_select
idx=find(idxKmeans==i);
NumberOfCells=length(idx);
subplot(length(GoodBetas_select),2,counter);plot(x,Cmap(i,:));title(num2str(NumberOfCells));axis([0 131 -0.04 0.1]);rectangle('FaceColor','r','Position',[11 -0.04 10 0.005]);rectangle('FaceColor','r','Position',[51 -0.04 10 0.005]);rectangle('FaceColor','r','Position',[91 -0.04 10 0.005]);rectangle('FaceColor','b','Position',[31 -0.04 10 0.005]);rectangle('FaceColor','b','Position',[71 -0.04 10 0.005]);rectangle('FaceColor','b','Position',[111 -0.04 10 0.005]);
subplot(length(GoodBetas_select),2,counter+1);imagesc(zscore(GoodCalcium(idx,:),1,2),[0 4]);colormap hot
%subplot(length(GoodBetas_select),3,counter+2);imagesc(zscore(GoodSpikes(idx,:),1,2),[0 4]);colormap hot
counter=counter+2;
end

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
counter=1;
for i=GoodBetas_select
    temp{counter}=find(idxKmeans==i);
    %tempidx=find(idxKmeans==idx);
    %temp{counter}=GoodClusters_goodmembers(counter).idx;
    counter=counter+1;    
end

Start=min(cellfun(@min, temp));Start=find(Numbers<Start,1,'last');
filename=MatFiles(Start).name;

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
            ROIsNb=[ROIsNb ; temp{k}(tempROIsNb)];
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
            name=strcat('_Kmeans_',imagename(4:end));
    imwrite(image3,name,'tif');
    end
    %image3=uint8(image3);

end
clearvars idx i temp tempFileNb fileNb AVG_files filename image counter Numbers image2 image3 k ROI ROIs ROIsNb Start tempROIsNb name imagename tempidx Raster


GoodClustersData=[];
for i=1:length(GoodBetas_select)
    GoodClustersData(i).DF=GoodCalcium(idxKmeans==GoodBetas_select(i),:)*100;
    GoodClustersData(i).Mean=mean(GoodClustersData(i).DF,1);
    GoodClustersData(i).STD=std(GoodClustersData(i).DF,1,1);
end

for i=1:numel(GoodClustersData)
    corr_temp=zeros(size(GoodClustersData(i).ZS,1),1);
    parfor j=1:size(GoodClustersData(i).ZS,1)
        temp=corrcoef(GoodClustersData(i).mean, GoodClustersData(i).ZS(j,:));
        corr_temp(j)=temp(1,2);
    end
    GoodClustersData(i).CorrCoef=corr_temp;
end

GoodClusters_goodmembers=[];Threshold=0.75;
idxKmeans_ZS_goodmembers=idxKmeans_ZS*0;
for i=1:length(GoodBetas_select)
%GoodClusters_goodmembers(i).Spikes=GoodClustersData(i).Spikes(find(GoodClustersData(i).CorrCoef>=0.5),:);
%GoodClusters_goodmembers(i).ZS=zscore(GoodClustersData(i).DF(find(GoodClustersData(i).CorrCoef>=0.5),:),1,2);
GoodClusters_goodmembers(i).ZS=GoodClustersData(i).ZS(find(GoodClustersData(i).CorrCoef>=Threshold),:);
temp=find(idxKmeans_ZS==GoodBetas_select(i));
GoodClusters_goodmembers(i).idx=temp(find(GoodClustersData(i).CorrCoef>=0.5));
GoodClusters_goodmembers(i).mean=mean(GoodClusters_goodmembers(i).ZS,1);
GoodClusters_goodmembers(i).STD=std(GoodClusters_goodmembers(i).ZS,1,1);
idx=find(idxKmeans_ZS==GoodBetas_select(i));
idx=idx(find(GoodClustersData(i).CorrCoef>=Threshold));
idxKmeans_ZS_goodmembers(idx)=GoodBetas_select(i);
%GoodClusters_goodmembers(i).Fish=idx_Fish(idx);
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);x = linspace(0.2,size(Cmap,2)/5,size(Cmap,2));
counter=1;counter2=1;xplot=floor(sqrt(length(GoodBetas_select)));yplot=ceil(length(GoodBetas_select)/xplot);
for i=1:length(GoodBetas_select)  
    if counter==3
        counter=counter+1;
    end
    subplot(3,3,counter);plot(x,GoodClusters_goodmembers(i).mean,'color',colors(counter2,:)/256);axis([0 131 -1 4]);rectangle('FaceColor','r','Position',[11 -1 10 0.25]);rectangle('FaceColor','r','Position',[51 -1 10 0.25]);rectangle('FaceColor','r','Position',[91 -1 10 0.25]);rectangle('FaceColor','b','Position',[31 -1 10 0.25]);rectangle('FaceColor','b','Position',[71 -1 10 0.25]);rectangle('FaceColor','b','Position',[111 -1 10 0.25]);
    counter=counter+1;
    counter2=counter2+1;
end

GoodBetas_select=GoodBetas_selectbackup;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);x = linspace(0.2,size(Cmap,2)/5,size(Cmap,2));
counter=1;counter2=1;xplot=floor(sqrt(length(GoodBetas_select)));yplot=ceil(length(GoodBetas_select)/xplot);
for i=1:length(GoodBetas_select)  
    subplot(length(GoodBetas_select),3,counter);plot(x,GoodClusters_goodmembers(i).mean,'color',colors(counter2,:)/256,'LineWidth',2);axis([0 131 -1 4]);rectangle('FaceColor','r','Position',[11 -1 10 0.25]);rectangle('FaceColor','r','Position',[51 -1 10 0.25]);rectangle('FaceColor','r','Position',[91 -1 10 0.25]);rectangle('FaceColor','b','Position',[31 -1 10 0.25]);rectangle('FaceColor','b','Position',[71 -1 10 0.25]);rectangle('FaceColor','b','Position',[111 -1 10 0.25]);
    subplot(length(GoodBetas_select),3,counter+1);imagesc(GoodClusters_goodmembers(i).ZS,[0 4]);colormap hot
    subplot(length(GoodBetas_select),3,counter+2);histogram(idx_Plane(GoodClusters_goodmembers(i).idx),[0:20:300]);
    counter=counter+3;
    counter2=counter2+1;
end

ClustMean=zeros(length(GoodBetas_select),length(GoodClusters_goodmembers(1).mean));
for i=1:length(GoodBetas_select)  
    ClustMean(i,:)=GoodClusters_goodmembers(i).mean;    
end

PrismTrans=zeros(3*numel(GoodClusters_goodmembers),size(GoodClusters_goodmembers(1).ZS,2));
counter=1;
for i=1:numel(GoodClusters_goodmembers)
    PrismTrans(counter,:)=GoodClusters_goodmembers(i).mean;
    PrismTrans(counter+1,:)=GoodClusters_goodmembers(i).STD;
    PrismTrans(counter+2,:)=ones(1,size(GoodClusters_goodmembers(1).ZS,2))*size(GoodClusters_goodmembers(i).ZS,1);
    counter=counter+3;
end

ZS=zeros(numel(GoodClusters_goodmembers),size(GoodClusters_goodmembers(1).ZS,2));
for i=1:numel(GoodClusters_goodmembers)
    ZS(i,:)=GoodClusters_goodmembers(i).mean;
end


startK=200;
endK=300;
AIC=zeros(1,endK-startK);
BIC=zeros(1,endK-startK);
for i=startK:endK
    [idxKmeans,ClustCenters,sumD]=kmeans(random_5k,i,'Distance','correlation','MaxIter',1000);
    [AICtemp,BICtemp] = Kmeans_AIC(ClustCenters,idxKmeans,sumD);
    AIC(i+1-startK)=AICtemp;
    BIC(i+1-startK)=BICtemp;
end
    
threshold=0.8;
nb_events=zeros(1,size(GoodCalcium,1));
parfor fluotrace=1:size(GoodCalcium,1)
    temp_fluo=GoodCalcium(fluotrace,:);
    corr_results=[];
    for time=1:(size(GoodCalcium,2)-length(GCaMP6))
        corr_temp=corrcoef(temp_fluo(time:time+length(GCaMP6)-1),GCaMP6);
        corr_results(time)=corr_temp(1,2);
    end
    pks = findpeaks(corr_results,'MinPeakProminence',threshold,'MinPeakDistance',length(GCaMP6))
    nb_events(fluotrace)=length(pks);
end

idx_events=find(nb_events>=3);
GoodCalcium_select=GoodCalcium(idx_events,:);
ZS=zscore(GoodCalcium_select,1,2);
imagesc(ZS,[0 4]);colormap hot
save('_Analysis_Flow.mat','-v7.3');
  


Merged_DF=ZS;%   ALL_DF_02=Dataset (neurons x time)
MergedIDX=cell(1,length(Merged_DF)); % Where the IDs of merged ROIs will be stored

MergedIDXb=MergedIDX;
Merged_DFb=Merged_DF;

counter=1;
progressbar % can be skipped if you haven't progressbar installed
while counter<=size(Merged_DF,1) % can also be for counter=1:size(Dataset,1)
    if isnan(Merged_DF(counter,1)) % if the timeserie has been merged, it's set to nan, so this prevents spending time on it
        max_corr=0;
        MergedIDX{counter}=[];
    else
        max_corr=1;    % hack to get the equivalent of a do... while loop (at least one iteration)
        if isempty(MergedIDX{counter})
            MergedIDX{counter}=counter;
        end
    end
    while max_corr>0.85 % standard threshold of Bianco et al
        corr_temp=zeros(size(Merged_DF,1),1);
        parfor j=1:size(Merged_DF,1)    % parallelized, can replace by for loop if not needed
            if j>counter && ~isnan(Merged_DF(j,1))
                temp=corrcoef(Merged_DF(counter,:), Merged_DF(j,:));
                corr_temp(j)=temp(1,2);
            end
        end
        [max_corr idx_max]=nanmax(corr_temp);
        if max_corr>0.85
            if max_corr>0.95 % shortcut to merge in one go all time series with >0.85 correlation, can be skipped or changed
                idx_max=find(corr_temp>0.95);
                for idx=1:length(idx_max)
                    if isempty(MergedIDX{idx_max(idx)})
                        MergedIDX{counter}=[MergedIDX{counter} idx_max(idx)];
                        MergedIDX{idx_max(idx)}=[];
                    else
                        MergedIDX{counter}=[MergedIDX{counter} MergedIDX{idx_max(idx)}];
                        MergedIDX{idx_max(idx)}=[];
                    end
                end
            else
                if isempty(MergedIDX{idx_max})
                    MergedIDX{counter}=[MergedIDX{counter} idx_max]; % merge things one by one, so it's long
                    MergedIDX{idx_max}=[];
                else
                    MergedIDX{counter}=[MergedIDX{counter} MergedIDX{idx_max}];
                    MergedIDX{idx_max}=[];
                end
            end
            Merged_DF(counter,:)=nanmean(ZS(MergedIDX{counter},:),1); %average from dataset, easier than adding one at a time with a factor
            Merged_DF(idx_max,:)=nan;
            corr_temp(idx_max)=nan;
        end
    end
    counter=counter+1;
    progressbar(counter/size(Merged_DF,1)); % can be skipped if you haven't progressbar installed
end
clearvars i j idx_max counter max_corr temp