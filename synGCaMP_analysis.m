filelist=dir('fish*.tif');
for File=1:length(filelist)
    [mixedsig, mixedfilters, CovEvals, covtrace, movm, movtm] = CellsortPCA(filelist(File).name,[],100,1,'C:\Temp\CROP');
    [ica_sig, ica_filters, ica_A, numiter] = CellsortICA(mixedsig, mixedfilters, CovEvals, [1:size(mixedsig,1)],  0.8,size(mixedsig,1),randn(size(mixedsig,1), size(mixedsig,1)),1e-6,50000);
    %PCA_ICA_results(File).ica_sig=ica_sig;
    %PCA_ICA_results(File).ica_filters=ica_filters;
    %[PCA_ICA_results(File).ica_segments, PCA_ICA_results(File).segmentlabel, PCA_ICA_results(File).segcentroid] = CellsortSegmentation(ica_filters, 2, 2, 20, 0)
    [ica_segments, segmentlabel, PCA_ICA_results(File).segcentroid] = CellsortSegmentation(ica_filters, 2, 2, 20, 0);    
    PCA_ICA_results(File).cell_sig = CellsortApplyFilter(filelist(File).name, ica_segments);
end
clearvars mixedsig ica_sig ica_filters ica_A numiter mixedsig mixedfilters CovEvals covtrace movm movtm
 
AllSegTraces=[];
idx=1;AllSegTraces=PCA_ICA_results(idx).cell_sig;
for idx=2:length(PCA_ICA_results)
    if size(PCA_ICA_results(idx).cell_sig,2)==750
        AllSegTraces=vertcat(AllSegTraces,PCA_ICA_results(idx).cell_sig(:,51:705));
    else
        AllSegTraces=vertcat(AllSegTraces,PCA_ICA_results(idx).cell_sig);
    end
end
DF=DeltaF2(AllSegTraces,11,21);
ZS=zscore(AllSegTraces,1,2);

GMModels = {};maxCompNb=80;
options = statset('MaxIter',500);
progressbar
for k = 76:maxCompNb
    GMModels{k} = fitgmdist(DF,k,'Options',options,'CovarianceType','diagonal','Options',options, 'Regularize', 1e-5);
    BIC(k)= GMModels{k}.BIC;
    GMModels{k}.BIC
    progressbar(k/maxCompNb);
end
[minBIC,numComponents] = min(BIC);
numComponents
BIC_smooth=smooth(BIC');
figure;plot(BIC_smooth);


options = statset('UseParallel',1); [idxKmeans Cmap]=kmeans(DF,50,'Options',options,'Distance','correlation','Replicates',5,'MaxIter',1000,'Display','final');
Threshold=0.2;
[Model_DF,GoodBetas]=Test_Regress(Cmap,flow,idxKmeans,Threshold);

parfor i=1:size(DF,1)
    mdl=stepwiselm(flow',DF(i,:),'Upper','linear','Intercept',false,'Criterion','bic','verbose',0);
    model_DF(i).coef=mdl.Coefficients;
    model_DF(i).MSE=mdl.MSE;
    model_DF(i).Fitted=mdl.Fitted;
    model_DF(i).rsquared=mdl.Rsquared.Adjusted;
end

parfor i=1:size(ZS,1)
    mdl=stepwiselm(flow',ZS(i,:),'Upper','linear','Intercept',false,'Criterion','bic','verbose',0);
    model_ZS(i).coef=mdl.Coefficients;
    model_ZS(i).MSE=mdl.MSE;
    model_ZS(i).Fitted=mdl.Fitted;
    model_ZS(i).rsquared=mdl.Rsquared.Adjusted;
end


Threshold=0.2;
idx_rsq=find([model_DF.rsquared]>Threshold);
figure;imagesc(DF(idx_rsq,:),[0 0.1]);colormap hot
options = statset('UseParallel',1); [idxKmeans2 Cmap2]=kmeans(DF(idx_rsq,:),5,'Options',options,'Distance','correlation','Replicates',5,'MaxIter',1000,'Display','final');
[Model_DF2,GoodBetas2]=Test_Regress(Cmap2,flow,idxKmeans2,Threshold);

Threshold=0.2;
idx_rsq_ZS=find([model_ZS.rsquared]>Threshold);
figure;imagesc(ZS(idx_rsq_ZS,:),[0 4]);colormap hot
options = statset('UseParallel',1); [idxKmeans_ZS Cmap_ZS]=kmeans(ZS(idx_rsq_ZS,:),5,'Options',options,'Distance','correlation','Replicates',5,'MaxIter',1000,'Display','final');
[Model_ZS,GoodBetas_ZS]=Test_Regress(Cmap_ZS,flow,idxKmeans_ZS,Threshold);

DF_rsq=DF(idx_rsq,:);
Merged_DF=DF(idx_rsq,:);%   ALL_DF_02=Dataset (neurons x time)
MergedIDX=cell(1,length(Merged_DF)); % Where the IDs of merged ROIs will be stored

Merged_DFb=Merged_DF;%   ALL_DF_02=Dataset (neurons x time)
MergedIDXb=MergedIDX; % Where the IDs of merged ROIs will be stored

counter=1;Threshold=0.75
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
    while max_corr>Threshold % standard threshold of Bianco et al
        corr_temp=zeros(size(Merged_DF,1),1);
        parfor j=1:size(Merged_DF,1)    % parallelized, can replace by for loop if not needed
            if j>counter && ~isnan(Merged_DF(j,1))
                temp=corrcoef(Merged_DF(counter,:), Merged_DF(j,:));
                corr_temp(j)=temp(1,2);
            end
        end
        [max_corr idx_max]=nanmax(corr_temp);
        if max_corr>Threshold
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
            Merged_DF(counter,:)=nanmean(DF_rsq(MergedIDX{counter},:),1); %average from dataset, easier than adding one at a time with a factor
            Merged_DF(idx_max,:)=nan;
            corr_temp(idx_max)=nan;
        end
    end
    counter=counter+1;
    progressbar(counter/size(Merged_DF,1)); % can be skipped if you haven't progressbar installed
end
clearvars i j idx_max counter max_corr temp

temp=find(numMerged>20);
MergedComps=[];
counter=1;
for i=temp'
    MergedComps(counter).mean=mean(DF_rsq(MergedIDX{i},:)*100,1);
    MergedComps(counter).STD=std(DF_rsq(MergedIDX{i},:)*100,1,1);
    MergedComps(counter).NB=ones(1,size(DF_rsq,2))*numel(MergedIDX{i});    
    counter=counter+1;
end

TempComps=MergedComps(1).mean;
for i=2:length(temp);
    TempComps=[TempComps; MergedComps(i).mean];
end

max_DF=max(DF_rsq,[],2);
min_DF=min(DF_rsq,[],2);
idx_Threshold_5to200=find(max_DF>0.05 & max_DF<2 & min_DF>-2);
idx_Threshold_5to200_bool=(max_DF>0.05 & max_DF<2 & min_DF>-2);
Select_DF_rsq=DF_rsq(idx_Threshold_5to200,:);
Threshold=0.1;

ClustMean=zeros(max(c),size(DF,2));
for i=1:max(c)
    ClustMean(i,:)=mean(Select_DF_rsq(find(c==i),:),1);
end

GMModels = {};
options = statset('MaxIter',500);
for k = 48:50
    GMModels{k} = fitgmdist(Select_DF_rsq,k,'Options',options,'CovarianceType','diagonal','Options',options, 'Regularize', 1e-5);
    BIC(k)= GMModels{k}.BIC;
end
[minBIC,numComponents] = min(BIC);
numComponents
BIC_smooth=smooth(BIC');
figure;plot(BIC_smooth);


for File=1:length(filelist)
    PCA_ICA_results(File).Numbers=size(PCA_ICA_results(File).cell_sig,1);
end

Numbers=cumsum([PCA_ICA_result.Numbers]);
Numbers=Number;
counter=1;
idx_Plane=nan(length(DF),1);
idx_fish=nan(length(DF),1);
name=strcat(filelist(1).name);
[Plane,~]=regexp(name,'_(\d+)','tokens','match');Plane=str2num(Plane{1}{1});
[Fish,~]=regexp(name,'Fish(\d)_','tokens','match');Fish=str2num(Fish{1}{1});
idx_Plane(1:Numbers(2))=Plane;
idx_fish(1:Numbers(2))=Fish;
for i=2:length(filelist)    
    name=strcat(filelist(i).name);    
    [Plane,~]=regexp(name,'_(\d+)','tokens','match');
    if isempty(Plane)
        if strfind(name,'Torus')
            Plane=150;
        elseif strfind(name,'nr 11')
            Plane=15;
        end
    else
        Plane=str2num(Plane{1}{1});
    end
    [Fish,~]=regexp(name,'Fish(\d)_','tokens','match');Fish=str2num(Fish{1}{1});
    idx_Plane(Numbers(i):Numbers(i+1))=Plane;
    idx_fish(Numbers(i):Numbers(i+1))=Fish;
end

counter=1;counter2=1;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
GoodClusterData=[];
%name=['ELO';'CLO'];
x = linspace(0.2,size(ZS,2)/5,size(ZS,2));
for i=GoodBetas2  
    idx=find(idxKmeans_ZS==i);
    GoodClusterData(counter).mean=mean(ZS_rsq(idx,:),1);
    GoodClusterData(counter).ZS=ZS_rsq(idx,:);
    GoodClusterData(counter).planes=idx_Plane(idx_rsq_ZS(idx));
    GoodClusterData(counter).fish=idx_fish(idx_rsq_ZS(idx));
    subplot(length(GoodBetas2),4,counter2);plot(x,GoodClusterData(counter).mean);title(num2str(length(idx)));ylim([-2 4]);xlim([0 130]);rectangle('FaceColor','r','Position',[10 -2 10 0.5]);rectangle('FaceColor','r','Position',[50 -2 10 0.5]);rectangle('FaceColor','r','Position',[90 -2 10 0.5]);rectangle('FaceColor','b','Position',[30 -2 10 0.5]);rectangle('FaceColor','b','Position',[70 -2 10 0.5]);rectangle('FaceColor','b','Position',[110 -2 10 0.5]);    
    subplot(length(GoodBetas2),4,counter2+1);imagesc(GoodClusterData(counter).ZS, [0 4]); colormap hot;title(num2str(length(idx)))    
    subplot(length(GoodBetas2),4,counter2+2);histogram(GoodClusterData(counter).planes);
    subplot(length(GoodBetas2),4,counter2+3);histogram(GoodClusterData(counter).fish);
    %subplot(length(GoodBetas_select),4,counter2+3);bar([sum(GoodClusterData(counter).State==3) sum(GoodClusterData(counter).State==4)]);set(gca,'xticklabel',name);
    counter2=counter2+4;
    counter=counter+1;
end



for i=1:numel(GoodClusterData)
    corr_temp=zeros(size(GoodClusterData(i).DF,1),1);
    parfor j=1:size(GoodClusterData(i).DF,1)
        temp=corrcoef(GoodClusterData(i).mean, GoodClusterData(i).DF(j,:));
        corr_temp(j)=temp(1,2);
    end
    GoodClusterData(i).CorrCoef=corr_temp;
end

GoodClusters_goodmembers=[];Threshold=0.5;
idxKmeans_ZS_goodmembers=idxKmeans_ZS*0;
for i=1:length(GoodBetas_ZS)
%GoodClusters_goodmembers(i).Spikes=GoodClusterData(i).Spikes(find(GoodClusterData(i).CorrCoef>=0.5),:);
%GoodClusters_goodmembers(i).ZS=zscore(GoodClusterData(i).DF(find(GoodClusterData(i).CorrCoef>=0.5),:),1,2);
GoodClusters_goodmembers(i).ZS=GoodClusterData(i).ZS(find(GoodClusterData(i).CorrCoef>=Threshold),:);
temp=find(idxKmeans_ZS==GoodBetas_ZS(i));
GoodClusters_goodmembers(i).idx=temp(find(GoodClusterData(i).CorrCoef>=Threshold));
GoodClusters_goodmembers(i).mean=mean(GoodClusters_goodmembers(i).ZS,1);
GoodClusters_goodmembers(i).STD=std(GoodClusters_goodmembers(i).ZS,1,1);
idx=find(idxKmeans_ZS==GoodBetas_ZS(i));
idx=idx(find(GoodClusterData(i).CorrCoef>=Threshold));
idxKmeans_ZS_goodmembers(idx)=GoodBetas_ZS(i);
%GoodClusters_goodmembers(i).Fish=idx_Fish(idx);
end

colors = [0         0    1.0000
         0    1    0
    1.0000         0         0
    1.0000    0.1034    0.7241
    1.0000    0.5000    0.3000
         0    0.7000    0.2000
    0.5000    0.5000         0
         0    0.5000    0.5000];
colors = colors*256;

GoodBetas=GoodBetas_ZS;
GoodBetas_select=GoodBetas2;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);x = linspace(0.2,size(ZS,2)/5,size(Z,2));
counter=1;counter2=1;xplot=floor(sqrt(length(GoodBetas_select)));yplot=ceil(length(GoodBetas_select)/xplot);
for i=1:length(GoodBetas_select)  
    subplot(length(GoodBetas_select),3,counter);plot(x,mean(zscore(GoodClusters_goodmembers(i).DF,1,2),1),'color',colors(counter2,:)/256,'LineWidth',2);axis([0 131 -1 4]);rectangle('FaceColor','r','Position',[11 -1 10 0.25]);rectangle('FaceColor','r','Position',[51 -1 10 0.25]);rectangle('FaceColor','r','Position',[91 -1 10 0.25]);rectangle('FaceColor','b','Position',[31 -1 10 0.25]);rectangle('FaceColor','b','Position',[71 -1 10 0.25]);rectangle('FaceColor','b','Position',[111 -1 10 0.25]);
    subplot(length(GoodBetas_select),3,counter+1);imagesc(zscore(GoodClusters_goodmembers(i).DF,1,2),[0 4]);colormap hot
    subplot(length(GoodBetas_select),3,counter+2);histogram(idx_Plane(GoodClusters_goodmembers(i).idx));
    counter=counter+3;
    counter2=counter2+1;
end

idxKmeans_final=zeros(1,size(DF,1));
idxKmeans_final(idx_rsq)=idxKmeans2;

File=1
load(strcat(filelist(File).name(1:end-3),'mat'), 'PCA_ICA_results');
Traces=PCA_ICA_results.Cell_sig;AllSegTraces=Traces;
Numbers=size(PCA_ICA_results.Cell_sig,1);
for File=2:length(filelist)
    load(strcat(filelist(File).name(1:end-3),'mat'), 'PCA_ICA_results');
    Traces=PCA_ICA_results.Cell_sig;
    Number=size(PCA_ICA_results.Cell_sig,1);
    Numbers=[Numbers Number];
    if size(Traces,2)==750
        AllSegTraces=vertcat(AllSegTraces,Traces(:,51:705));
    else
        AllSegTraces=vertcat(AllSegTraces,Traces);
    end
end

colors = distinguishable_colors(length(GoodBetas_ZS),[1 1 1; 0 0 0]);
colors = colors(GoodBetas,:);
Number=[1 Number];
idxKmeans_final(idx_rsq_ZS)=idxKmeans_ZS_goodmembers;
for File=1:length(filelist)
    load(strcat(filelist(File).name(1:end-3),'mat'), 'PCA_ICA_results');
    ROIs=PCA_ICA_results.ROIs;
    idx_temp=idxKmeans_final(Number(File):Number(File+1));
    imagename=strcat('AVG_',filelist(File).name);
    image=double(imread(imagename));image(image==max(max(image)))=0;image=image/max(prctile(image,95));image=image*64;
    image=uint8(image);    
    image3=repmat(image,1,1,3);
    for i=GoodBetas_ZS
        idx_ROI=find(idx_temp==i);
        image2=squeeze(sum(ROIs(idx_ROI,:,:),1));
        image2=(image2/max(max(image2)))*150;image2=uint8(image2);
        for j=1:3
            image3(:,:,j)=image3(:,:,j)+image2*colors(i,j);
        end
    end
    name=strcat('Kmeans_good',imagename(4:end));
    imwrite(image3,name,'tif');
end

GoodBetas=[1 3 5 2 4];
x = linspace(0.2,size(ZS,2)/5,size(ZS,2));
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;xplot=floor(sqrt(length(GoodBetas)));yplot=ceil(length(GoodBetas)/xplot);
for i=GoodBetas    
    NumberOfCells=length(find(idxKmeans_final==i));
    %NumberOfCells=length(find(idxKmeans_ZS_goodmembers==i));
    %subplot(5,1,counter);plot(x,Cmap_s(i,:),'color',colors(counter,:));title(num2str(NumberOfCells))
    subplot(5,1,counter);plot(x,GoodClusters_goodmembers(i).mean,'color',colors(counter,:),'LineWidth',2);axis([0 131 -2 4]);rectangle('FaceColor','r','Position',[11 -2 10 0.25]);rectangle('FaceColor','r','Position',[51 -2 10 0.25]);rectangle('FaceColor','r','Position',[91 -2 10 0.25]);rectangle('FaceColor','b','Position',[31 -2 10 0.25]);rectangle('FaceColor','b','Position',[71 -2 10 0.25]);rectangle('FaceColor','b','Position',[111 -2 10 0.25]);;title(num2str(NumberOfCells))
%    xlim([0 size(Cmap_ZS,2)])
    counter=counter+1;
end
clearvars idx i temp tempFileNb fileNb AVG_files filename image counter Numbers image2 image3 k ROI ROIs ROIsNb Start tempROIsNb name imagename tempidx Raster

