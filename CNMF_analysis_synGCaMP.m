MatFiles=dir('*matlab*.mat');
name=strcat(MatFiles(1).name);
Calcium=load(name, 'DenoisedTraces');
Calcium=Calcium.DenoisedTraces;
MatFiles(1).number=size(Calcium,1);
Spikes=load(name, 'Spikes');
Spikes=Spikes.Spikes;
Noise=load(name, 'Noise');
Noise=Noise.Noise;
%DF=load(name, 'dFonF');
%DF=DF.dFonF;
Fitness=load(name, 'idx_components');
Fitness=Fitness.idx_components+1;
GoodCalcium=Calcium(Fitness,:);
GoodSpikes=Spikes(Fitness,:);
GoodNoise=Noise(Fitness,:);
%GoodDF=DF(Fitness,:);
MatFiles(1).GoodNumber=length(Fitness);
for i = 2:length(MatFiles)
    name=strcat(MatFiles(i).name);
    C=load(name, 'DenoisedTraces');
    C=C.DenoisedTraces;
    %     if i==3
    %         C=[C(:,1) C(:,1) C(:,1:58)];
    %     end
    S=load(name, 'Spikes');
    S=S.Spikes;
    N=load(name, 'Noise');
    N=N.Noise;
    F=load(name, 'idx_components');
    F=F.idx_components+1;
    %D=load(name, 'dFonF');
    %D=D.dFonF;
    GC=C(F,:);
    GS=S(F,:);
    %GD=D(F,:);
    
    GN=N(F,:);
    if size(C,2)==750
        Calcium=vertcat(Calcium,C(:,51:705));
        Spikes=vertcat(Spikes,S(:,51:705));
        Fitness=horzcat(Fitness,F);
        Noise=vertcat(Noise,N(:,51:705));
        GoodCalcium=vertcat(GoodCalcium,GC(:,51:705));
        GoodNoise=vertcat(GoodNoise,GN(:,51:705));
        %GoodDF=vertcat(GoodDF,GD);
        GoodSpikes=vertcat(GoodSpikes,GS(:,51:705));
    else
        Calcium=vertcat(Calcium,C);
        Spikes=vertcat(Spikes,S);
        Fitness=horzcat(Fitness,F);
        Noise=vertcat(Noise,N);
        GoodCalcium=vertcat(GoodCalcium,GC);
        GoodNoise=vertcat(GoodNoise,GN);
        %GoodDF=vertcat(GoodDF,GD);
        GoodSpikes=vertcat(GoodSpikes,GS);
    end
    
    MatFiles(i).number=size(Calcium,1);
    MatFiles(i).GoodNumber=MatFiles(i-1).GoodNumber+length(F);
end
clearvars GC C S F N name i GS GN N;
ZS=zscore(GoodCalcium,1,2);

parfor i=1:size(ZS,1)
    mdl=stepwiselm(flow',ZS(i,:),'Upper','linear','Intercept',false,'Criterion','bic','verbose',0);
    model_ZS(i).coef=mdl.Coefficients;
    model_ZS(i).MSE=mdl.MSE;
    model_ZS(i).Fitted=mdl.Fitted;
    model_ZS(i).rsquared=mdl.Rsquared.Adjusted;
end

Threshold=0.15;
idx_rsq=find([model_ZS.rsquared]>Threshold);

options = statset('UseParallel',1); [idxKmeans_s Cmap_s]=kmeans(ZS_select,5,'Options',options,'Replicates',5,'MaxIter',1000,'Display','final');
[Model_ZS,GoodBetas]=Test_Regress(Cmap_s,flow,idxKmeans_s,0.2);

Numbers=[MatFiles.GoodNumber];
Numbers=[0 Numbers];
counter=1;
idx_Plane=nan(length(ZS),1);
idx_fish=nan(length(ZS),1);
name=strcat(MatFiles(1).name);
[Plane,~]=regexp(name,'_(\d+)','tokens','match');Plane=str2num(Plane{1}{1});
[Fish,~]=regexp(name,'Fish(\d)_','tokens','match');Fish=str2num(Fish{1}{1});
idx_Plane(1:Numbers(2))=Plane;
idx_fish(1:Numbers(2))=Fish;
for i=2:length(MatFiles)      
    name=strcat(MatFiles(i).name);    
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
for i=GoodBetas  
    idx=find(idxKmeans_s==i);
    GoodClusterData(counter).mean=mean(ZS_select(idx,:),1);
    GoodClusterData(counter).ZS=ZS_select(idx,:);
    GoodClusterData(counter).planes=idx_Plane(idx_rsq(idx));
    GoodClusterData(counter).fish=idx_fish(idx_rsq(idx));
    subplot(length(GoodBetas),4,counter2);plot(GoodClusterData(counter).mean,'color',colors(counter,:));title(num2str(length(idx))),xlim([0 size(ZS,2)]),ylim([-1 4])    
    subplot(length(GoodBetas),4,counter2+1);imagesc(GoodClusterData(counter).ZS, [0 4]); colormap hot;title(num2str(length(idx)))    
    subplot(length(GoodBetas),4,counter2+2);histogram(GoodClusterData(counter).planes);
    subplot(length(GoodBetas),4,counter2+3);histogram(GoodClusterData(counter).fish);
    %subplot(length(GoodBetas_select),4,counter2+3);bar([sum(GoodClusterData(counter).State==3) sum(GoodClusterData(counter).State==4)]);set(gca,'xticklabel',name);
    counter2=counter2+4;
    counter=counter+1;
end

for i=1:numel(GoodClusterData)
    corr_temp=zeros(size(GoodClusterData(i).ZS,1),1);
    parfor j=1:size(GoodClusterData(i).ZS,1)
        temp=corrcoef(GoodClusterData(i).mean, GoodClusterData(i).ZS(j,:));
        corr_temp(j)=temp(1,2);
    end
    GoodClusterData(i).CorrCoef=corr_temp;
end

GoodClusters_goodmembers=[];Threshold=0.5;
idxKmeans_ZS_goodmembers=idxKmeans_s*0;
for i=1:length(GoodBetas)
%GoodClusters_goodmembers(i).Spikes=GoodClusterData(i).Spikes(find(GoodClusterData(i).CorrCoef>=0.5),:);
%GoodClusters_goodmembers(i).ZS=zscore(GoodClusterData(i).DF(find(GoodClusterData(i).CorrCoef>=0.5),:),1,2);
GoodClusters_goodmembers(i).ZS=GoodClusterData(i).ZS(find(GoodClusterData(i).CorrCoef>=Threshold),:);
temp=find(idxKmeans_s==GoodBetas(i));
GoodClusters_goodmembers(i).idx=temp(find(GoodClusterData(i).CorrCoef>=Threshold));
GoodClusters_goodmembers(i).mean=mean(GoodClusters_goodmembers(i).ZS,1);
GoodClusters_goodmembers(i).STD=std(GoodClusters_goodmembers(i).ZS,1,1);
idx=find(idxKmeans_s==GoodBetas(i));
idx=idx(find(GoodClusterData(i).CorrCoef>=Threshold));
idxKmeans_ZS_goodmembers(idx)=GoodBetas(i);
%GoodClusters_goodmembers(i).Fish=idx_Fish(idx);
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
%GoodBetas_select=[1 2];
idxKmeans_final=zeros(size(idxKmeans));
for i=GoodBetas
    idxKmeans_final(idx_rsq(find(idxKmeans_ZS_goodmembers==i)))=i;
    temp{counter}=find(idxKmeans_final==i);
    %tempidx=find(idxKmeans==idx);
    %temp{counter}=Select(tempidx)';
    counter=counter+1;    
end

Start=min(cellfun(@min, temp));Start=find(Numbers<Start,1,'last');
filename=MatFiles(Start).name;

%colors = [1,0,0;0,1,0;0,0,1;1,0.600000000000000,0;1,0,0.7];
colors = distinguishable_colors(length(GoodBetas),[1 1 1; 0 0 0]);
%colors = colors*256;
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
        imagename=regexp(filename,'_output_analysis_matlab.mat','split');
        imagename=strcat('AVG_',imagename{1},'.tif');
        image=double(imread(imagename));image(image==max(max(image)))=0;image=image/max(prctile(image,95));image=image*100;
        image=uint8(image);
        image2=zeros(size(image(:,:,1)));
        image3=repmat(image,1,1,3);
        ROIs=All_ROIs{idx};       
        ROIsNb=ROIsNb-Numbers(idx);
        ROIs=ROIs(:,ROIsNb);
        for k = 1 : size(ROIs,2)
            image2=zeros(size(image(:,:,1)));
            ROI=full(reshape(ROIs(:,k),size(image,1),size(image,2)));            
            image2=image2+ROI;image2=(image2/max(max(image2)))*120;image2=uint8(image2);
            for j=1:3
                image3(:,:,j)=image3(:,:,j)+image2*colors(ClusterNb(k),j);
            end
        end
        %image3(:,:,3)=image;
    end
    %image3=uint8(image3);
    name=strcat('Kmeans_good',imagename(4:end));
    imwrite(image3,name,'tif');
end

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

