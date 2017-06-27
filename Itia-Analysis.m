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
Noise=vertcat(Noise,N);
Calcium=vertcat(Calcium,C);
%DF=vertcat(DF,D);
Spikes=vertcat(Spikes,S);
Fitness=horzcat(Fitness,F);
GoodCalcium=vertcat(GoodCalcium,GC);
%GoodDF=vertcat(GoodDF,GD);
GoodSpikes=vertcat(GoodSpikes,GS);
MatFiles(i).number=size(Calcium,1);
MatFiles(i).GoodNumber=MatFiles(i-1).GoodNumber+length(F);
end
clearvars GC C S F N name i GS;\
ZS=zscore(GoodCalcium,1,2);

Vestibular=zeros(2,size(GoodCalcium,2));
GCaMP6=[-0.104392135015146,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
%GCaMP6=[0.000256990000000000;0.00850739000000000;0.0654158300000000;0.0784609000000000;0.0764130100000000;0.0665958600000000;0.0579028900000000;0.0467942900000000;0.0232079800000000;0.0144564400000000;0.00695772000000000;0.00526551000000000;0.00299500000000000;0.00198520000000000;0.00128512000000000;0.00134175000000000;0.000403170000000000;0];
%back=[14 114 214];
%back_off=[24 124 224];
GCaMP6=interp(GCaMP6,4);
back=[30 230 430];
%back=[30/5 230/5 430/5];
back_off=[50 250 450];
%back_off=[50/5 250/5 450/5];
Vestibular(1,back(1):back(1)+size(GCaMP6,1)-1)=GCaMP6';
Vestibular(1,back(2):back(2)+size(GCaMP6,1)-1)=GCaMP6';
%Vestibular(1,back(3):back(3)+size(GCaMP6,1)-1)=GCaMP6';
Vestibular(1,back(3):back(3)+size(GCaMP6,1)-1)=GCaMP6';
Vestibular(2,back_off(1):back_off(1)+size(GCaMP6,1)-1)=GCaMP6';
Vestibular(2,back_off(2):back_off(2)+size(GCaMP6,1)-1)=GCaMP6';
%Vestibular(2,back_off(3):back_off(3)+size(GCaMP6,1)-1)=GCaMP6';
Vestibular(2,back_off(3):back_off(3)+size(GCaMP6,1)-1)=GCaMP6';
%Vestibular=Vestibular(:,1:60);
clearvars GCaMP6 back back_off fwd fwd_off;





options = statset('UseParallel',1); [idxKmeans_ZS Cmap_ZS]=kmeans(ZS,100,'Options',options,'Distance','correlation','Replicates',5,'MaxIter',1000,'Display','final');

[Model_GC,GoodBetas_GC]=Test_Regress(Cmap_ZS,Vestibular,idxKmeans_ZS,0.5);
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);

GoodBetas_select=GoodBetas_GC;
% GoodBetas_select_DF=GoodBetas_DF([1 3 5 6]);
%GoodBetas_select=GoodBetas_GC([1 3 5 7 8 9 10 11 12 15]);


% x = linspace(1,size(Cmap_ZS,2),size(Cmap_ZS,2));
% Fighandle=figure;
% set(Fighandle, 'Position', [100, 100, 1300, 900]);
% counter=1;xplot=floor(sqrt(length(GoodBetas_select)));yplot=ceil(length(GoodBetas_select)/xplot);
% for i=GoodBetas_select    
%     NumberOfCells=length(find(idxKmeans_ZS==i));
%     %subplot(5,1,counter);plot(x,Cmap_ZS(i,:),x,Model_DF(i).Fitted);title(num2str(NumberOfCells))
%     subplot(xplot,yplot,counter);plot(Cmap_ZS(i,:));title(num2str(NumberOfCells))
%     %subplot(xplot,yplot,counter);imagesc(DF(find(idxKmeans_ZS==i),:),[0 0.3]);colormap hot;title(num2str(NumberOfCells))
%     xlim([0 size(Cmap_ZS,2)])
%     counter=counter+1;
% end

Numbers=[0 [MatFiles.GoodNumber]];
counter=1;
idx_Plane_CTRL=nan(length(GoodCalcium),1);
name=strcat(MatFiles(1).name);
%[Plane,~]=regexp(name,'_(\d+)mA_','tokens','match');Plane=str2num(Plane{1}{1});
[Plane,~]=regexp(name,'_(\d+)','tokens','match');Plane=str2num(Plane{1}{1});
if strfind(name,'CLO');
    CENTRE=strfind(name,'CLO');CENTRE=double(~isempty(CENTRE));
elseif strfind(name,'DLO');    
    CENTRE=strfind(name,'DLO');CENTRE=double(~isempty(CENTRE))*2;
elseif strfind(name,'ELO');
    CENTRE=strfind(name,'ELO');CENTRE=double(~isempty(CENTRE))*3;
elseif strfind(name,'OLO');
    CENTRE=strfind(name,'OLO');CENTRE=double(~isempty(CENTRE))*4;    
end
idx_Plane_CTRL(1:Numbers(2))=Plane;
idx_State(1:Numbers(2))=CENTRE;
for i=2:length(MatFiles)
    %[Fish,~]=regexp(files{i},'(\d+)_','tokens','match');Fish=str2num(Fish{1}{1});
    name=strcat(MatFiles(i).name);
    %[Plane,~]=regexp(name,'_(\d+)mA_','tokens','match');Plane=str2num(Plane{1}{1});
    [Plane,~]=regexp(name,'_(\d+)','tokens','match');Plane=str2num(Plane{1}{1});
    if strfind(name,'mA')
        Plane=Plane/10;
    end
    if strfind(name,'CLO');
        CENTRE=strfind(name,'CLO');CENTRE=double(~isempty(CENTRE));
    elseif strfind(name,'DLO');
        CENTRE=strfind(name,'DLO');CENTRE=double(~isempty(CENTRE))*2;
    elseif strfind(name,'ELO');
        CENTRE=strfind(name,'ELO');CENTRE=double(~isempty(CENTRE))*3;
    elseif strfind(name,'OLO');
        CENTRE=strfind(name,'OLO');CENTRE=double(~isempty(CENTRE))*4;
    end
    idx_Plane_CTRL(Numbers(i):Numbers(i+1))=Plane;
    idx_State(Numbers(i):Numbers(i+1))=CENTRE;
end

% load('Itia2-dict_allindex_Maskb.mat');
% index=[1; index];
% for i=1:length(files)
%     name=files{i};
%     [Plane,~]=regexp(name,'fish1_(\d+)','tokens','match');Plane=str2num(Plane{1}{1});    
%     if strfind(name,'CENTRE');
%         CENTRE=strfind(name,'CENTRE');CENTRE=double(~isempty(CENTRE));
%     elseif strfind(name,'OUT');
%         strfind(name,'OUT');
%         CENTRE=strfind(name,'OUT');CENTRE=double(~isempty(CENTRE))*2;
%     elseif strfind(name,'ELO');
%         CENTRE=strfind(name,'ELO');CENTRE=double(~isempty(CENTRE))*3;
%     end
%     idx_Plane_Morph(index(i):index(i+1))=Plane;
%     idx_State_Morph(index(i):index(i+1))=CENTRE;
% end
% idx_Plane_Morph=idx_Plane_Morph(idx_Threshold_5to200);
% idx_State_Morph=idx_State_Morph(idx_Threshold_5to200);

%GoodBetas_select=GoodBetas_GM;
counter=1;counter2=1;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
name=['CLO';'DLO';'ELO';'OLO'];
GoodClusterData=[];
%name=['ELO';'CLO'];
for i=GoodBetas_select  
    idx=find(idxKmeans_ZS==i);
    GoodClusterData(counter).mean=mean(ZS(idx,:),1);
    GoodClusterData(counter).ZS=ZS(idx,:);
    GoodClusterData(counter).planes=idx_Plane_CTRL(idx);
    GoodClusterData(counter).State=idx_State(idx);
    subplot(length(GoodBetas_select),4,counter2);plot(GoodClusterData(counter).mean);title(num2str(length(idx)))    
    subplot(length(GoodBetas_select),4,counter2+1);imagesc(GoodClusterData(counter).ZS, [0 4]); colormap hot;title(num2str(length(idx)))    
    subplot(length(GoodBetas_select),4,counter2+2);histogram(GoodClusterData(counter).planes);
    subplot(length(GoodBetas_select),4,counter2+3);bar([sum(GoodClusterData(counter).State==1) sum(GoodClusterData(counter).State==2) sum(GoodClusterData(counter).State==3) sum(GoodClusterData(counter).State==4)]);set(gca,'xticklabel',name);
    %subplot(length(GoodBetas_select),4,counter2+3);bar([sum(GoodClusterData(counter).State==3) sum(GoodClusterData(counter).State==4)]);set(gca,'xticklabel',name);
    counter2=counter2+4;
    counter=counter+1;
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
imagesc(ZS, [0 4]); colormap hot;

imagesc(Cmap_GM, [0 4]); colormap hot;
% counter=1;counter2=1;
% Fighandle=figure;
% set(Fighandle, 'Position', [100, 100, 1300, 900]);
% for i=GoodBetas_select_DF  
%     idx=find(idxKmeans_DF==i);
%     GoodClusterData(counter).mean=mean(DF_5to200(idx,:),1);
%     GoodClusterData(counter).DF_5to200=DF_5to200(idx,:);
%     GoodClusterData(counter).planes=idx_Plane_Morph(idx);
%     GoodClusterData(counter).State=idx_State_Morph(idx);
%     subplot(length(GoodBetas_select),3,counter2);plot(GoodClusterData(counter).mean);title(num2str(length(idx)))
%     subplot(length(GoodBetas_select),3,counter2+1);histogram(GoodClusterData(counter).planes,[0:50:250]);
%     subplot(length(GoodBetas_select),3,counter2+2);bar([sum(GoodClusterData(counter).State==1) sum(GoodClusterData(counter).State==2) sum(GoodClusterData(counter).State==3)]);
%     counter2=counter2+3;
%     counter=counter+1;
% end

ModelResultsSeg_ZS=[];
parfor i=1:length(ZS)
    %mdl=stepwiselm(Regressor',DataSet(i,:),'linear','Criterion','adjrsquared','Intercept',false,'Upper','interactions','Verbose',0);
    mdl=stepwiselm(Vestibular',ZS(i,:),'linear','Criterion','adjrsquared','Upper','linear','Verbose',0);
    ModelResultsSeg_ZS(i).coef=mdl.Coefficients;
    ModelResultsSeg_ZS(i).MSE=mdl.MSE;
    ModelResultsSeg_ZS(i).Fitted=mdl.Fitted;
    ModelResultsSeg_ZS(i).rsquared=mdl.Rsquared.Adjusted;
end
idx_rsq_ZS=find([ModelResultsSeg_ZS.rsquared]>0.3);
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);
imagesc(ZS(idx_rsq_ZS,:),[0 5]);colormap hot

for i=idx_rsq_ZS  
    idx=find(idxKmeans_ZS==i);
    GoodClusterData_rsq(counter).mean=mean(ZS(idx,:),1);
    GoodClusterData_rsq(counter).ZS=ZS(idx,:);
    GoodClusterData_rsq(counter).planes=idx_Plane_CTRL(idx);
    GoodClusterData_rsq(counter).State=idx_State(idx);
    %subplot(length(GoodBetas_select),4,counter2);plot(GoodClusterData_rsq(counter).mean);title(num2str(length(idx)))    
    %subplot(length(GoodBetas_select),4,counter2+1);imagesc(GoodClusterData_rsq(counter).ZS, [0 4]); colormap hot;title(num2str(length(idx)))    
    %subplot(length(GoodBetas_select),4,counter2+2);histogram(GoodClusterData_rsq(counter).planes);
    %subplot(length(GoodBetas_select),4,counter2+3);bar([sum(GoodClusterData_rsq(counter).State==1) sum(GoodClusterData_rsq(counter).State==2) sum(GoodClusterData_rsq(counter).State==3) sum(GoodClusterData_rsq(counter).State==4)]);set(gca,'xticklabel',name);
    %subplot(length(GoodBetas_select),4,counter2+3);bar([sum(GoodClusterData(counter).State==3) sum(GoodClusterData(counter).State==4)]);set(gca,'xticklabel',name);
    counter2=counter2+4;
    counter=counter+1;
end
figure;bar([sum([GoodClusterData_rsq.State]==1) sum([GoodClusterData_rsq.State]==2) sum([GoodClusterData_rsq.State]==3) sum([GoodClusterData_rsq.State]==4)]);set(gca,'xticklabel',name);

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
    temp{counter}=find(idxKmeans_ZS==i);
    %tempidx=find(idxKmeans==idx);
    %temp{counter}=Select(tempidx)';
    counter=counter+1;    
end

Start=min(cellfun(@min, temp));Start=find(Numbers<Start,1,'last');
filename=MatFiles(Start).name;

%colors = [1,0,0;0,1,0;0,0,1;1,0.600000000000000,0;1,0,0.7];
colors = distinguishable_colors(length(GoodBetas_select),[1 1 1; 0 0 0]);
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
        image=double(imread(imagename));image(image==max(max(image)))=0;image=image/max(prctile(image,95));image=image*64;
        image=uint8(image);
        image2=zeros(size(image(:,:,1)));
        image3=repmat(image,1,1,3);
        ROIs=All_ROIs{idx};       
        ROIsNb=ROIsNb-Numbers(idx);
        ROIs=ROIs(:,ROIsNb);
        for k = 1 : size(ROIs,2)
            image2=zeros(size(image(:,:,1)));
            ROI=full(reshape(ROIs(:,k),size(image,1),size(image,2)));            
            image2=image2+ROI;image2=(image2/max(max(image2)))*128;image2=uint8(image2);
            for j=1:3
                image3(:,:,j)=image3(:,:,j)+image2*colors(ClusterNb(k),j);
            end
        end
        %image3(:,:,3)=image;
    end
    %image3=uint8(image3);
    name=strcat('Kmeans',imagename(4:end));
    imwrite(image3,name,'tif');
end

x = linspace(1,size(Cmap_ZS,2),size(Cmap_ZS,2));
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;xplot=floor(sqrt(length(GoodBetas_select)));yplot=ceil(length(GoodBetas_select)/xplot);
for i=GoodBetas_select    
    NumberOfCells=length(find(idxKmeans_ZS==i));
    %subplot(5,1,counter);plot(x,Cmap_ZS(i,:),x,Model_DF(i).Fitted);title(num2str(NumberOfCells))
    subplot(xplot,yplot,counter);plot(Cmap_ZS(i,:),'color',colors(counter,:));title(num2str(NumberOfCells))
    xlim([0 size(Cmap_ZS,2)])
    counter=counter+1;
end
clearvars idx i temp tempFileNb fileNb AVG_files filename image counter Numbers image2 image3 k ROI ROIs ROIsNb Start tempROIsNb name imagename tempidx Raster

for i=1:numel(GoodClusterData)
corr_temp=zeros(size(GoodClusterData(i).ZS,1),1);
parfor j=1:size(GoodClusterData(i).ZS,1)
temp=corrcoef(GoodClusterData(i).mean, GoodClusterData(i).ZS(j,:));
corr_temp(j)=temp(1,2);
end
GoodClusterData(i).CorrCoef=corr_temp;
end
GoodClusters_goodmembers=[];
for i=1:length(GoodBetas_select)
%GoodClusters_goodmembers(i).Spikes=GoodClusterData(i).Spikes(find(GoodClusterData(i).CorrCoef>=0.5),:);
GoodClusters_goodmembers(i).ZS=GoodClusterData(i).ZS(find(GoodClusterData(i).CorrCoef>=0.5),:);
temp=find(idxKmeans_ZS==GoodBetas_select(i));
GoodClusters_goodmembers(i).idx=temp(find(GoodClusterData(i).CorrCoef>=0.5));
GoodClusters_goodmembers(i).mean=mean(GoodClusters_goodmembers(i).ZS,1);
GoodClusters_goodmembers(i).STD=std(GoodClusters_goodmembers(i).ZS,1,1);
idx=find(idxKmeans_ZS==GoodBetas_select(i));
idx=idx(find(GoodClusterData(i).CorrCoef>=0.5));
%GoodClusters_goodmembers(i).Fish=idx_Fish(idx);
end

figure;counter=1;
for i=1:numel(GoodClusterData)
    for state=1:4
       idx=find(GoodClusterData(i).State==state);
       subplot(numel(GoodClusterData),4,counter);plot(mean(GoodClusterData(i).ZS(idx,:),1));xlim([0 size(Cmap_ZS,2)]);
       counter=counter+1;
    end
end