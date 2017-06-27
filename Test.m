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
clearvars GC C S F N name i GS;

ZS=zscore(GoodCalcium,1,2);
options = statset('UseParallel',1); [idxKmeans_ZS Cmap_ZS]=kmeans(ZS,50,'Options',options,'Distance','correlation','Replicates',5,'MaxIter',1000,'Display','final');

[Model_GM,GoodBetas_GM]=Test_Regress(Cmap_GM,Stimuli,idxKmeans_GM,0.1);
[Model_GM_off,GoodBetas_GM_off]=Test_Regress(Cmap_GM,Stimuli_off,idxKmeans_GM,0.1);

temp=GoodCalcium;
temp(:,5591:5600)=repmat(GoodCalcium(:,5590),1,10);
temp2=temp;
x = linspace(1,size(GoodCalcium,2),size(GoodCalcium,2));
parfor idx=1:size(temp,1)
    p = polyfit(x, temp(idx,:), 4);
    temp2(idx,:) = temp(idx,:) - polyval(p,x);
end
ZS=zscore(temp2,1,2);

ZS_std=std(ZS,1,2);
idx_std=find(ZS_std<0.4);
ZS(idx_std,:)=NaN;

options = statset('UseParallel',1); [idxKmeans_ZS Cmap_ZS]=kmeans(ZS,50,'Options',options,'Distance','correlation','Replicates',5,'MaxIter',1000,'Display','final');
[Model_ZS,GoodBetas_ZS]=Test_Regress(Cmap_ZS,Stimuli,idxKmeans_ZS,0.1);
[Model_ZS_off,GoodBetas_ZS_off]=Test_Regress(Cmap_ZS,Stimuli_off,idxKmeans_ZS,0.1);

idx_on=find([ModelResults_check.rsquared]>0.1);
figure;plot(mean(ZS(idx_on,:),1));

x = linspace(1,size(Cmap_ZS,2),size(Cmap_ZS,2));
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;xplot=floor(sqrt(length(idx_on)));yplot=ceil(length(idx_on)/xplot);
for i=idx_on
    subplot(xplot,yplot,counter);plot(ZS(i,:));
    counter=counter+1;
end

idx_off=find([ModelResults_check.rsquared_off]>0.1);
figure;plot(mean(ZS(idx_off,:),1));

counter=1;counter2=1;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
OffCells=[];
idx=idx_off;
OffCells.mean=mean(ZS(idx,:),1);
OffCells.ZS=ZS(idx,:);
OffCells.planes=idx_Plane(idx);
OffCells.age=idx_age(idx);
subplot(2,2,counter2);plot(OffCells.mean(1:5550));title(num2str(length(idx)))
subplot(2,2,counter2+1);imagesc(OffCells.ZS, [0 4]); colormap hot;
subplot(2,2,counter2+2);histogram(OffCells.planes);
subplot(2,2,counter2+3);histogram(OffCells.age);
%subplot(length(GoodBetas_select),4,counter2+3);bar([sum(GoodClusterData(counter).State==3) sum(GoodClusterData(counter).State==4)]);set(gca,'xticklabel',name);

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
ONCells=[];
idx=idx_on;
ONCells.mean=mean(ZS(idx,:),1);
ONCells.ZS=ZS(idx,:);
ONCells.planes=idx_Plane(idx);
ONCells.age=idx_age(idx);
subplot(2,2,1);plot(ONCells.mean(1:5550));title(num2str(length(idx)))
subplot(2,2,2);imagesc(ONCells.ZS, [0 4]); colormap hot;
subplot(2,2,3);histogram(ONCells.planes);
subplot(2,2,4);histogram(ONCells.age);

idx_off_crap=find([ModelResults_check.rsquared_off]>0.05 & [ModelResults_check.rsquared_off]<0.1);
counter=1;counter2=1;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
OffCrapCells=[];
idx=idx_off_crap;
OffCrapCells.mean=mean(ZS(idx,:),1);
OffCrapCells.ZS=ZS(idx,:);
OffCrapCells.planes=idx_Plane(idx);
OffCrapCells.age=idx_age(idx);
subplot(2,2,counter2);plot(OffCrapCells.mean(1:5550));title(num2str(length(idx)))
subplot(2,2,counter2+1);imagesc(OffCrapCells.ZS, [0 4]); colormap hot;
subplot(2,2,counter2+2);histogram(OffCrapCells.planes);
subplot(2,2,counter2+3);histogram(OffCrapCells.age);
%subplot(length(GoodBetas_select),4,counter2+3);bar([sum(GoodClusterData(counter).State==3) sum(GoodClusterData(counter).State==4)]);set(gca,'xticklabel',name);

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);counter=1;
Stimulus_train=[1 2 3 3 2 1 1 3 2 2 1 3 1 3 2 2 1 3 3 1 2 3 2 1 1 2 3 2 3 1]; %1 = loom, 2=dim,3=check

for age=[6 9 12]
    idx=find([OffCells.age]==age);
    subplot(3,1,counter);plot(mean(OffCells.ZS(idx,:),1));ylim([-1, 3]);title(num2str(age));
    start=120;
    for i=Stimulus_train
        switch i
            case 1
                rectangle('FaceColor','r','Position',[start,-1,80,0.2]);
            case 2
                rectangle('FaceColor','b','Position',[start,-1,80,0.2]);
            case 3
                rectangle('FaceColor','g','Position',[start,-1,80,0.2]);
        end
        start=start+186;
    end
    counter=counter+1;
end

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);counter=1;
for age=[6 9 12]
    idx=find([OffCrapCells.age]==age);
    subplot(3,1,counter);plot(mean(OffCrapCells.ZS(idx,:),1));ylim([-1, 3]);title(num2str(age));
    start=120;
    for i=Stimulus_train
        switch i
            case 1
                rectangle('FaceColor','r','Position',[start,-1,80,0.2]);
            case 2
                rectangle('FaceColor','b','Position',[start,-1,80,0.2]);
            case 3
                rectangle('FaceColor','g','Position',[start,-1,80,0.2]);
        end
        start=start+186;
    end
    counter=counter+1;
end

figure;imagesc(OffCells.ZS(287:310,:), [0 4]); colormap hot;


CheckerBoardTiming=[];
Stimulus_train=[1 2 3 3 2 1 1 3 2 2 1 3 1 3 2 2 1 3 3 1 2 3 2 1 1 2 3 2 3 1]; %1 = loom, 2=dim,3=check
idx=187;
for i=Stimulus_train
    switch i
        case 1
            idx=idx+186;
        case 2            
            idx=idx+186;
        case 3
            CheckerBoardTiming=[CheckerBoardTiming idx];
            idx=idx+186;
    end    
    
end

SlicedZS=[];
for i=CheckerBoardTiming
    SlicedZS=[SlicedZS ZS(:,i-10:i+90)];
end

Stimuli_off(i,idx:idx+size(GCaMP6,1)-1)=GCaMP6';

GCaMP6=[-0.104392135015146,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
GCaMP6=interp(GCaMP6,2);
CheckerBoard_stimuli=zeros(1,size(SlicedZS,2));
for i=1:10
    CheckerBoard_stimuli(1,10+((i-1)*100):10+((i-1)*100)+size(GCaMP6,1)-1)=GCaMP6';
end


Coeff_check=[];
parfor i=1:size(SlicedZS,1)
    temp=corrcoef(SlicedZS(i,:), CheckerBoard_stimuli);
    Coeff_check(i)=temp(1,2);
end
idx=find(Coeff_check>0.3);
CheckCells=[];
CheckCells.mean=mean(ZS(idx,:),1);
CheckCells.ZS=ZS(idx,:);
CheckCells.planes=idx_Plane(idx);
CheckCells.age=idx_age(idx);

counter=1;counter2=1;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
subplot(2,2,counter2);plot(mean(ZS(idx,:),1));title(num2str(length(idx)))
subplot(2,2,counter2+1);imagesc(ZS(idx,:), [0 4]); colormap hot;
subplot(2,2,counter2+2);histogram(idx_Plane(idx));
subplot(2,2,counter2+3);histogram(idx_age(idx));
%subplot(length(GoodBetas_select),4,counter2+3);bar([sum(GoodClusterData(counter).State==3) sum(GoodClusterData(counter).State==4)]);set(gca,'xticklabel',name);

Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);counter=1;

for age=[6 9 12]
    idx=find([CheckCells.age]==age);
    subplot(3,1,counter);plot(mean(CheckCells.ZS(idx,:),1));ylim([-1, 3]);title(num2str(age));
    start=120;
    for i=Stimulus_train
        switch i
            case 1
                rectangle('FaceColor','r','Position',[start,-1,80,0.2]);
            case 2
                rectangle('FaceColor','b','Position',[start,-1,80,0.2]);
            case 3
                rectangle('FaceColor','g','Position',[start,-1,80,0.2]);
        end
        start=start+186;
    end
    counter=counter+1;
end

GoodBetas_select=[1 2];
idxKmeans_ZS_off=zeros(size(ZS,1));
idxKmeans_ZS_off(idx_off)=1;
idxKmeans_ZS_off(idx_off_crap)=2;
counter=1;
for i=GoodBetas_select
    temp{counter}=find(idxKmeans_ZS_off==i);
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
    NumberOfCells=length(find(idxKmeans_ZS_off==i));
    %subplot(5,1,counter);plot(x,Cmap_ZS(i,:),x,Model_DF(i).Fitted);title(num2str(NumberOfCells))
    subplot(xplot,yplot,counter);plot(mean(ZS(find(idxKmeans_ZS_off==i),:),1),'color',colors(counter,:));title(num2str(NumberOfCells))
    xlim([0 size(Cmap_ZS,2)])
    counter=counter+1;
end
clearvars idx i temp tempFileNb fileNb AVG_files filename image counter Numbers image2 image3 k ROI ROIs ROIsNb Start tempROIsNb name imagename tempidx Raster
