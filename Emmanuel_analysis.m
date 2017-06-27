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
Noise=vertcat(Noise,N);
GN=N(F,:);
Calcium=vertcat(Calcium,C);
%DF=vertcat(DF,D);
Spikes=vertcat(Spikes,S);
Fitness=horzcat(Fitness,F);
GoodCalcium=vertcat(GoodCalcium,GC);
GoodNoise=vertcat(GoodNoise,GN);
%GoodDF=vertcat(GoodDF,GD);
GoodSpikes=vertcat(GoodSpikes,GS);
MatFiles(i).number=size(Calcium,1);
MatFiles(i).GoodNumber=MatFiles(i-1).GoodNumber+length(F);
end
clearvars GC C S F N name i GS GN N;
ZS=zscore(GoodCalcium,1,2);
temp=GoodCalcium;
%temp(:,5591:5600)=repmat(GoodCalcium(:,5590),1,10);
temp(:,9291:9300)=repmat(GoodCalcium(:,9290),1,10);
temp(:,1:30)=repmat(GoodCalcium(:,31),1,30);
ZS=zscore(temp,1,2);

Merged_DF=ZS;%   ALL_DF_02=Dataset (neurons x time)
MergedIDX=cell(1,length(Merged_DF)); % Where the IDs of merged ROIs will be stored

Merged_DFb=Merged_DF;
MergedIDXb=MergedIDX;

counter=1;threshold=0.8;
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
            Merged_DF(counter,:)=nanmean(ZS(MergedIDX{counter},:),1); %average from dataset, easier than adding one at a time with a factor                
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
temp=find(numMerged>50);

Model_merged=[];counter=1;
for i=temp'
    %mdl=stepwiselm(Regressor',DataSet(i,:),'linear','Criterion','adjrsquared','Intercept',false,'Upper','interactions','Verbose',0);
    mdl=stepwiselm(Stimuli',Merged_DF(i,:),'linear','Criterion','adjrsquared','Upper','interactions','Verbose',0);
    Model_merged(counter).coef=mdl.Coefficients;
    Model_merged(counter).MSE=mdl.MSE;
    Model_merged(counter).Fitted=mdl.Fitted;
    Model_merged(counter).rsquared=mdl.Rsquared.Adjusted;
    mdl=stepwiselm(Stimuli_off',Merged_DF(i,:),'linear','Criterion','adjrsquared','Upper','interactions','Verbose',0);
    Model_merged(counter).coef_off=mdl.Coefficients;
    Model_merged(counter).MSE_off=mdl.MSE;
    Model_merged(counter).Fitted_off=mdl.Fitted;
    Model_merged(counter).rsquared_off=mdl.Rsquared.Adjusted;
    counter=counter+1;
end
SelectedSamples=find([Model_merged.rsquared]>0.01);
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1500, 1000]);
counter=1;xplot=floor(sqrt(length(SelectedSamples)));yplot=ceil(length(SelectedSamples)/xplot);
for i=SelectedSamples    
    NumberOfCells=numel(MergedIDX{temp(i)});
    subplot(xplot,yplot,counter);plot(x,Merged_DF(temp(i),:));title(num2str(NumberOfCells))
    xlim([0 size(Merged_DF,2)])
    counter=counter+1;
end

SelectedSamples_off=find([Model_merged.rsquared_off]>0.01);
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1500, 1000]);
counter=1;xplot=floor(sqrt(length(SelectedSamples_off)));yplot=ceil(length(SelectedSamples_off)/xplot);
for i=SelectedSamples_off 
    NumberOfCells=numel(MergedIDX{temp(i)});
    subplot(xplot,yplot,counter);plot(x,Merged_DF(temp(i),:));title(num2str(NumberOfCells))
    xlim([0 size(Merged_DF,2)])
    counter=counter+1;
end

x = linspace(1,size(GoodCalcium,2),size(GoodCalcium,2));
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;xplot=floor(sqrt(length(temp)));yplot=ceil(length(temp)/xplot);
for i=temp'
    NumberOfCells=numel(MergedIDX{i});
    %subplot(5,1,counter);plot(x,Cmap_ZS(i,:),x,Model_DF(i).Fitted);title(num2str(NumberOfCells))
    %subplot(xplot,yplot,counter);plot(Merged_DF(i,:));title(num2str(NumberOfCells))
    subplot(length(temp),2,counter);plot(Merged_DF(i,:));title(num2str(NumberOfCells))
    subplot(length(temp),2,counter+1);imagesc(ZS(MergedIDX{i},:));colormap hot;title(num2str(NumberOfCells))
    xlim([0 size(GoodCalcium,2)])
    counter=counter+2;
end



Stimuli=zeros(3,size(GoodCalcium,2));
GCaMP6=[-0.104392135015146,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
GCaMP6=interp(GCaMP6,2);
Stimulus_train=[1 2 3 3 2 1 1 3 2 2 1 3 1 3 2 2 1 3 3 1 2 3 2 1 1 2 3 2 3 1]; %1 = loom, 2=dim,3=check
idx=115;
for i=Stimulus_train
    switch i
        case 1
            Stimuli(i,idx:idx+size(GCaMP6,1)-1)=GCaMP6';
        case 2
            Stimuli(i,idx:idx+size(GCaMP6,1)-1)=GCaMP6';
        case 3
            Stimuli(i,idx:idx+size(GCaMP6,1)-1)=GCaMP6';
    end    
    idx=idx+186;
end

Stimuli_off=zeros(3,size(GoodCalcium,2));
GCaMP6=[-0.104392135015146,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
GCaMP6=interp(GCaMP6,2);
Stimulus_train=[1 2 3 3 2 1 1 3 2 2 1 3 1 3 2 2 1 3 3 1 2 3 2 1 1 2 3 2 3 1]; %1 = loom, 2=dim,3=check
idx=187;
for i=Stimulus_train
    switch i
        case 1
            Stimuli_off(i,idx:idx+size(GCaMP6,1)-1)=GCaMP6';
        case 2
            Stimuli_off(i,idx:idx+size(GCaMP6,1)-1)=GCaMP6';
        case 3
            Stimuli_off(i,idx:idx+size(GCaMP6,1)-1)=GCaMP6';
    end    
    idx=idx+186;
end
Stimuli_off=Stimuli_off(1:3,1:size(GoodCalcium,2));

options = statset('UseParallel',1); [idxKmeans_ZS Cmap_ZS]=kmeans(ZS,50,'Options',options,'Distance','correlation','Replicates',5,'MaxIter',1000,'Display','final');
[Model_ZS,GoodBetas_ZS]=Test_Regress(Cmap_ZS,Stimuli,idxKmeans_ZS,0.1);
[Model_ZS_off,GoodBetas_ZS_off]=Test_Regress(Cmap_ZS,Stimuli_off,idxKmeans_ZS,0.1);


ModelResults=[];
parfor i=1:length(temp)
    %mdl=stepwiselm(Regressor',DataSet(i,:),'linear','Criterion','adjrsquared','Intercept',false,'Upper','interactions','Verbose',0);
    mdl=stepwiselm(Stimuli',Merged_DF(temp(i),:),'linear','Criterion','adjrsquared','Upper','interactions','Verbose',0);
    ModelResults(i).coef=mdl.Coefficients;
    ModelResults(i).MSE=mdl.MSE;
    ModelResults(i).Fitted=mdl.Fitted;
    ModelResults(i).rsquared=mdl.Rsquared.Adjusted;
end


figure;
counter=1;
for i=[7 8 9 12]
    subplot(4,1,counter);plot(Merged_DF(temp(i),:)*5);hold on;plot(Stimuli(1,:),'color',[1 0 0]);hold on;plot(Stimuli(2,:),'color',[0 1 0]);hold on;plot(Stimuli(3,:),'color',[0 0 1]);
    counter=counter+1;
end

ModelResults=[];
parfor i=1:size(ZS,1)
    %mdl=stepwiselm(Regressor',DataSet(i,:),'linear','Criterion','adjrsquared','Intercept',false,'Upper','interactions','Verbose',0);
    mdl=stepwiselm(Stimuli',ZS(i,:),'linear','Criterion','adjrsquared','Upper','interactions','Verbose',0);
    ModelResults(i).coef=mdl.Coefficients;        
    ModelResults(i).rsquared=mdl.Rsquared.Adjusted;
    mdl=stepwiselm(Stimuli_off',ZS(i,:),'linear','Criterion','adjrsquared','Upper','interactions','Verbose',0);
    ModelResults(i).coef_off=mdl.Coefficients;        
    ModelResults(i).rsquared_off=mdl.Rsquared.Adjusted;
end

Stimuli_long=zeros(2,size(GoodCalcium,2));
GCaMP6=[-0.104392135015146,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
GCaMP6=interp(GCaMP6,2);
%Stimulus_train=[1 2 3 3 2 1 1 3 2 2 1 3 1 3 2 2 1 3 3 1 2 3 2 1 1 2 3 2 3 1]; %1 = loom, 2=dim,3=check
idx=120;
for i=1:50
    Stimuli_long(1,idx:idx+size(GCaMP6,1)-1)=GCaMP6';    
    Stimuli_long(2,idx+60:idx+60+size(GCaMP6,1)-1)=GCaMP6';
    idx=idx+186;
end
Stimuli_long=Stimuli_long(:,1:size(GoodCalcium,2));
[Model_ZS_long2,GoodBetas_ZS_long2]=Test_Regress(Cmap_ZS2,Stimuli_long,idxKmeans_ZS2,0.05);

Numbers=[0 [MatFiles.GoodNumber]];
counter=1;
idx_Plane=nan(length(GoodCalcium),1);
idx_age=nan(length(GoodCalcium),1);
name=strcat(MatFiles(1).name);
[Plane,~]=regexp(name,'_depth(\d)','tokens','match');Plane=str2num(Plane{1}{1});
[Age,~]=regexp(name,'_age(\d+)','tokens','match');Age=str2num(Age{1}{1});
idx_Plane(1:Numbers(2))=Plane;
idx_age(1:Numbers(2))=Age;
for i=2:length(MatFiles)    
    name=strcat(MatFiles(i).name);    
    [Plane,~]=regexp(name,'_depth(\d)','tokens','match');Plane=str2num(Plane{1}{1});
    [Age,~]=regexp(name,'_age(\d+)','tokens','match');Age=str2num(Age{1}{1});
    idx_Plane(Numbers(i):Numbers(i+1))=Plane;
    idx_age(Numbers(i):Numbers(i+1))=Age;
end

GoodBetas_select=GoodBetas_ZS_off;
GoodBetas_select=GoodBetas_ZS_long([3 5 7 10 13 14]);

counter=1;counter2=1;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
GoodClusterData=[];
%name=['ELO';'CLO'];
for i=GoodBetas_select  
    idx=find(idxKmeans_ZS==i);
    GoodClusterData(counter).mean=mean(ZS(idx,:),1);
    GoodClusterData(counter).ZS=ZS(idx,:);
    GoodClusterData(counter).planes=idx_Plane(idx);
    GoodClusterData(counter).age=idx_age(idx);
    subplot(length(GoodBetas_select),4,counter2);plot(GoodClusterData(counter).mean(1:5550));title(num2str(length(idx)))    
    subplot(length(GoodBetas_select),4,counter2+1);imagesc(GoodClusterData(counter).ZS, [0 4]); colormap hot;title(num2str(length(idx)))    
    subplot(length(GoodBetas_select),4,counter2+2);histogram(GoodClusterData(counter).planes);
    subplot(length(GoodBetas_select),4,counter2+3);histogram(GoodClusterData(counter).age);
    %subplot(length(GoodBetas_select),4,counter2+3);bar([sum(GoodClusterData(counter).State==3) sum(GoodClusterData(counter).State==4)]);set(gca,'xticklabel',name);
    counter2=counter2+4;
    counter=counter+1;
end

counter=1;counter2=1;
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
for idx=1:length(GoodClusterData)
    idx_age_temp=find([GoodClusterData(idx).age]==6);
    subplot(length(GoodBetas_select),2,counter2);plot(mean(GoodClusterData(idx).ZS(idx_age_temp,1:5550)));title(strcat('Age 6, nb cells=',num2str(length(idx_age_temp))))
    idx_age_temp=find(GoodClusterData(idx).age==12);
    subplot(length(GoodBetas_select),2,counter2+1);plot(mean(GoodClusterData(idx).ZS(idx_age_temp,1:5550)));title(strcat('Age 12, nb cells=',num2str(length(idx_age_temp))))
    counter2=counter2+2;
    counter=counter+1;
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
for i=GoodBetas_select
    %temp{counter}=find(idxKmeans_ZS==i);
    temp{counter}=find(idxKmeans_ZS_goodmembers==i);
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
            image2=image2+ROI;image2=(image2/max(max(image2)))*68;image2=uint8(image2);
            for j=1:3
                image3(:,:,j)=image3(:,:,j)+image2*colors(ClusterNb(k),j);
            end
        end
        %image3(:,:,3)=image;
    end
    %image3=uint8(image3);
    name=strcat('Kmeans_good_',imagename(4:end));
    imwrite(image3,name,'tif');
end

x = linspace(1,size(Cmap,2),size(Cmap,2));
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1300, 900]);
counter=1;xplot=floor(sqrt(length(GoodBetas_select)));yplot=ceil(length(GoodBetas_select)/xplot);
for i=GoodBetas_select    
    NumberOfCells=length(find(idxKmeans==i));
    %NumberOfCells=length(find(idxKmeans_ZS_goodmembers==i));
    subplot(5,1,counter);plot(x,Cmap(i,:),'color',colors(counter,:));title(num2str(NumberOfCells))
    %subplot(xplot,yplot,counter);plot(GoodClusters_goodmembers(counter).mean(10:1000),'color',colors(counter,:));title(num2str(NumberOfCells))
%    xlim([0 size(Cmap_ZS,2)])
    counter=counter+1;
end
clearvars idx i temp tempFileNb fileNb AVG_files filename image counter Numbers image2 image3 k ROI ROIs ROIsNb Start tempROIsNb name imagename tempidx Raster


[Model_ZS_old,GoodBetas_ZS_old]=Test_Regress(Cmap_ZS_old,Stimuli,idxKmeans_ZS_old,0.1);
[Model_ZS_off_old,GoodBetas_ZS_off_old]=Test_Regress(Cmap_ZS_old,Stimuli_off,idxKmeans_ZS_old,0.1);

opol = 6;idx=1;x=x(1:9299);
[p,s,mu] = polyfit(x,Cmap_ZS(idx,1:9299),opol);
f_y = polyval(p,x,[],mu);
dt_ecgnl = Cmap_ZS(idx,1:9299) - f_y;
subplot(2,1,1)
plot(x,Cmap_ZS(idx,1:9299)), grid
title 'Detrended ECG Signals', ylabel 'Voltage (mV)'

opol = 6;idx=1;x=x(1:9299);
[p,s,mu] = polyfit(x,Cmap_ZS(idx,1:9299),opol);
f_y = polyval(p,x,[],mu);
dt_ecgnl = Cmap_ZS(idx,1:9299) - f_y;
subplot(2,1,1)
plot(x,Cmap_ZS(idx,1:9299)), grid
title 'Detrended ECG Signals', ylabel 'Voltage (mV)'

subplot(2,1,2)
plot(x,dt_ecgnl), grid
xlabel Sample, ylabel 'Voltage (mV)'

temp2=temp;
x = linspace(1,size(GoodCalcium,2),size(GoodCalcium,2));
parfor idx=1:size(temp,1)
p = polyfit(x, temp(idx,:), 4);
temp2(idx,:) = temp(idx,:) - polyval(p,x);
end
ZS=zscore(temp2,1,2);