filelist=dir('*.tif');
for File=1:length(filelist)
    [mixedsig, mixedfilters, CovEvals, covtrace, movm, movtm] = CellsortPCA(filelist(File).name,[],100,1,'C:\Temp\CROP');
    [ica_sig, ica_filters, ica_A, numiter] = CellsortICA(mixedsig, mixedfilters, CovEvals, [1:size(mixedsig,1)],  0.8,size(mixedsig,1),randn(size(mixedsig,1), size(mixedsig,1)),1e-6,50000);
    PCA_ICA_results(File).ica_sig=ica_sig;
    PCA_ICA_results(File).ica_filters=ica_filters;
end
%save('PCA_ICA_results_Lucy_SynGCaMP_prelim.mat');

% idx=1;AllTraces=PCA_ICA_results(idx).ica_sig;
% for idx=2:length(PCA_ICA_results)    
%     AllTraces=vertcat(AllTraces,PCA_ICA_results(idx).ica_sig);
% end
%     
% ZS=zscore(AllTraces,1,2);
% options = statset('UseParallel',1); [idxKmeans_ZS Cmap_ZS]=kmeans(ZS,20,'Options',options,'Distance','correlation','Replicates',5,'MaxIter',1000,'Display','final');
% [Model_ZS,GoodBetas_ZS]=Test_Regress(Cmap_ZS,Short_Stim,idxKmeans_ZS,0.1);
% 
% Fighandle=figure;
% set(Fighandle, 'Position', [100, 100, 1400, 900]);
% imagesc(ZS,[0 4]);colormap hot
% 
% GoodBetas_select=GoodBetas_ZS;
% GoodClustersData=[];
% for i=1:length(GoodBetas_select)
% GoodClustersData(i).DF=AllTraces(idxKmeans_ZS==GoodBetas_select(i),:)*100;
% GoodClustersData(i).Mean=mean(GoodClustersData(i).DF,1);
% GoodClustersData(i).STD=std(GoodClustersData(i).DF,1,1);
% end
% 
% for i=1:numel(GoodClustersData)
% corr_temp=zeros(size(GoodClustersData(i).DF,1),1);
% parfor j=1:size(GoodClustersData(i).DF,1)
% temp=corrcoef(GoodClustersData(i).Mean, GoodClustersData(i).DF(j,:));
% corr_temp(j)=temp(1,2);
% end
% GoodClustersData(i).CorrCoef=corr_temp;
% end
% 
% GoodClusters_goodmembers=[];
% for i=1:length(GoodBetas_select)
% %GoodClusters_goodmembers(i).Spikes=GoodClustersData(i).Spikes(find(GoodClustersData(i).CorrCoef>=0.5),:);
% GoodClusters_goodmembers(i).ZS=zscore(GoodClustersData(i).DF(find(GoodClustersData(i).CorrCoef>=0.5),:),1,2);
% temp=find(idxKmeans_ZS==GoodBetas_select(i));
% GoodClusters_goodmembers(i).idx=temp(find(GoodClustersData(i).CorrCoef>=0.5));
% GoodClusters_goodmembers(i).mean=mean(GoodClusters_goodmembers(i).ZS,1);
% GoodClusters_goodmembers(i).STD=std(GoodClusters_goodmembers(i).ZS,1,1);
% idx=find(idxKmeans_ZS==GoodBetas_select(i));
% idx=idx(find(GoodClustersData(i).CorrCoef>=0.5));
% %GoodClusters_goodmembers(i).Fish=idx_Fish(idx);
% end
% 
% Fighandle=figure;
% set(Fighandle, 'Position', [100, 100, 1400, 900]);x = linspace(0.2,size(Cmap_ZS,2)/5,size(Cmap_ZS,2));
% counter=1;counter2=1;xplot=floor(sqrt(length(GoodBetas_select)));yplot=ceil(length(GoodBetas_select)/xplot);
% for i=1:length(GoodBetas_select)  
%     if counter==3
%         counter=counter+1;
%     end
%     subplot(3,3,counter);plot(x,GoodClusters_goodmembers(i).mean,'color',colors(counter2,:)/256);axis([0 131 -1 4]);rectangle('FaceColor','r','Position',[11 -1 10 0.25]);rectangle('FaceColor','r','Position',[51 -1 10 0.25]);rectangle('FaceColor','r','Position',[91 -1 10 0.25]);rectangle('FaceColor','b','Position',[31 -1 10 0.25]);rectangle('FaceColor','b','Position',[71 -1 10 0.25]);rectangle('FaceColor','b','Position',[111 -1 10 0.25]);
%     counter=counter+1;
%     counter2=counter2+1;
% end
% 
% Long_Stim=zeros(4,size(AllTraces,2));
% check=[45,105,465,1005,1184];
% dim=[167,348,407,587,787];
% loom=[225,525,825,1065,1125];
% circle=[285,645,705,775,945];
% GCaMP6=[-0.104392135015146,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
% GCaMP6=interp(GCaMP6,2);
% for i=1:length(check)
%     Long_Stim(1,check(i):check(i)+size(GCaMP6,1)-1)=GCaMP6';
%     Long_Stim(1,check(i)+1500:check(i)+1500+size(GCaMP6,1)-1)=GCaMP6';
%     Long_Stim(2,dim(i):dim(i)+size(GCaMP6,1)-1)=GCaMP6';
%     Long_Stim(2,dim(i)+1500:dim(i)+1500+size(GCaMP6,1)-1)=GCaMP6';
%     Long_Stim(3,loom(i):loom(i)+size(GCaMP6,1)-1)=GCaMP6';
%     Long_Stim(3,loom(i)+1500:loom(i)+1500+size(GCaMP6,1)-1)=GCaMP6';
%     Long_Stim(4,circle(i):circle(i)+size(GCaMP6,1)-1)=GCaMP6';
%     Long_Stim(4,circle(i)+1500:circle(i)+1500+size(GCaMP6,1)-1)=GCaMP6';
% end
% clearvars GCaMP6 back back_off fwd fwd_off check loom dim circle idx i k
% 
% [Model_ZS,GoodBetas_ZS]=Test_Regress(Cmap_ZS,Long_Stim,idxKmeans_ZS,0.05);


Short_Stim=zeros(2,size(AllTraces,2));
GCaMP6=[-0.104392135015146,1.69644104899772,5.13796058542217,8.27886020152244,10.3756715204800,11.8173714529814,12.2425184714093,10.8571417354877,8.80831829681196,6.91339112244670,5.46959264663869,4.30868766622567,3.42533619066766,2.75378443486879,2.18017250852183,1.72816235135824,1.32732537295463,1.00684435500268,0.730210038304555,0.530242444093118,0.362253250339685,0.227668255288566,0.0869242416152502,0.000718266708050853,-0.0828334873368325]';
GCaMP6=interp(GCaMP6,2);
stim_present1 = [96,296,496,696,896,1096,1296,1496,1696,1896];
stim_present2 = [56,256,456,656,856,1056,1256,1456,1656,1856];    
for i=1:length(stim_present1)
    Short_Stim(1,stim_present1(i):stim_present1(i)+size(GCaMP6,1)-1)=GCaMP6';
    Short_Stim(2,stim_present2(i):stim_present2(i)+size(GCaMP6,1)-1)=GCaMP6';
end

% [Model_ZS,GoodBetas_ZS]=Test_Regress(Cmap_ZS,Short_Stim,idxKmeans_ZS,0.2);
% [Model_GM,GoodBetas_GM]=Test_Regress(Cmap_GM,Short_Stim,idxKmeans_GM,0.2);
% 

% ModelResults=[];
% parfor i=1:length(ZS)
%     %mdl=stepwiselm(Regressor',DataSet(i,:),'linear','Criterion','adjrsquared','Intercept',false,'Upper','interactions','Verbose',0);
%     mdl=stepwiselm(Short_Stim',ZS(i,:),'linear','Criterion','adjrsquared','Upper','linear','Verbose',0);
%     ModelResults(i).coef=mdl.Coefficients;
%     ModelResults(i).MSE=mdl.MSE;
%     ModelResults(i).Fitted=mdl.Fitted;
%     ModelResults(i).rsquared=mdl.Rsquared.Adjusted;
% end
% idx_rsq=find([ModelResults.rsquared]>0.1);
% Fighandle=figure;
% set(Fighandle, 'Position', [100, 100, 1400, 900]);
% 
% imagesc(ZS(idx_rsq,:),[0 4]);colormap hot
% 
% ModelResults=[];
% parfor i=1:size(ZS,1)
%     %mdl=stepwiselm(Regressor',DataSet(i,:),'linear','Criterion','adjrsquared','Intercept',false,'Upper','interactions','Verbose',0);
%     mdl=stepwiselm(Long_Stim',ZS(i,:),'linear','Criterion','adjrsquared','Upper','linear','Verbose',0);
%     ModelResults(i).coef=mdl.Coefficients;
%     ModelResults(i).MSE=mdl.MSE;
%     ModelResults(i).Fitted=mdl.Fitted;
%     ModelResults(i).rsquared=mdl.Rsquared.Adjusted;
% end
% idx_rsq=find([ModelResults.rsquared]>0.2);
% GoodModels=ModelResults(idx_rsq);
% 
% fps=10;
% x = linspace(1,size(ZS,2),size(ZS,2));
% Fighandle=figure;
% set(Fighandle, 'Position', [100, 100, 1300, 900]);
% counter=1;xplot=floor(sqrt(length(idx_rsq)));yplot=ceil(length(idx_rsq)/xplot);
% for i=idx_rsq
%     subplot(xplot,yplot,counter);plot(x,ZS(i,:),x,GoodModels(counter).Fitted);legend('ZS trace','Fitted LR','Location','northeast')    
%     xlim([0 size(ZS,2)])
%     counter=counter+1;
% end
% 
% coefficients_above02=[];
% for idx=1:length(idx_rsq)
%     coef=[GoodModels(idx).coef];
%     for coef_idx=2:height(coef)
%         if coef.pValue(coef_idx)<0.05
%             coefficients_above02{idx,coef_idx-1}=coef.tStat(coef_idx);
%         end
%     end    
% end
% idxempty=cellfun('isempty',coefficients_above02);
% coefficients_above02(idxempty)={0};
% clearvars idxempty idx coef_idx coef
% coefficients_above02=cell2mat(coefficients_above02);
% [max_coef,idx_max]=max(coefficients_above02,[],2);

% nb_components=100;counter=1;
% for i=idx_rsq
%     idx=floor(i/nb_components)+1;
%     name=filelist(idx).name;
%     imagename=strcat('AVG_',name);
%     image=double(imread(imagename));image=image/max(prctile(image,95));image=image*64;
%     image=uint8(image);image3=repmat(image,1,1,3);
%     ROI=squeeze(PCA_ICA_results(idx).ica_filters(mod(i,nb_components),:,:));ROI(ROI<max(prctile(ROI,40)))=0;ROI=ROI/max(max(ROI));ROI=ROI*250;
%     image3(:,:,1)=image+uint8(ROI);
%     filename=strcat('GoodTraces_',num2str(counter),'_',name);
%     imwrite(image3,filename,'tif');
%     counter=counter+1;
% end


% fps=10;
% x = linspace(1,size(ZS,2),size(ZS,2));
% Fighandle=figure;
% set(Fighandle, 'Position', [100, 100, 1300, 900]);
% counter=1;xplot=floor(sqrt(length(idx_rsq)));yplot=ceil(length(idx_rsq)/xplot);
% for i=idx_rsq
%     subplot(xplot,yplot,counter);plot(x,ZS(i,:)+15,x,Long_Stim(4,:));legend('ZS trace','Circles','Location','northeast')    
%     xlim([0 size(ZS,2)])
%     counter=counter+1;
% end
% 
% [ica_segments, segmentlabel, segcentroid] = CellsortSegmentation(ica_filters, smwidth, thresh, arealims, plotting)

for File=1:length(filelist)    
    [PCA_ICA_results(File).ica_segments, PCA_ICA_results(File).segmentlabel, PCA_ICA_results(File).segcentroid] = CellsortSegmentation(PCA_ICA_results(File).ica_filters, 2, 2, 20, 0)
    PCA_ICA_results(File).cell_sig = CellsortApplyFilter(filelist(File).name, PCA_ICA_results(File).ica_segments)
end

AllSegTraces=[];
idx=1;AllSegTraces=PCA_ICA_results(idx).cell_sig;
for idx=2:length(PCA_ICA_results)    
    AllSegTraces=vertcat(AllSegTraces,PCA_ICA_results(idx).cell_sig);
end

DF=DeltaF2(AllSegTraces,11,21);
ZS=zscore(AllSegTraces,1,2);
% ModelResultsSeg=[];
% parfor i=1:length(DF)
%     %mdl=stepwiselm(Regressor',DataSet(i,:),'linear','Criterion','adjrsquared','Intercept',false,'Upper','interactions','Verbose',0);
%     mdl=stepwiselm(Short_Stim',DF(i,:),'linear','Criterion','adjrsquared','Upper','linear','Verbose',0);
%     ModelResultsSeg(i).coef=mdl.Coefficients;
%     ModelResultsSeg(i).MSE=mdl.MSE;
%     ModelResultsSeg(i).Fitted=mdl.Fitted;
%     ModelResultsSeg(i).rsquared=mdl.Rsquared.Adjusted;
% end
% idx_rsq_Seg=find([ModelResultsSeg.rsquared]>0.1);
% Fighandle=figure;
% set(Fighandle, 'Position', [100, 100, 1400, 900]);
% imagesc(DF(idx_rsq_Seg,:)*100,[0 40]);colormap hot

for i=1:length(PCA_ICA_results)
    PCA_ICA_results(i).nbComp=size(PCA_ICA_results(i).cell_sig,1);
end
Numbers=[0 cumsum([PCA_ICA_results.nbComp])];

Final_data=[];
for c=1:4
    idx=[];idx2=[];
    for Fish=1:6
        idx=[idx find(Numbers(c+((Fish-1)*4))<idx_rsq_Seg & idx_rsq_Seg<Numbers((c+1)+((Fish-1)*4)))];
        idx2=[idx2 find(Numbers(c+((Fish-1)*4))<idx_rsq_ZS & idx_rsq_ZS<Numbers((c+1)+((Fish-1)*4)))];
    end
    Final_data(c).DF=DF(idx_rsq_Seg(idx),:);
    Final_data(c).ZS=ZS(idx_rsq_ZS(idx2),:);
end



temp=idx_rsq_Seg;counter=1;
for i=1:length(filelist)
    name=filelist(i).name;
    imagename=strcat('AVG_',name);
    image=double(imread(imagename));image=image/max(prctile(image,95));image=image*64;
    image=uint8(image);image3=repmat(image,1,1,3);
    idx=find(temp<Numbers(i+1));
    if idx
        temp(idx)=max(Numbers)+1;
        idx=idx_rsq_Seg(idx)-Numbers(i);
        ROIs=zeros(size(image));
        for comp=idx
            ROI=squeeze(PCA_ICA_results(i).ica_segments(comp,:,:));ROI(ROI<max(prctile(ROI,40)))=0;ROI=ROI/max(max(ROI));ROI=ROI*250;
            ROIs=ROIs+ROI;
        end
        image3(:,:,1)=image+uint8(ROIs);
        filename=strcat('AllTraces_',num2str(counter),'_',name);
        imwrite(image3,filename,'tif');
        counter=counter+1;
    end
end

temp=idx_rsq_Seg;counter=1;MeanFinal_data=[];
for i=1:length(filelist)    
    idx=find(temp<Numbers(i+1));    
    if idx
        Traces=mean(DF(temp(idx),:),1);
        temp(idx)=max(Numbers)+1;        
        MeanFinal_data{counter}=Traces;
        counter=counter+1;
    end    
end
temp=cell2mat(MeanFinal_data');

temp=idx_rsq_Seg;counter=1;MeanFinal_data=[];
for i=1:length(filelist)    
    idx=find(temp<Numbers(i+1));    
    if idx
        Traces=mean(AllSegTraces(temp(idx),:),1);
        temp(idx)=max(Numbers)+1;        
        MeanRAW_data{counter}=Traces;
        counter=counter+1;
    end    
end


ZS=zscore(AllSegTraces,1,2);
ModelResultsSeg_ZS=[];
for i=1:length(ZS)
    %mdl=stepwiselm(Regressor',DataSet(i,:),'linear','Criterion','adjrsquared','Intercept',false,'Upper','interactions','Verbose',0);
    mdl=stepwiselm(NewFlow',ZS(i,:),'linear','Criterion','adjrsquared','Upper','linear','Verbose',0);
    ModelResultsSeg_ZS(i).coef=mdl.Coefficients;
    ModelResultsSeg_ZS(i).MSE=mdl.MSE;
    ModelResultsSeg_ZS(i).Fitted=mdl.Fitted;
    ModelResultsSeg_ZS(i).rsquared=mdl.Rsquared.Adjusted;
end
idx_rsq_ZS=find([ModelResultsSeg_ZS.rsquared]>0.2);
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);
imagesc(ZS(idx_rsq_ZS,:),[0 5]);colormap hot

temp=idx_rsq_ZS;counter=1;
for i=1:length(filelist)
    name=filelist(i).name;
    imagename=strcat('AVG_',name);
    image=double(imread(imagename));image=image/max(prctile(image,95));image=image*64;
    image=uint8(image);image3=repmat(image,1,1,3);
    idx=find(temp<Numbers(i+1));
    if idx
        temp(idx)=max(Numbers)+1;
        idx=idx_rsq_ZS(idx)-Numbers(i);
        ROIs=zeros(size(image));
        for comp=idx
            ROI=squeeze(PCA_ICA_results(i).ica_segments(comp,:,:));ROI(ROI<max(prctile(ROI,40)))=0;ROI=ROI/max(max(ROI));ROI=ROI*250;
            ROIs=ROIs+ROI;
        end
        image3(:,:,1)=image+uint8(ROIs);
        filename=strcat('ZSTraces_',num2str(counter),'_',name);
        imwrite(image3,filename,'tif');
        counter=counter+1;
    end
end

temp=idx_rsq_ZS;counter=1;MeanFinal_data=[];
for i=1:length(filelist)    
    idx=find(temp<Numbers(i+1));    
    if idx
        Traces=mean(ZS(temp(idx),:),1);
        temp(idx)=max(Numbers)+1;        
        MeanFinal_dataZS{counter}=Traces;
        counter=counter+1;
    end    
end
temp=cell2mat(MeanFinal_dataZS');