HDBscanCluster=zeros(max(hdbscan_labels),size(GoodCalcium,2));
for i=1:max(hdbscan_labels)
    idx=find(hdbscan_labels==i);
    HDBscanCluster(i,:)=mean(zscore(GoodCalcium(idx,:),1,2),1);
end

ZS_spikes=zscore(GoodSpikes,1,2);
temp=ZS_spikes;
temp(temp<2)=0;
temp=boolean(temp);
temp2=sum(temp,2);

ZS=zscore(GoodCalcium,1,2);

GoodCalcium_backup=GoodCalcium;
GoodCalcium(GoodCalcium<1)=0;

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



x = linspace(1,size(Cmap,2),size(Cmap,2));
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);
counter=1;xplot=floor(sqrt(length(GoodBetas_select)));yplot=ceil(length(GoodBetas_select)/xplot);
for i=GoodBetas_select    
    NumberOfCells=length(find(idxKmeans==i));
    %subplot(5,1,counter);plot(x,Cmap(i,:),x,Model_DF(i).Fitted);title(num2str(NumberOfCells))
    subplot(xplot,yplot,counter);plot(x,Cmap(i,:));title(num2str(NumberOfCells))
    xlim([0 size(Cmap,2)])
    counter=counter+1;
end


x = linspace(0.2,size(Cmap,2)/5,size(Cmap,2));
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1400, 900]);
i=GoodBetas_select(1);NumberOfCells=length(find(idxKmeans==i));
subplot(3,4,1);plot(x,Cmap(i,:));title(num2str(NumberOfCells));axis([0 131 -0.04 0.1]);rectangle('FaceColor','r','Position',[10 -0.04 10 0.005]);rectangle('FaceColor','r','Position',[50 -0.04 10 0.005]);rectangle('FaceColor','r','Position',[90 -0.04 10 0.005]);rectangle('FaceColor','b','Position',[30 -0.04 10 0.005]);rectangle('FaceColor','b','Position',[70 -0.04 10 0.005]);rectangle('FaceColor','b','Position',[110 -0.04 10 0.005]);
i=GoodBetas_select(9);NumberOfCells=length(find(idxKmeans==i));
subplot(3,4,2);plot(x,Cmap(i,:));title(num2str(NumberOfCells));axis([0 131 -0.04 0.1]);rectangle('FaceColor','r','Position',[10 -0.04 10 0.005]);rectangle('FaceColor','r','Position',[50 -0.04 10 0.005]);rectangle('FaceColor','r','Position',[90 -0.04 10 0.005]);rectangle('FaceColor','b','Position',[30 -0.04 10 0.005]);rectangle('FaceColor','b','Position',[70 -0.04 10 0.005]);rectangle('FaceColor','b','Position',[110 -0.04 10 0.005]);
i=GoodBetas_select(11);NumberOfCells=length(find(idxKmeans==i));
subplot(3,4,4);plot(x,Cmap(i,:));title(num2str(NumberOfCells));axis([0 131 -0.04 0.1]);rectangle('FaceColor','r','Position',[10 -0.04 10 0.005]);rectangle('FaceColor','r','Position',[50 -0.04 10 0.005]);rectangle('FaceColor','r','Position',[90 -0.04 10 0.005]);rectangle('FaceColor','b','Position',[30 -0.04 10 0.005]);rectangle('FaceColor','b','Position',[70 -0.04 10 0.005]);rectangle('FaceColor','b','Position',[110 -0.04 10 0.005]);
i=GoodBetas_select(4);NumberOfCells=length(find(idxKmeans==i));
subplot(3,4,5);plot(x,Cmap(i,:));title(num2str(NumberOfCells));axis([0 131 -0.04 0.1]);rectangle('FaceColor','r','Position',[10 -0.04 10 0.005]);rectangle('FaceColor','r','Position',[50 -0.04 10 0.005]);rectangle('FaceColor','r','Position',[90 -0.04 10 0.005]);rectangle('FaceColor','b','Position',[30 -0.04 10 0.005]);rectangle('FaceColor','b','Position',[70 -0.04 10 0.005]);rectangle('FaceColor','b','Position',[110 -0.04 10 0.005]);
i=GoodBetas_select(10);NumberOfCells=length(find(idxKmeans==i));
subplot(3,4,6);plot(x,Cmap(i,:));title(num2str(NumberOfCells));axis([0 131 -0.04 0.1]);rectangle('FaceColor','r','Position',[10 -0.04 10 0.005]);rectangle('FaceColor','r','Position',[50 -0.04 10 0.005]);rectangle('FaceColor','r','Position',[90 -0.04 10 0.005]);rectangle('FaceColor','b','Position',[30 -0.04 10 0.005]);rectangle('FaceColor','b','Position',[70 -0.04 10 0.005]);rectangle('FaceColor','b','Position',[110 -0.04 10 0.005]);
i=GoodBetas_select(6);NumberOfCells=length(find(idxKmeans==i));
subplot(3,4,7);plot(x,Cmap(i,:));title(num2str(NumberOfCells));axis([0 131 -0.04 0.1]);rectangle('FaceColor','r','Position',[10 -0.04 10 0.005]);rectangle('FaceColor','r','Position',[50 -0.04 10 0.005]);rectangle('FaceColor','r','Position',[90 -0.04 10 0.005]);rectangle('FaceColor','b','Position',[30 -0.04 10 0.005]);rectangle('FaceColor','b','Position',[70 -0.04 10 0.005]);rectangle('FaceColor','b','Position',[110 -0.04 10 0.005]);
i=GoodBetas_select(2);NumberOfCells=length(find(idxKmeans==i));
subplot(3,4,8);plot(x,Cmap(i,:));title(num2str(NumberOfCells));axis([0 131 -0.04 0.1]);rectangle('FaceColor','r','Position',[10 -0.04 10 0.005]);rectangle('FaceColor','r','Position',[50 -0.04 10 0.005]);rectangle('FaceColor','r','Position',[90 -0.04 10 0.005]);rectangle('FaceColor','b','Position',[30 -0.04 10 0.005]);rectangle('FaceColor','b','Position',[70 -0.04 10 0.005]);rectangle('FaceColor','b','Position',[110 -0.04 10 0.005]);
i=GoodBetas_select(5);NumberOfCells=length(find(idxKmeans==i));
subplot(3,4,9);plot(x,Cmap(i,:));title(num2str(NumberOfCells));axis([0 131 -0.04 0.1]);rectangle('FaceColor','r','Position',[10 -0.04 10 0.005]);rectangle('FaceColor','r','Position',[50 -0.04 10 0.005]);rectangle('FaceColor','r','Position',[90 -0.04 10 0.005]);rectangle('FaceColor','b','Position',[30 -0.04 10 0.005]);rectangle('FaceColor','b','Position',[70 -0.04 10 0.005]);rectangle('FaceColor','b','Position',[110 -0.04 10 0.005]);
i=GoodBetas_select(7);NumberOfCells=length(find(idxKmeans==i));
subplot(3,4,10);plot(x,Cmap(i,:));title(num2str(NumberOfCells));axis([0 131 -0.04 0.1]);rectangle('FaceColor','r','Position',[10 -0.04 10 0.005]);rectangle('FaceColor','r','Position',[50 -0.04 10 0.005]);rectangle('FaceColor','r','Position',[90 -0.04 10 0.005]);rectangle('FaceColor','b','Position',[30 -0.04 10 0.005]);rectangle('FaceColor','b','Position',[70 -0.04 10 0.005]);rectangle('FaceColor','b','Position',[110 -0.04 10 0.005]);
i=GoodBetas_select(8);NumberOfCells=length(find(idxKmeans==i));
subplot(3,4,12);plot(x,Cmap(i,:));title(num2str(NumberOfCells));axis([0 131 -0.04 0.1]);rectangle('FaceColor','r','Position',[10 -0.04 10 0.005]);rectangle('FaceColor','r','Position',[50 -0.04 10 0.005]);rectangle('FaceColor','r','Position',[90 -0.04 10 0.005]);rectangle('FaceColor','b','Position',[30 -0.04 10 0.005]);rectangle('FaceColor','b','Position',[70 -0.04 10 0.005]);rectangle('FaceColor','b','Position',[110 -0.04 10 0.005]);

GoodBetas=GoodBetas_select([1 9 4 10 6 2 5 7 8]);
counter=1;GoodClusters=[];
for i=GoodBetas
    idx=find(idxKmeans==i);
    GoodClusters(counter,:)=mean(GoodCalcium(idx,:),1);
    counter=counter+1;
end