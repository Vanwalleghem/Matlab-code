function Draw_GoodBetas(GoodBetas,Cmap_BF,idxKmeans)
Fighandle=figure;
set(Fighandle, 'Position', [100, 100, 1500, 1000]);
x = linspace(1,size(Cmap_BF,1),size(Cmap_BF,1));
counter=1;xplot=floor(sqrt(length(GoodBetas)));yplot=ceil(length(GoodBetas)/xplot);
colors = distinguishable_colors(length(GoodBetas));
for i=GoodBetas
    NumberOfCells=length(find(idxKmeans==i));
    subplot(xplot,yplot,counter);plot(x,Cmap_BF(:,i),'Color',colors(counter,:));title(num2str(NumberOfCells));
    xlim([0 size(Cmap_BF,1)])
    counter=counter+1;
end
clearvars NumberOfCells