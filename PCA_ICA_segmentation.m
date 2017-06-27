filelist=dir('fish*.tif');
for File=1:length(filelist)
    [mixedsig, mixedfilters, CovEvals, covtrace, movm, movtm] = CellsortPCA(filelist(File).name,[],100,1,'C:\Temp\CROP');
    [ica_sig, ica_filters, ica_A, numiter] = CellsortICA(mixedsig, mixedfilters, CovEvals, [1:size(mixedsig,1)],  0.8,size(mixedsig,1),randn(size(mixedsig,1), size(mixedsig,1)),1e-6,50000);
    %PCA_ICA_results(File).ica_sig=ica_sig;
    %PCA_ICA_results(File).ica_filters=ica_filters;
    %[PCA_ICA_results(File).ica_segments, PCA_ICA_results(File).segmentlabel, PCA_ICA_results(File).segcentroid] = CellsortSegmentation(ica_filters, 2, 2, 20, 0)
    [PCA_ICA_results.ROIs, segmentlabel, PCA_ICA_results.SegCentroid] = CellsortSegmentation(ica_filters, 2, 2, 20, 0);    
    PCA_ICA_results.Cell_sig = CellsortApplyFilter(filelist(File).name, PCA_ICA_results.ROIs);
    save(strcat(filelist(File).name(1:length(filelist(File).name)-4),'.mat'),'PCA_ICA_results','-v7.3');
end
clearvars mixedsig ica_sig ica_filters ica_A numiter mixedsig mixedfilters CovEvals covtrace movm movtm
 