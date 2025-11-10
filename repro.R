
 n_trs<-200;n_voxels<-50;n_trials<-40;tr<-2
 onsets<-seq(10,n_trs*tr-20,length.out=n_trials)
 conditions<-rep(c('A','B'),each=n_trials/2)
 event_table<-data.frame(onset=onsets,condition=conditions,duration=1)
 print('A')
 X_sim<-model.matrix(~ condition - 1,data=event_table)
 print('B')
 signal<-X_sim[,1]-X_sim[,2]
 print('C')
 hrf_basis<-fmrihrf::gen_hrf(1:20,hrf=fmrihrf::getHRF('spmg2'))
 print('D')
 hrf_vec<-as.numeric(hrf_basis %*% c(1,0))
 print('E')
 signal_conv<-stats::convolve(signal,rev(hrf_vec),type='open')[1:n_trs]
 print('F')
 Y_train<-matrix(rnorm(n_trs*n_voxels,sd=1),n_trs,n_voxels)
 print('G')
 true_pattern<-c(rep(1,n_voxels/2),rep(0,n_voxels/2))
 print('H')
 for(v in 1:n_voxels){Y_train[,v]<-Y_train[,v]+signal_conv*true_pattern[v]*0.5}
 print('I')
 print(dim(Y_train))
