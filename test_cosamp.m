% Runs CoSaMP on MIT-BIH Arrhythmia dataset (mitbih.mat)
% All data (.mat) files will be found in a location to be specified later

N=1024;
DSET="mitbih";

load('BernoulliSample.mat');

fprintf('>> Loading %s\n', DSET);
load(DSET + ".mat");

% Generate basis
[ww]=dwtmtx( N,'db2',5);
Psi=[ww];

for CR=[2, 4, 8]
  logfile = fopen(sprintf('./results/%s/%s_rec_%d.txt', DSET, DSET, CR), 'w');
  fprintf('>> Starting CR%d\n',CR)
  M=fix(N / CR);

  data_rec = zeros(size(data));

  % Sensing matrix
  Phi=BernoulliSample(1:M,:);

  toc_sum=0;
  PRD_sum=0;
  SNR_sum=0;
  for num_ECG=1:10;
    fprintf('loading ECG%d\n',num_ECG);
    x=data(num_ECG,:)';
    y=Phi*x;
    T=Phi*Psi';
    tic
    hat_s=cs_cosamp(y,T,N);
    hat_x=real(Psi'*hat_s.');
    time_end=toc;
    toc_sum=toc_sum+time_end;
    PRD=norm(x-hat_x)/norm(x)*100;
    PRD_sum=PRD_sum+PRD;
    SNR_sum=SNR_sum+snr(x,x-hat_x);
  end
  fprintf('>> CR%d, finished \n',CR);
  data_rec(num_ECG,:)=hat_x;
  PRD_aver=PRD_sum/num_ECG;
  SNR_aver=SNR_sum/num_ECG;
  time_aver=toc_sum/num_ECG;

  A(CR)=PRD_aver;
  B(CR)=toc_sum;
  C(CR)=SNR_aver;

  fprintf('>> CR%d SNR: %f PRD: %f Time: %f\n', CR, SNR_aver, PRD_aver, time_aver);
  fprintf(logfile, '>> CR%d SNR: %f PRD: %f Time: %f\n', CR, SNR_aver, PRD_aver, time_aver);

  % Save reconstruction results
  save(sprintf('./results/%s/%s_rec_%d.mat', DSET, DSET, CR), 'data_rec');
  fclose(logfile);
end

fprintf('>> All Completed <<\n');
