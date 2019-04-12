%% 1.18  ����С����
% �Ա�52��С��������� �ͺ�ʱ
% �ع��㷨�̶�bp����������̶���Ŭ������
% 19��ѹ���� ����100��ECGȡƽ��ֵ

%%                                
N=1024;   
load('BernoulliSample.mat');   %����һ���̶��������� 1024*1024
for num_CR=1:19          %19��CR
    fprintf('>>ѹ����%d����ʼ����\n',num_CR)
    M=fix((1-num_CR*5*0.01)*N);             %����CR��Ӧ������  CR��95%��90% ~ 5%����ӦM=972��921~51   
    Phi=BernoulliSample(1:M,:);                %ÿ��CR��Ӧһ�ֲ�����M����Ӧһ�ֲ�������M*N
   % Phi=BernoulliMtx( M,N );        
         for num_wave=1:52        % 52��С���� 
                   fprintf('>ѹ����%d-С����%d����ʼ����\n',num_CR,num_wave)
                   toc_sum=0;      %��ʱ���͹���
                   PRD_sum=0;       %��������  ��Ϊ�˼���һ���㷨��һ��CR��  ��100���ĵ�������ƽ��ֵ
                    
                   for num_ECG=1:100;    %100��ECG
                            fprintf('loading ECG%d\n',num_ECG)
                            ecgstr=['ecg',num2str(num_ECG)];                 
                            x=cell2mat(struct2cell(load(ecgstr)));        %ÿ�ε���һ��ecg  ��ecg1~ecg100
                            y=Phi*x; 
                            tic       %��ʼ��ʱ   �ؽ��źź�ʱ
                            [ww]=dwtmtxx( N, num_wave,5);    %С������ѡ��  �ں���dwtmtxx.m��
                            Psi=[ww];                                           
                            T=Phi*Psi'; 
                            hat_s=cs_bp(y,T,N);                     %��׷���㷨BP
                            hat_x=real(Psi'*hat_s);
                            time_end=toc;   %������ʱ   �ؽ��źź�ʱ
                            toc_sum=toc_sum+time_end;          %�ۼ�ʱ��
                            PRD=norm(x-hat_x)/norm(x)*100;     %���㱾�����
                            PRD_sum=PRD_sum+PRD;                %�ۼ����      
                    end           %  100��ECG�������
                    
                    fprintf('>ѹ����%d-С����%d����ɲ��� \n',num_CR,num_wave)
                    PRD_aver=PRD_sum/num_ECG;  %100��ECG�����ƽ��ֵ
                    time_aver=toc_sum/num_ECG;   %100��ECG�ĺ�ʱƽ��ֵ
                    A(num_wave,num_CR)=PRD_aver;     %A�� �������
                    B(num_wave,num_CR)=toc_sum;     %B�� ��ʱ����         
          end        %һ��ѹ���ȶ�Ӧ��С������������
    
          fprintf('>>ѹ����%d��������в��� \n',num_CR)
end          % 19��CR ��������
fprintf('>>>������в��ԣ��������ݴ�����... \n')

xlswrite('sparse_19CR_100ecg.xlsx',A,'Sheet1','B3');
xlswrite('sparse_19CR_100ecg.xlsx',B,'Sheet2','B3');
fprintf('>>>>>>>> �ѱ������ݵ����<<<<<<<<\n')
sprintf('>>All Completed<<\n')







