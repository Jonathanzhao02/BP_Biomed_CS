%% 1.18  �����ع��㷨��
%�Ա�5���㷨�����ͺ�ʱ
%ϡ����̶�С����db2����������̶���Ŭ������ 
%19��ѹ����  ��100��ECG���ȡƽ��ֵ

%%                                
N=1024;   
load('BernoulliSample.mat');   %����һ���̶��������� 1024*1024
for num_CR=1:19          %19��CR
    fprintf('>>ѹ����%d����ʼ����\n',num_CR)
    M=fix((1-num_CR*5*0.01)*N);             %����CR��Ӧ������  CR��95%��90% ~ 5%����ӦM=972��921~51   
    Phi=BernoulliSample(1:M,:);                %ÿ��CR��Ӧһ�ֲ�����M����Ӧһ�ֲ�������M*N
   % Phi=BernoulliMtx( M,N );        
         for num_algo=1:5       % 5���ع��㷨 
                   fprintf('>ѹ����%d-�㷨%d����ʼ����\n',num_CR,num_algo)
                   toc_sum=0;      %��ʱ���͹���
                   PRD_sum=0;       %��������  ��Ϊ�˼���һ���㷨��һ��CR��  ��100���ĵ�������ƽ��ֵ
                    for num_ECG=1:100;    %100��ECG
                            fprintf('loading ECG%d\n',num_ECG)
                            ecgstr=['ecg',num2str(num_ECG)];                 
                            x=cell2mat(struct2cell(load(ecgstr)));        %ÿ�ε���һ��ecg  ��ecg1~ecg100
                            y=Phi*x;                                                     
                            [ww]=dwtmtx( N,'db2',5);     %ѡ��С����
                            Psi=[ww];                                  
                            T=Phi*Psi'; 
                            tic       %��ʼ��ʱ   �ؽ��źź�ʱ
                                        switch num_algo                                         %�ع��㷨ѡ��  K����
                                            case 1
                                                 hat_s=cs_omp(y,T,N);
                                                 hat_x=real(Psi'*hat_s.');                  
                                            case 2
                                                hat_s=cs_bp(y,T,N);    
                                                hat_x=real(Psi'*hat_s);        
                                            case 3
                                                 hat_s=cs_cosamp(y,T,N);
                                                 hat_x=real(Psi'*hat_s.');         
                                            case 4
                                                 hat_s=cs_irls(y,T,N);  
                                                 hat_x=real(Psi'*hat_s);       
                                            case 5
                                                 hat_s=cs_sp(y,T,N); 
                                                 hat_x=real(Psi'*hat_s.');       
                                        end
                                time_end=toc;   %������ʱ   �ؽ��źź�ʱ
                                toc_sum=toc_sum+time_end;          %�ۼ�ʱ��
                                PRD=norm(x-hat_x)/norm(x)*100;     %���㱾�����
                                PRD_sum=PRD_sum+PRD;                %�ۼ����      
                                                                   
                    end           %  100��ECG�������
                    fprintf('>ѹ����%d-�㷨%d����ɲ��� \n',num_CR,num_algo)
                    PRD_aver=PRD_sum/num_ECG;  %100��ECG�����ƽ��ֵ
                    time_aver=toc_sum/num_ECG;   %100��ECG�ĺ�ʱƽ��ֵ

                    A(num_algo,num_CR)=PRD_aver;     %A�� �������
                    B(num_algo,num_CR)=toc_sum;     %B�� ��ʱ����
                   
          end        %�ع��㷨��������
          fprintf('>>ѹ����%d����ɲ��� \n',num_CR)
end          % 19��CR ��������
fprintf('>>>������в��ԣ��������ݴ�����... \n')

xlswrite('algo_19CR_100ecg.xlsx',A,'Sheet1','B3');  %�����ع��㷨���Աȵı��
xlswrite('algo_19CR_100ecg.xlsx',B,'Sheet2','B3');  %�����ع��㷨��ʱ�Աȵı��
fprintf('>>>>>>>> �ѱ������ݵ����<<<<<<<<\n')
sprintf('>>All Completed<<\n')







