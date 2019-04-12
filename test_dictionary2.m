%% 1.18  ����KSVD��JBHI���걸�ֵ�
% �Ա�KSVD��JBHI���������ͺ�ʱ 
%19��ѹ���� ����100��ECGȡƽ��ֵ
% �ع��㷨�̶�bp����������̶�һ����Ŭ������
% JBHI������QRS��ⷽ�������Ľ�
%%                                
load D_S.mat    %�����ֵ�
load D_NQp.mat  %����QRS�����ֵ�
load D_Qp.mat      %����ʹ�����ֵ䣨12����
N=192; 
q=12;
% threshold=1.5;     %�ж�QRS����ֵ
load('BernoulliSample.mat');   %����һ���̶��������� 1024*1024

for num_CR=1:19          %19��CR
    fprintf('>>ѹ����%d����ʼ����\n',num_CR)
    M=fix((1-num_CR*5*0.01)*N);             %����CR��Ӧ������  CR��95%��90% ~ 5%����ӦM=972��921~51   
    Phi=BernoulliSample(1:M,1:N);                %ÿ��CR��Ӧһ�ֲ�����M����Ӧһ�ֲ�������M*N     
      
                  
                  fprintf('>ѹ����%d����ʼ����\n',num_CR)
                  toc_sum=0;      %��ʱ���͹���
                  PRDsum_KSVD=0;       %��������  ��Ϊ�˼���һ���㷨��һ��CR��  ��100���ĵ�������ƽ��ֵ
                  PRDsum_JBHI=0; 
                  TIMEsum_KSVD=0;
                  TIMEsum_JBHI=0;
                  
                   for num_ECG=1:100;   
                            fprintf('loading ECG%d\n',num_ECG)
                            ecgstr=['ecg',num2str(num_ECG)];                 
                            x=cell2mat(struct2cell(load(ecgstr)));        %ÿ�ε���һ��ecg  ��ecg1~ecg100
                            x_dic=x(1:960);
                         
                            %%   KSVD����
                            L=5;       %�ź�����
                            A=reshape(x_dic,192,5);
                            X=[ ];                      %�ع��� ����С�����ӵ�һά�ź� 
                            tic;
                            for i=1:L
                                x=A(:,i); 
                                %Phi=BernoulliMtx( M,N );                                          
                                y=Phi*x;                                                             
                                Psi=D_S;
                                T=Phi*Psi; 
                                % hat_s=cs_irls(y,T,N);                     %�����ؼ�Ȩ
                                hat_s=cs_bp(y,T,N);                         %��׷��
                                hat_x=real(Psi*hat_s);
                                X=[X hat_x']; 
                            end
                            X_s=X';                                           %X_s��ʹ�ñ�׼�ֵ���ع��ź�
                            TIME_KSVD=toc;
                            PRD_KSVD=norm(x_dic-X_s)/norm(x_dic)*100;     % ��׼�ֵ���ع����
                            
                           %% JBHI����        
region_n=12;            %����������
region_l=N/region_n; %������           
win_w=N/q;    %���ڿ��
Q=cell(1,q);             %Ԫ������q��
X_sub=[];
num=5;
tic;
[R,V,R_no]=QRS_crest(x_dic); 
R_num=length(R);
for i=1:num
%% ���źŽ���QRS�����ж�
    p_start=N*(i-1)+1;
    p_end=N*i;
    x=x_dic(p_start:p_end);
    flag_qrs=0;
        for ii=1:R_num
            if R(ii)>=p_start&&R(ii)<=p_end
                ww=ceil((R(ii)-p_start+1)/region_l);  %wΪ�������ڴ������� λ����1~q�е�һ��
                dictionary=D_Qp{ww};   %ʹ�ú�QRS���ֵ�
                flag_qrs=1;    
            end     
        end
        
        if  flag_qrs==0
            dictionary=D_NQp;    %ʹ�ò���QRS���ֵ�
        end
        
% Phi=BernoulliMtx( M,N );                                          
y=Phi*x;                                                             
Psi=dictionary;
T=Phi*Psi; 
% hat_s=cs_irls(y,T,N);                     %�����ؼ�Ȩ
hat_s=cs_bp(y,T,N);                         %��׷��
hat_x=real(Psi*hat_s);
X_sub=[X_sub hat_x'];                   %�ֶ��ź�ƴ��
end
X_sub=X_sub';
TIME_JBHI=toc;
TIME_JBHI=TIME_JBHI+TIME_KSVD;
PRD_JBHI=norm(x_dic-X_sub)/norm(x_dic)*100;  % ���ֵ���ع����  
%% ͳ����������ʱ��
PRDsum_KSVD=PRDsum_KSVD+PRD_KSVD;
PRDsum_JBHI=PRDsum_JBHI+PRD_JBHI;
TIMEsum_KSVD=TIMEsum_KSVD+TIME_KSVD;
TIMEsum_JBHI=TIMEsum_JBHI+TIME_JBHI;

                    end           %  100��ECG�������
                    
                    PRDaver_KSVD=PRDsum_KSVD/num_ECG;
                    PRDaver_JBHI=PRDsum_JBHI/num_ECG;             
                    TIMEaver_KSVD=TIMEsum_KSVD/num_ECG;
                    TIMEaver_JBHI= TIMEsum_JBHI/num_ECG;
                    
                    Sheet1(1,num_CR)=PRDaver_KSVD;     % KSVD ���
                    Sheet1(2,num_CR)=PRDaver_JBHI;      % JBHI  ���
                    Sheet2(1,num_CR)=TIMEaver_KSVD;     %KSVD ��ʱ     
                    Sheet2(2,num_CR)=TIMEaver_JBHI;     %JBHI  ��ʱ

                    fprintf('>ѹ����%d����ɲ��� \n',num_CR)
                   
          fprintf('>>ѹ����%d��������в��� \n',num_CR)
end          % 19��CR ��������
fprintf('>>>������в��ԣ��������ݴ�����... \n')

xlswrite('KSVD&JBHI_19CR_100ecg.xlsx',Sheet1,'Sheet1','B3');
xlswrite('KSVD&JBHI_19CR_100ecg.xlsx',Sheet2,'Sheet2','B3');
fprintf('>>>>>>>> �ѱ������ݵ����<<<<<<<<\n')
sprintf('>>All Completed<<\n')


% figure(1);
% hold on;
% plot(Sheet1(1,:),'r')
% plot(Sheet1(2,:),'g') 
% legend('KSVD','JBHI')  
% 
% 
% figure(2);
% plot(Sheet2(1,:),'r')
% plot(Sheet2(2,:),'g') 
% legend('KSVD','JBHI')  





