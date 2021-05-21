%% 1.18  ����KSVD��JBHI���걸�ֵ�
% �Ա�KSVD��JBHI���������ͺ�ʱ 
%19��ѹ���� ����100��ECGȡƽ��ֵ
% �ع��㷨�̶�bp����������̶�һ����Ŭ������
% JBHI������QRS��ⷽ�������Ľ�
%%                                
load data/D_S.mat    %�����ֵ�
load data/D_NQp.mat  %����QRS�����ֵ�
load data/D_Qp.mat      %����ʹ�����ֵ䣨12����
N=192; 
L=5;       %�ź�����
q=12;
% threshold=1.5;     %�ж�QRS����ֵ
load('BernoulliSample.mat');   %����һ���̶��������� 1024*1024

for num_CR=1:19          %19��CR
    fprintf('>>starting CR%d\n',num_CR)
    M=fix((1-num_CR*5*0.01)*N);             %����CR��Ӧ������  CR��95%��90% ~ 5%����ӦM=972��921~51   
    Phi=BernoulliSample(1:M,1:N);                %ÿ��CR��Ӧһ�ֲ�����M����Ӧһ�ֲ�������M*N     
    toc_sum=0;      %��ʱ���͹���
    PRDsum_KSVD=0;       %��������  ��Ϊ�˼���һ���㷨��һ��CR��  ��100���ĵ�������ƽ��ֵ
    TIMEsum_KSVD=0;
    
    for num_ECG=1:10;   
        fprintf('loading ECG%d\n',num_ECG)
        ecgstr=['data/ecg',num2str(num_ECG),'.mat'];                 
        x=cell2mat(struct2cell(load(ecgstr)));        %ÿ�ε���һ��ecg  ��ecg1~ecg100
        x_dic=x(1:N*L);
        
        %%   KSVD����
        A=reshape(x_dic,N,L);
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

        PRDsum_KSVD=PRDsum_KSVD+PRD_KSVD;
        TIMEsum_KSVD=TIMEsum_KSVD+TIME_KSVD;
                            


    end           %  100��ECG�������
                    
    PRDaver_KSVD=PRDsum_KSVD/num_ECG;
    TIMEaver_KSVD=TIMEsum_KSVD/num_ECG;

    fprintf('CR%d PRD: %f\n', num_CR, PRDaver_KSVD)   
end          % 19��CR ��������

fprintf('>>All Completed<<\n')


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





