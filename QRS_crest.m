%% Tompkins��ַ����QRS�����ҵ�R�����壬����ÿ������ĺ�����

function [R,V,R_no]=QRS_crest(x)
%RΪ��������ꣻVΪ���������ꣻR_noΪQRS������

m=length(x);   % =SAMPLES2READ
y1 = [];%һ�ײ��
for n = 1:(m-1)
    y1(n) = x(n+1) - x(n);
end
y1(m)=y1(m-1);

y2 = [];%���ײ��
for n = 2:(m-1)
    y2(n) = x(n+1)-2*x(n)+x(n-1);
end
y2(1)=y2(2);
y2(m)=y2(m-1);

y=[]; %Tompkins��ַ���ȡ�ĵ��źŵ�һ������ײ�ֵ�ƽ������ΪQRS��Ⱥ�����ǵ������ź�
for k=1:1:m;
    y(k)=y1(k)^2+y2(k)^2;
end

 th=0.02;  % ����ֵ,�ɻ��ݼ��,������ֵ--������R����ֵ����
 QRS=find(y-th>0);
 qrs={};  % QRS��Ⱥ
 
 r=1;
 j=2;
 qrs{1}(1)=QRS(1);%qrs{r}���ڼ�¼R��ʱ������
 
for i=2:1:length(QRS);
    if QRS(i)-QRS(i-1)<=100;%ʱ���������100������Ϊû�е��¸�QRS���Σ������QRS��ʱ���������qrs{r}
        qrs{r}(j)=QRS(i);
            j=j+1;
    else
        qrs{r+1}(1)=QRS(i);%ʱ���������100������Ϊ���¸�QRS���Σ������QRS��ʱ���������qrs{r+1}
        r=r+1;
        R_no=r;%��¼QRS������
        j=2;   
    end 
end

Rsite=zeros(R_no,1);%QRS����Ϊ������10�У�ÿ��QRS���μ�¼10������


for r=1:1:R_no;
    n=1;
    for i=2:1:length(qrs{r});
    if qrs{r}(i)-qrs{r}(i-1)>1;%����˼�¼���ʱ�������ǰһ����1����
        Rsite(r)=qrs{r}(i-1)+2;%��ǰһ����¼���ʱ�������2��¼��Rsite
        break;
    end
    end
    if Rsite(r)==0;
        Rsite(r)=min(qrs{r});
    end
end

R=zeros(R_no,1);%R������ĺ�����
V=zeros(R_no,1);%R�������������

for i=1:R_no;
[V(i),R(i)]=max(x(Rsite(i)-20:Rsite(i)+20));
R(i)=R(i)+Rsite(i)-21;
end
