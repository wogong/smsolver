
%Structure Mechanics
%matrix displacement method
%�����������������
%wogong

%�����б���������ܣ�
%NE����Ԫ����
%BL����Ԫ�˳�
%EI����Ԫ����ն�
%JD����Ԫ��λ����
%FM����Ԫ�̶����
%i����Ԫ�߸ն�
%N����Ԫת��δ֪������
%P(N)���ȱ�ʾֱ�ӽ����أ����ʾ��Ч������
%P0(2,NE):��Ԫ�ǽ�����
%KE(N,N)���ṹ����նȾ���
%DELTA(N)����Ԫδ֪λ��
%FJ(2,NE)����Ԫ�˶����
%I��ѭ����������


%%%%%%%%%%%%%%%%����ԭʼ���ݲ����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear,clc;
!echo ��ӭʹ��������������������밴�����������
pause;
!echo ����ԭʼ���ݡ�
!echo ���ϸ�������������ʽ�������ݣ�
!echo ��Ԫ���� Number of element:6
!echo ��Ԫ�˳� Element length:[4 6 6 8 4 6]
!echo ��Ԫ����ն� Bending stiffness:sym('[EI 1.5*EI EI 2*EI EI 1.5*EI]')
!echo ��Ԫ��λ���� Element localization vector:[0 1 2 3 4 5;1 2 3 4 5 0]
!echo ��Ԫ�ǽ����أ�[0 0.5 0.5 0 0.5 0;0 12 8 0 10 0]
!echo ֱ�ӽڵ����� Straight joint moment:[0;-8;0;10;0]
pause;
NE=input('��Ԫ���� Number of element:');
BL=input('��Ԫ�˳� Element length:');
EI=input('��Ԫ����ն� Bending stiffness:');
JD=input('��Ԫ��λ���� Element localization vector:');
P0=input('��Ԫ�ǽ����أ�');%����P0(1,NE)��ΪNE��Ԫ��������ʶ��������Ϊ�����������ֵΪ1����Ϊ���к������ֵ��ʾ���к�������λ�þ�˶�1�ľ����뵥Ԫ�ܳ��ı�ֵ�����˵�Ԫ���޷ǽ����أ����ֵΪ0��
                                  %����P0(2,NE)��ΪNE��Ԫ���ش�С��ʶ������Ϊ�����������ֵ��ʾ�������ؼ��ȣ���Ϊ���к������ֵ��ʾ���к��ش�С�����˵�Ԫ�޺��أ����ֵ���������ã�Ĭ��Ϊ0.
for I=1:1:NE
    if(P0(2,I)~=0)
        if(P0(1,I)==1)  %�ж��Ƿ�Ϊ��������
            FM(:,I)=[-P0(2,I)*BL(I)^2/12;P0(2,I)*BL(I)^2/12];%�������ع̶���ؼ��㹫ʽ
        else            %�ж��Ƿ�Ϊ���к���
            FM(:,I)=[-P0(2,I)*P0(1,I)*(1-P0(1,I))^2*BL(I);P0(2,I)*P0(1,I)^2*(1-P0(1,I))*BL(I)];%���к��ع̶���ؼ��㹫ʽ
        end
    else
        FM(:,I)=[0;0];  %�˵�Ԫ�޷ǽ�����
    end
end%�̶��������
N=JD(1,NE); %�ڵ�ת��δ֪������
if JD(2,NE)~=0
    N=JD(2,NE);
end
P=input('ֱ�ӽڵ����� Straight joint moment:');%�����������أ�ά��N
i=EI./BL;%��Ԫ�߸ն�
!echo ��������������ԭʼ�������£�
!echo ��Ԫ����:
NE
!echo ��Ԫ�˳���
BL
!echo ��Ԫ����նȣ�
EI
!echo ��Ԫ��λ������
JD
!echo ��Ԫ�̶���أ�
FM
!echo ���ת��δ֪��������
N
!echo ֱ�ӽ�����أ�
P
!echo ���������ʼ���㡣
pause;


%%%%%%%%%%%%%%%%�γɽ���������%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%���õ�Ԫ��λ�������ɽ���������
for I=1:1:NE
		I1=JD(1,I);
		I2=JD(2,I);
		if I1~=0
            P(I1)=P(I1)-FM(1,I);
       end
		if I2~=0
            P(I2)=P(I2)-FM(2,I);
        end
end


%%%%%%%%%%%%%%%%�γ�����նȾ���%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%���õ�Ԫ��λ����װ������նȾ���
KE=zeros(N,N);%��Ԫ�նȾ�������
KE=sym([KE]);%����Ϊ���ž���
for I=1:1:NE
    I1=JD(1,I);
	I2=JD(2,I);
    if I1~=0
        KE(I1,I1)=vpa(KE(I1,I1)+4*i(I));
        if I2~=0
            KE(I1,I2)=vpa(KE(I1,I2)+2*i(I));
            KE(I2,I1)=vpa(KE(I2,I1)+2*i(I));
        end
    end
    if I2~=0
        KE(I2,I2)=vpa(KE(I2,I2)+4*i(I));
    end
end


%%%%%%%%%%%%%%%%�ⷽ��%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DELTA=KE\P;


%%%%%%%%%%%%%%%%��˶���ز����%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FJ=zeros(2,N);%��Ԫ�˶���س�ʼ��
FJ=sym(FJ);%����Ϊ���ž���
for I=1:1:NE
    I1=JD(1,I);
	I2=JD(2,I);
    DZ=sym('[0;0]');
    if I1~=0
        DZ(1)=DELTA(I1);
    end
    if I2~=0
        DZ(2)=DELTA(I2);
    end
    FJ(1,I)=4*i(I)*DZ(1)+2*i(I)*DZ(2)+FM(1,I);
    FJ(2,I)=2*i(I)*DZ(1)+4*i(I)*DZ(2)+FM(2,I);
end%����˶����
!echo ������������������
pause;
!echo ���ת�ǣ�
DELTA
!echo ��Ԫ�˶���ؼ�������
FJ