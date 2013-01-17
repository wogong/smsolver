
%Structure Mechanics
%matrix displacement method
%连续梁内力计算程序
%wogong

%程序中变量含义汇总：
%NE：单元总数
%BL：单元杆长
%EI：单元抗弯刚度
%JD：单元定位向量
%FM：单元固端弯矩
%i：单元线刚度
%N：单元转角未知量总数
%P(N)：先表示直接结点荷载，后表示等效结点荷载
%P0(2,NE):单元非结点荷载
%KE(N,N)：结构整体刚度矩阵
%DELTA(N)：单元未知位移
%FJ(2,NE)：单元杆端弯矩
%I：循环计数变量


%%%%%%%%%%%%%%%%输入原始数据并输出%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear,clc;
!echo 欢迎使用连续梁内力计算程序，请按任意键继续。
pause;
!echo 输入原始数据。
!echo 请严格按照如下算例格式输入数据：
!echo 单元总数 Number of element:6
!echo 单元杆长 Element length:[4 6 6 8 4 6]
!echo 单元抗弯刚度 Bending stiffness:sym('[EI 1.5*EI EI 2*EI EI 1.5*EI]')
!echo 单元定位向量 Element localization vector:[0 1 2 3 4 5;1 2 3 4 5 0]
!echo 单元非结点荷载：[0 0.5 0.5 0 0.5 0;0 12 8 0 10 0]
!echo 直接节点力矩 Straight joint moment:[0;-8;0;10;0]
pause;
NE=input('单元总数 Number of element:');
BL=input('单元杆长 Element length:');
EI=input('单元抗弯刚度 Bending stiffness:');
JD=input('单元定位向量 Element localization vector:');
P0=input('单元非结点荷载：');%其中P0(1,NE)作为NE单元荷载类型识别量，若为均布荷载则此值为1，若为集中荷载则此值表示集中荷载作用位置距杆端1的距离与单元总长的比值，若此单元上无非结点荷载，则此值为0；
                                  %其中P0(2,NE)作为NE单元荷载大小标识量，若为均布荷载则此值表示均布荷载集度，若为集中荷载则此值表示集中荷载大小，若此单元无何载，则此值不会起作用，默认为0.
for I=1:1:NE
    if(P0(2,I)~=0)
        if(P0(1,I)==1)  %判断是否为均布荷载
            FM(:,I)=[-P0(2,I)*BL(I)^2/12;P0(2,I)*BL(I)^2/12];%均布荷载固端弯矩计算公式
        else            %判断是否为集中荷载
            FM(:,I)=[-P0(2,I)*P0(1,I)*(1-P0(1,I))^2*BL(I);P0(2,I)*P0(1,I)^2*(1-P0(1,I))*BL(I)];%集中荷载固端弯矩计算公式
        end
    else
        FM(:,I)=[0;0];  %此单元无非结点荷载
    end
end%固端弯矩生成
N=JD(1,NE); %节点转角未知量总数
if JD(2,NE)~=0
    N=JD(2,NE);
end
P=input('直接节点力矩 Straight joint moment:');%非零编码结点荷载，维数N
i=EI./BL;%单元线刚度
!echo 所计算连续梁的原始数据如下：
!echo 单元总数:
NE
!echo 单元杆长：
BL
!echo 单元抗弯刚度：
EI
!echo 单元定位向量：
JD
!echo 单元固端弯矩：
FM
!echo 结点转角未知量总数：
N
!echo 直接结点力矩：
P
!echo 按任意键开始计算。
pause;


%%%%%%%%%%%%%%%%形成结点荷载向量%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%利用单元定位向量生成结点荷载向量
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


%%%%%%%%%%%%%%%%形成整体刚度矩阵%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%利用单元定位向量装配整体刚度矩阵
KE=zeros(N,N);%单元刚度矩阵置零
KE=sym([KE]);%定义为符号矩阵
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


%%%%%%%%%%%%%%%%解方程%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

DELTA=KE\P;


%%%%%%%%%%%%%%%%求杆端弯矩并输出%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FJ=zeros(2,N);%单元杆端弯矩初始化
FJ=sym(FJ);%定义为符号矩阵
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
end%计算杆端弯矩
!echo 按任意键输出计算结果。
pause;
!echo 结点转角：
DELTA
!echo 单元杆端弯矩计算结果：
FJ
