!***************
module NumKind
!***************
   implicit none
   integer (kind(1)),parameter :: ikind = kind(1), rkind = kind(0.D0)
   real (rkind),parameter :: Zero = 0._rkind, One  = 1._rkind, Two   = 2._rkind, &
                             Three= 3._rkind, Four = 4._rkind, Five  = 5._rkind, &
                             Six  = 6._rkind, Seven= 7._rkind, Eight = 8._rkind, &
                             Nine = 9._rkind, Ten  =10._rkind, Twelve=12._rkind
end module NumKind

!***************
module TypeDef
!**************
   use NumKind
   implicit none

   type :: typ_Joint !�ڵ�����
      real (rkind)      :: x,y
      integer (ikind)   :: GDOF(3)
   end type typ_Joint

   type :: typ_Element !��Ԫ����
      integer (ikind)   :: JointNo(2),GlbDOF(6)
      real (rkind)      :: Length,CosA,SinA,EI,EA,mass
   end type typ_Element

   type :: typ_JointLoad !�ڵ��������
      integer (ikind)   :: JointNo,LodDOF
      real (rkind)      :: LodVal
   end type typ_JointLoad

   type :: typ_ElemLoad  !��Ԫ��������
      integer (ikind)   :: ElemNo,Indx
      real (rkind)      :: Pos,LodVal
   end type typ_ElemLoad

   contains

   !===================================
   subroutine SetElemProp(Elem,Joint) !�趨��Ԫ����
   !===================================
      type (typ_Element),intent(in out) :: Elem(:)
      type (typ_Joint),intent(in) :: Joint(:)
      integer(ikind) :: i,NElm
      real(rkind) :: x1,x2,y1,y2
      NElm=size(Elem,dim=1)
      do i=1,NElm
         x1=Joint(Elem(i)%JointNo(1))%x
         x2=Joint(Elem(i)%JointNo(2))%x
         y1=Joint(Elem(i)%JointNo(1))%y
         y2=Joint(Elem(i)%JointNo(2))%y
         Elem(i)%Length=sqrt((x1-x2)**2+(y1-y2)**2)
         Elem(i)%CosA=(x2-x1)/Elem(i)%Length
         Elem(i)%SinA=(y2-y1)/Elem(i)%Length
         Elem(i)%GLbDOF(1:3)=Joint(Elem(i)%JointNo(1))%GDOF(:)
         Elem(i)%GLbDOF(4:6)=Joint(Elem(i)%JointNo(2))%GDOF(:)
      end do
      return
   end subroutine SetElemProp

   !======================================
   subroutine TransMatrix (ET, CosA,SinA)
   !======================================
      real (rkind),intent(out) :: ET(:,:)
      real (rkind),intent(in) :: CosA,SinA
      ET=Zero
      ET(1,1:2) = (/ CosA, SinA /)
      ET(2,1:2) = (/-SinA, CosA /)
      if (size(ET,1) > 2) ET(3,3) = One
      if (size(ET,1) > 3) ET(4:6,4:6) = ET(1:3,1:3)
      return
   end subroutine TransMatrix

end module TypeDef

!***************
module BandMat
!*************
   use NumKind
   use TypeDef,only : typ_Element  ! ���ø�ģ���е�typ_Element
   implicit none

   private  !  Ĭ�����е����ݺ͹���Ϊ˽�У���ǿ��װ��
   public   :: SetMatBand, DelMatBand, VarBandSolv

   type,public :: typ_Kcol
      real (rkind),pointer :: row(:)
   end type typ_Kcol

   contains

   !==================================
   subroutine SetMatBand (Kcol, Elem)
   !==================================
      ! ...[6-4-2]

      type (typ_KCol),intent(in out) :: Kcol(:)
      type (typ_Element),intent(in)  :: Elem(:)
      integer (ikind)                :: minDOF,ELocVec(6)
      integer (ikind)                :: ie,j,NGlbDOF,NElem
      integer (ikind)                :: row1(size(Kcol,dim=1))
      ! row1���Զ����飬�ӳ���������Զ��ͷ��ڴ�ռ�
      NGlbDOF = size(Kcol,1)
      NElem   = size(Elem,1)
      row1    = NGlbDOF       !  ����ʼ����Ϊ����
      ! ȷ������ʼ���룬��������row1(:)��
      do ie=1,NElem
         ELocVec = Elem(ie)%GlbDOF
         minDOF  = minval (ELocVec,mask = ELocVec > 0)
         where (ELocVec > 0)
            row1(ELocVec) = min(row1(ELocVec), minDOF)
         end where
      end do
      ! Ϊ���еİ�������ռ䲢��ʼ��
      do j=1,NGlbDOF
         allocate ( Kcol(j)%row(row1(j):j) )
         Kcol(j)%row = Zero   ! ����
      end do
      return
   end subroutine SetMatBand

   !===================================
   subroutine DelMatBand (Kcol)
   !===================================
      !...[6-5-5]
      type (typ_KCol), intent(in out)  :: Kcol(:)
      integer (ikind)                  :: j,NGlbDOF
      NGlbDOF = size(Kcol,1)
      do j=1,NGlbDOF
         deallocate ( Kcol(j)%row )
      end do
      return
   end subroutine DelMatBand

   !===========================================
   subroutine VarBandSolv ( Disp, Kcol, GLoad )
   !===========================================
      type (typ_KCol),intent(in out)    :: Kcol( : )
      real (rkind),intent(out)          :: Disp( : )
      real (rkind),intent(in)           :: GLoad( : )

      integer (ikind) :: i,j,row1j,row_1,NCol
      real (rkind)     :: s,Diag(size(Kcol,1))
      ! ...[6 - 5 - 2]
      NCol = size(Kcol,1)
      Diag(:)=(/(Kcol(j)%row(j),j=1,NCol)/)
      do j=2,NCol
         row1j=lbound(Kcol(j)%row,1)
         do i=row1j,j-1
            row_1=max(row1j,lbound(Kcol(i)%row,1))
            ! s=sum(Diag(row1j:i-1)*Kcol(i)%row(row_1:i-1)*Kcol(j)%row(row_1:i-1))
            s=sum(Kcol(i)%row(row_1:i-1)*Kcol(j)%row(row_1:i-1))
            ! Kcol(j)%row(i)=(Kcol(j)%row(i)-s)/Diag(i)
            Kcol(j)%row(i)=Kcol(j)%row(i)-s
         end do
         Kcol(j)%row(row1j:j-1)=Kcol(j)%row(row1j:j-1)/Diag(row1j:j-1)
         s=sum(Diag(row1j:j-1)*Kcol(j)%row(row1j:j-1)**2)
         Diag(j)=Diag(j)-s
      end do
      Disp(:)=GLoad(:)
      ! ...[ 6 - 5 - 3�ڵĴ��룺����GP��ΪDisp ]
      do j=2,NCol
         row1j=lbound(Kcol(j)%row,1)
         Disp(j)=Disp(j)-sum(Kcol(j)%row(row1j:j-1)*Disp(row1j:j-1))
      end do
      ! ...[ 6 - 5 - 4�ڵĴ��룺����GP��ΪDisp ]
      Disp(:)=Disp(:)/Diag(:)
      do j=NCol,1,-1
         row1j=lbound(Kcol(j)%row,1)
         Disp(row1j:j-1)=Disp(row1j:j-1)-Disp(j)*Kcol(j)%row(row1j:j -1)
      end do
      return
   end subroutine VarBandSolv
   
end module BandMat

!***************
module DispMethod
!****************
   use NumKind
   use TypeDef
   use BandMat
   implicit none

   contains

   !==================================================
   subroutine SolveDisp (Disp, Elem,Joint,JLoad,ELoad)
   !==================================================
      real (rkind),intent(out)         :: Disp(:)
      type (typ_Element),intent(in)  :: Elem(:)
      type (typ_Joint),intent(in)  :: Joint(:)
      type (typ_JointLoad),intent(in)  :: JLoad(:)
      type (typ_ElemLoad),intent(in)  :: ELoad(:)
      integer (ikind) ::NGlbDOF
      type (typ_Kcol),pointer :: Kcol(:)
      real (rkind),pointer :: Gload(:)

      NGlbDOF = size(Disp,1)
      allocate (GLoad(NGlbDOF))
      allocate (Kcol(NGlbDOF))

      call SetMatBand (Kcol, Elem) !�õ��նȾ������
      call GLoadVec (GLoad, Elem,JLoad,ELoad,Joint) !�õ������������
      call GStifMat (Kcol, Elem) !�õ�����նȾ���
      call VarBandSolv (Disp, Kcol,GLoad) !���λ��
      call DelMatBand(Kcol)
      return
   end subroutine SolveDisp

   !==================================================
   subroutine GStifMat (Kcol, Elem) !����նȾ���
   !================================
      type (typ_Kcol),intent(in out)  :: Kcol(:)
      type (typ_Element),intent(in)   :: Elem(:)
      ! ...[6-4-3]
      integer (ikind)                 :: ie,j,JGDOF,NElem
      real (rkind)                    :: EK(6,6),ET(6,6) 
      integer (ikind)                 :: ELocVec(6) 
      NElem=size(Elem,1) 
      do ie=1,NElem 
         call EStifMat(EK,Elem(ie)%Length,Elem(ie)%EI,Elem(ie)%EA) 
         call TransMatrix(ET,Elem(ie)%CosA,Elem(ie)%SinA) 
         EK=matmul(transpose(ET),matmul(EK,ET)) 
         ELocVec=Elem(ie)%GlbDOF 
         do j=1,6 
            JGDOF=ELocVec(j) 
            if(JGDOF==0) cycle 
            where (ELocVec<=JGDOF.and.ELocVec>0) 
               Kcol(JGDOF)%row(ELocVec)=Kcol(JGDOF)%row(ELocVec)+EK(:,j) 
            end where 
         end do 
      end do 
      return
   end subroutine GStifMat

   !==================================================
   subroutine GLoadVec (GLoad,Elem,JLoad,ELoad,Joint) !�����������
   !===================================================
      real (rkind),intent(out) :: Gload(:)
      type (typ_Element),  intent(in)  :: Elem(:)
      type (typ_Joint),  intent(in)  :: Joint(:)
      type (typ_JointLoad),  intent(in)  :: JLoad(:)
      type (typ_ElemLoad),  intent(in)  :: ELoad(:)
      
      integer (ikind) ::i,DOFj,NJLoad,NELoad
      real (rkind) :: F0(6),ET(6,6)

      NJLoad=size(JLoad,dim=1)
      NELoad=size(ELoad,dim=1)

      Gload=Zero

      do i=1,NJLoad
         DOFj=Joint(JLoad(i)%JointNo)%GDOF(JLoad(i)%LodDOF)
         GLoad(DOFj)=GLoad(DOFj)+JLoad(i)%LodVal
      end do 
      
      do i=1,NELoad 
         call EFixendF(F0,ELoad(i)%Indx,ELoad(i)%Pos,ELoad(i)%LodVal,Elem(ELoad(i)%ElemNo))
         call TransMatrix(ET,Elem(ELoad(i)%ElemNo)%CosA,Elem(ELoad(i)%ElemNo)%SinA)
         where (Elem(ELoad(i)%ElemNo)%GlbDOF>0)
            GLoad(Elem(ELoad(i)%ElemNo)%GlbDOF)=GLoad(Elem(ELoad(i)%ElemNo)%GlbDOF)-matmul(transpose(ET),F0(:))
         end where         
      end do
      return
   end subroutine GLoadVec
   
   !==================================================
   subroutine EStifMat (EK,ELen,EI,EA) !��Ԫ�նȾ���
   !===================================
      real (rkind),intent(out) :: EK(:,:)
      real (rkind),intent(in) :: Elen,EI,EA
      real (rkind) :: a1,a2,a3,a4
      a1=EA/ELen
      a2=Twelve*EI/(ELen**3)
      a3=Six*EI/(ELen**2)
      a4=Four*EI/ELen
      EK=Zero
      EK(1,1)=a1
      EK(4,1)=-a1
      EK(1,4)=-a1
      EK(4,4)=a1
      EK(2,2)=a2
      EK(5,5)=a2
      EK(2,5)=-a2
      EK(5,2)=-a2
      EK(2,3)=a3
      EK(3,2)=a3
      EK(2,6)=a3
      EK(6,2)=a3
      EK(3,5)=-a3
      EK(5,3)=-a3
      EK(5,6)=-a3
      EK(6,5)=-a3
      EK(3,3)=a4
      EK(6,6)=a4
      EK(3,6)=a4/Two
      EK(6,3)=a4/Two
      return
   end subroutine EStifMat

   !==================================================
   subroutine EFixendF (F0,Indx,a,q,Elem) !��Ԫ�̶���
   !======================================
      real (rkind),intent(in out) :: F0(:)
      real (rkind),intent(in) :: a,q
      integer (ikind),intent(in) :: Indx
      type (typ_Element),intent(in) :: Elem
      real (rkind) :: l,ra,rb
      l=Elem%Length
      F0=Zero
      select case(Indx)
        case(1)
            ra=a*l
            F0(2)=(-q*ra/Two)*(Two-Two*(ra**2)/(l**2)+(ra**3)/(l**3))
            F0(5)=(-q*(ra**3)/(Two*(l**2)))*(2-ra/l)
            F0(3)=(-q*(ra**2)/Twelve)*(Six-Eight*ra/l+Three*(ra**2)/(l**2))
            F0(6)=(q*(ra**3)/(Twelve*l))*(Four-Three*ra/l)
        case(2)
            ra=a*l
            rb=l-ra
            F0(2)=(-q*(rb**2)/(l**2))*(One+Two*ra/l)
            F0(5)=(-q*(ra**2)/(l**2))*(One+Two*rb/l)
            F0(3)=-q*ra*(rb**2)/(l**2)
            F0(6)=q*(ra**2)*rb/(l**2)
        case(3)
            if(a<One/Two) then
                F0(1)=Elem%EA*q/l
                F0(4)=-F0(1)
            else
                F0(1)=-Elem%EA*q/l
                F0(4)=-F0(1)
            end if
        case(4)
            if(a<One/Two) then
                F0(2)=Twelve*Elem%EI*q/(l**3)
                F0(3)=Six*Elem%EI*q/(l**2)
                F0(5)=-F0(2)
                F0(6)=F0(3)
            else
                F0(2)=-Twelve*Elem%EI*q/(l**3)
                F0(3)=-Six*Elem%EI*q/(l**2)
                F0(5)=-F0(2)
                F0(6)=F0(3)
            end if
        end select
        return
    end subroutine EFixendF

   !==================================================
   subroutine ElemDisp (EDisp,Disp,Elem)     ! Good: ie is not necessary !
   !=========================================
      real (rkind),intent(out) :: EDisp(:)
      real (rkind),intent(in) :: Disp(:)
      type (typ_Element),intent(in) :: Elem
      integer (ikind) :: i
      do i=1,6
         if(Elem%GlbDOF(i)==0) then
             EDisp(i)=Zero
         else
             EDisp(i)=Disp(Elem%GlbDOF(i))
         end if
      end do
      return
   end subroutine ElemDisp

   !==================================================
   subroutine ElemForce (EForce,ie,Disp,Elem,ELoad)
   !==================================================
   real (rkind),intent(out) :: EForce(:)
   real (rkind),intent(in) :: Disp(:)
   type (typ_Element),intent(in) :: Elem(:)
   type (typ_ElemLoad),intent(in) ::ELoad(:)
   integer (ikind),intent(in) :: ie

   real (rkind) ::EK(6,6),ET(6,6),EP(6),F0(6),EDisp(6)
   integer (ikind) :: NELoad,i

   NELoad=size(ELoad,1)
   EForce=Zero
   EP(:)=Zero
   do i=1,NELoad
      if(ELoad(i)%ElemNo==ie) then
          call EFixendF(F0,ELoad(i)%Indx,Eload(i)%Pos,ELoad(i)%LodVal,Elem(ELoad(i)%ElemNo))
          EP(:)=EP(:)-F0(:)
      end if
   end do
   call EStifMat(EK,Elem(ie)%Length,Elem(ie)%EI,Elem(ie)%EA)
   call TransMatrix(ET,Elem(ie)%CosA,Elem(ie)%SinA)
   call ElemDisp(EDisp,Disp,Elem(ie))
   EForce(:)=matmul(EK,matmul(ET,EDisp))-EP(:)
   return
   end subroutine ElemForce

end module DispMethod

!***************
module Vibration
!***************
   use NumKind
   use TypeDef 
   use BandMat
   implicit none
   real (rkind) :: PI=3.141592653589793

   contains

   !==================================================
   subroutine GetFreq(Freq,Elem,Kcol,SFreq,NFreq,Tol)!���ַ���ȡƵ��
   !==================================================
      real (rkind),intent(in out) :: Freq(:)
      type (typ_Element),intent(in) :: Elem(:)
      type (typ_Kcol),intent(in out) :: Kcol(:)
      integer (ikind),intent(in) :: SFreq,NFreq
      real (rkind),intent(in) :: Tol  
      real (rkind) :: freq_l,freq_u,freqm
      integer (ikind) :: k    
      integer (ikind) :: J01,J02,J_l,J_u,Jk,J0,Jm,MFreq

      MFreq=1   
      freq_l=One
      freq_u=Ten        
      do k=SFreq,SFreq+NFreq-1
         if(MFreq>1) then !��Ƶ������1
            MFreq=MFreq-1
            Freq(k-SFreq+1)=Freq(k-SFreq)
            cycle
         end if
         if(k>SFreq)  freq_l=freq(k-SFreq)
         do
            call Calculate_J0(J01,freq_l,Elem)
            call Calculate_Jk(Jk,freq_l,Kcol,Elem)
            J_l=J01+Jk
            if(J_l<k) then 
                freq_u=freq_l !���Ч��
                exit
            end if
            freq_l = freq_l/Two
         end do
         do
            call Calculate_J0(J02,freq_u,Elem)
            call Calculate_Jk(Jk,freq_u,Kcol,Elem)
            J_u=J02+Jk
            if(J_u>=k) exit
            freq_l=freq_u
            freq_u=freq_u*Two
         end do
         do !���ַ�
            freqm=(freq_l+freq_u)/Two
            call Calculate_J0(J0,freqm,Elem)
            call Calculate_Jk(Jk,freqm,Kcol,Elem)
            Jm=J0+Jk
            if(Jm>=k)then
               freq_u=freqm
            else
               freq_l=freqm
            end if
            if((freq_u-freq_l)<Tol*(One+freq_l)) then
               call Calculate_J0(J0,freq_u,Elem)
               call Calculate_Jk(Jk,freq_u,Kcol,Elem)
               Jm=J0+Jk
               MFreq=Jm-k+1 !�õ���Ƶ��
               exit
            end if
         end do
         Freq(k+1-SFreq)=(freq_u+freq_l)/Two
      end do
      return
   end subroutine GetFreq

   !==================================================
   subroutine Mode(Elem,Kcol,SFreq,NFreq,Tol)!��ȡ����
   !==================================================
      type (typ_Element),intent(in) :: Elem(:)
      type (typ_Kcol),intent(in out) :: Kcol(:)
      real (rkind),intent(in) :: Tol  
      integer (ikind),intent(in) :: SFreq,NFreq
      real (rkind),allocatable :: delta0(:),delta(:),GLoad(:),diag(:),GLoad2(:),Freq(:),Mdelta(:,:) !��¼��Ƶ����
      integer (ikind),allocatable :: MFreq(:) !��¼��Ƶ�ѳ��ִ���
      real (rkind) :: de,s,freq0,alpha,EDisp(6)
      integer (ikind) :: NGlbDOF,NElem,i,j,row1j,row_1,NCol,k,t,r
    
      NGlbDOF=size(Kcol,1)
      NCol=size(Kcol,1)
      NElem=size(Elem,1)
      r=min(NFreq,NGlbDOF)

      allocate(delta0(NGlbDOF))
      allocate(delta(NGlbDOF))
      allocate(GLoad(NGlbDOF))
      allocate(diag(NGlbDOF))
      allocate(MFreq(NFreq+1))
      allocate(GLoad2(NGlbDOF))
      allocate(Freq(NFreq))

      MFreq=Zero !��Ƶ����ʼ��

      call GetFreq(Freq,Elem,Kcol,SFreq,NFreq,Tol)
      
      do i=1,NFreq-1 !������Ƶ��Ӧ���ѳ��ִ���
         if((Freq(i+1)-Freq(i))<Two*Tol*(1+Freq(i+1))) MFreq(i+1)=MFreq(i)+1
      end do

      allocate(Mdelta(NGlbDOF,sum(MFreq))) 
      write(8,*) NFreq

      do k=1,r
         call random_number(delta(:)) !�������������ʼ����������
         delta0=-One

         do
            if(MFreq(k)>0) then !��Ƶ����������
                do t=1,MFreq(k)
                    call EGLoadVec(GLoad,Elem,delta,Freq(k))
                    call EGLoadVec(GLoad2,Elem,Mdelta(1:NGlbDOF,t),Freq(k))
                    alpha=dot_product(Mdelta(1:NGlbDOF,t),GLoad)/dot_product(Mdelta(1:NGlbDOF,t),GLoad2)
                    delta=delta-Mdelta(1:NGlbDOF,t)*alpha
                end do
            end if

            !���ݵ���
            call SetMatBand(Kcol,Elem)
            call GStifMat(Kcol,Elem,Freq(k))
            
            ! ...[6 - 5 - 2]
            diag(:)=(/(Kcol(j)%row(j),j=1,NCol)/) 
            do j=2,NCol
               row1j=lbound(Kcol(j)%row,1)
               do i=row1j,j-1
                  row_1=max(row1j,lbound(Kcol(i)%row,1))
                  s=sum(diag(row_1:i-1)*Kcol(i)%row(row_1:i-1)*Kcol(j)%row(row_1:i-1))
                  Kcol(j)%row(i)=(Kcol(j)%row(i)-s)/diag(i)
               end do
               s=sum(diag(row1j:j-1)*Kcol(j)%row(row1j:j-1)**2)
               diag(j)=diag(j)-s
            end do      
            
            do j=2,NCol
               row1j=lbound(Kcol(j)%row,1) 
               GLoad(j)=GLoad(j)-sum(Kcol(j)%row(row1j:j-1)*GLoad(row1j:j-1)) 
            end do 
            GLoad(:)=GLoad(:)/Diag(:)

            do j=NCol,1,-1 
               row1j=lbound(Kcol(j)%row,1) 
               GLoad(row1j:j-1)=GLoad(row1j:j-1)-GLoad(j)*Kcol(j)%row(row1j:j-1) 
            end do
            
            de=GLoad(1)
            do i=2,NGlbDOF !��׼��
               if(abs(GLoad(i))>abs(de)) de=GLoad(i)     
            end do
            do i=1,NGlbDOF
               GLoad(i)=GLoad(i)/de
            end do
            
            if(dot_product(GLoad-delta0,GLoad-delta0)<Tol*Tol) then
               if(MFreq(k+1)>0) Mdelta(1:NGlbDOF,MFreq(k+1))=GLoad !��¼��Ƶ���ͣ����ں���������
               exit
            end if
            
            delta0=GLoad
            freq0=Freq(k)
            Freq(k)=Freq(k)-1/de !Newton���Ƶ��
            if(Freq(k)>freq0+Tol*(One+freq0)) Freq(k)=freq0 
            if(Freq(k)<freq0-Tol*(One+freq0)) Freq(k)=freq0
         end do
         
         write(8,*) Freq(k),0                  
         do t=1,NElem                      
            call ElemDisp2(EDisp,GLoad,Elem(t))
            write(8,*) EDisp       
         end do
      end do

      EDisp=0
      
      do k=r+1,NFreq
         write(8,*) Freq(k),0
         do t=1,NElem
            write(8,*) EDisp
         end do
      end do
      
      deallocate(diag)
      deallocate(delta0)
      deallocate(GLoad)
      call DelMatBand (Kcol)
      return

   end subroutine Mode

   !==================================================
   subroutine Calculate_J0(J0,freq,Elem)!�������freq�ĵ�Ԫ�̶�Ƶ����J0
   !==================================================
      integer (ikind),intent(out) :: J0
      real (rkind),intent(in) :: freq
      type (typ_Element),intent(in) :: Elem(:)
      integer (ikind) :: NElem,Ja,Jb,i,j
      real (rkind) :: nu,lambda,sg

      NElem=size(Elem,dim=1)
      J0=0
      do i=1,NElem
         nu    = freq*Elem(i)%Length*sqrt(Elem(i)%mass/Elem(i)%EA)
         Ja    = int(nu/PI)
         lambda= Elem(i)%Length*sqrt(sqrt(freq**2*Elem(i)%mass/Elem(i)%EI))
         sg    = sign(One,One-cosh(lambda)*cos(lambda))
         j     = int(lambda/PI)
         Jb    = j-(1-(-1)**j*sg)/2
         J0    = J0+Ja+Jb
      end do
      return
   end subroutine Calculate_J0

   !==================================================
   subroutine Calculate_Jk(Jk,freq,Kcol,Elem)!�������freq������Ƶ����Jk
   !==================================================
      integer (ikind),intent(out) :: Jk
      real (rkind),intent(in) :: freq
      type (typ_Kcol),intent(in out) :: Kcol(:)
      type (typ_Element),intent(in) :: Elem(:)
      integer (ikind) :: NGlbDOF
      real (rkind),allocatable :: diag(:)
      integer (ikind) :: i,j,row1j,row_1,NCol
      real (rkind) :: s
      
      NGlbDOF=size(Kcol,dim=1)
      allocate(diag(NGlbDOF))

      call SetMatBand(Kcol,Elem)
      call GStifMat(Kcol,Elem,freq)

      ! ...[6 - 5 - 2]
      NCol=size(Kcol,1)
      diag(:)=(/(Kcol(j)%row(j),j=1,NCol)/) 
      do j=2,NCol
         row1j=lbound(Kcol(j)%row,1)
         do i=row1j,j-1
            row_1=max(row1j,lbound(Kcol(i)%row,1))
            s=sum(diag(row_1:i-1)*Kcol(i)%row(row_1:i-1)*Kcol(j)%row(row_1:i-1))
            Kcol(j)%row(i)=(Kcol(j)%row(i)-s)/diag(i)
         end do
         s=sum(diag(row1j:j-1)*Kcol(j)%row(row1j:j-1)**2)
         diag(j)=diag(j)-s
      end do      

      Jk=count(mask=diag<Zero,dim=1)
      deallocate(diag)
      return
   end subroutine Calculate_Jk
   
   !==================================================
   subroutine GStifMat(Kcol,Elem,freq)!���嶯���նȾ���
   !==================================================
      type (typ_Kcol),intent(in out) :: Kcol(:)
      type (typ_Element),intent(in) :: Elem(:)
      real (rkind),intent(in) :: freq
      integer (ikind) :: ie,j,JGDOF,NElem
      real (rkind) :: EK(6,6),ET(6,6)
      integer (ikind) :: ELocVec(6)

      NElem=size(Elem,1)
      do ie=1,NElem
         call EStifMat(EK,Elem,ie,freq) 
         call TransMatrix(ET,Elem(ie)%CosA,Elem(ie)%SinA)
         EK=matmul(transpose(ET),matmul(EK,ET)) 
         ELocVec=Elem(ie)%GlbDOF 
         do j=1,6
            JGDOF=ELocVec(j) 
            if(JGDOF==0) cycle 
            where (ELocVec<=JGDOF.and.ELocVec>0) 
               Kcol(JGDOF)%row(ELocVec)=Kcol(JGDOF)%row(ELocVec)+EK(:,j) 
            end where 
         end do 
      end do
      return
   end subroutine  GStifMat

   !==================================================
      subroutine EGLoadVec(GLoad,Elem,delta,freq)!��Ч��������������������ݵ���
   !==================================================
         real (rkind),intent(out) :: GLoad(:)
         type (typ_Element),intent(in) :: Elem(:)
         real (rkind),intent(in) :: freq
         real (rkind),intent(in) :: delta(:)
         integer (ikind) :: ie,NElem,i
         real (rkind) :: F0(6),EK(6,6)

         GLoad(:)=Zero
         NElem=size(Elem,1)
         do ie=1,NElem 
            F0=Zero
            do i=1,6
               if(Elem(ie)%GlbDOF(i)>0) F0(i)=delta(Elem(ie)%GlbDOF(i))
            end do
            call DEStifMat(EK,Elem,ie,freq) !��Ԫ�����նȾ���ĵ�������
            where (Elem(ie)%GlbDOF>0)
               GLoad(Elem(ie)%GlbDOF)=GLoad(Elem(ie)%GlbDOF)+matmul(EK,F0)
            end where
         end do
         return
    end subroutine EGLoadVec

   !==================================================
   subroutine EStifMat(EK,Elem,ie,freq)!��Ԫ�����նȾ���
   !==================================================
      real (rkind),intent(out) :: EK(6,6)
      type (typ_Element),intent(in) :: Elem(:)
      integer (ikind),intent(in) :: ie
      real (rkind),intent(in) :: freq
      real (rkind) :: EI,EA,Length,m
      real (rkind) :: phi
      real (rkind) :: B1,B2,T,R,Q,H,S,C
      real (rkind) :: nu,lambda
      integer (ikind) :: i,j

      EK    = Zero        
      EI    = Elem(ie)%EI
      EA    = Elem(ie)%EA
      Length= Elem(ie)%Length
      m     = Elem(ie)%mass
      nu    = freq*Length*sqrt(m/EA)
      lambda= Length*(freq**2*m/EI)**0.25

      phi = One-cosh(lambda)*cos(lambda)
      B1  =nu/tan(nu)
      B2  =nu/sin(nu)
      T   = lambda**3*(sin(lambda)*cosh(lambda)+cos(lambda)*sinh(lambda))/phi
      R   = lambda**3*(sin(lambda)+sinh(lambda))/phi
      Q   = lambda**2*(sin(lambda)*sinh(lambda))/phi
      H   = lambda**2*(cosh(lambda)-cos(lambda))/phi
      S   = lambda*(sin(lambda)*cosh(lambda)-cos(lambda)*sinh(lambda))/phi
      C   = lambda*(sinh(lambda)-sin(lambda))/phi

      EK(1,1) = B1*EA/Length
      EK(1,4) = -B2*EA/Length
      EK(2,2) = T*EI/Length**3
      EK(2,3) = Q*EI/Length**2
      EK(2,5) = -R*EI/Length**3
      EK(2,6) = H*EI/Length**2
      EK(3,3) = S*EI/Length
      EK(3,5) = -H*EI/Length**2
      EK(3,6) = C*EI/Length
      EK(4,4) = B1*EA/Length
      Ek(5,5) = T*EI/Length**3
      EK(5,6) = -Q*EI/Length**2
      EK(6,6) = S*EI/Length

      do j=1,6
         do i=j+1,6
            EK(i,j)=EK(j,i)
         end do
      end do
      return
   end subroutine EStifMat

   !==================================================
   subroutine DEStifMat(EK,Elem,ie,freq)!��Ԫ�����նȾ���ĵ�������
   !==================================================
      real (rkind),intent(out) :: EK(6,6)
      type (typ_Element),intent(in) :: Elem(:)
      integer (ikind),intent(in) :: ie
      real (rkind),intent(in) :: freq
      real (rkind) :: EI,EA,Length,m
      real (rkind) :: phi
      real (rkind) :: T,R,Q,H,S,C,B11,B21,T1,R1,Q1,H1,S1,C1
      real (rkind) :: nu,lambda
      integer (ikind) :: i,j

      EK    = Zero        
      EI    = Elem(ie)%EI
      EA    = Elem(ie)%EA
      Length= Elem(ie)%Length
      m     = Elem(ie)%mass
      nu    = freq*Length*sqrt(m/EA)
      lambda= Length*(freq**2*m/EI)**0.25

      phi  = One-cosh(lambda)*cos(lambda)
      B11=nu/freq/tan(nu)-nu*nu/freq*(One/tan(nu)/tan(nu)+1)
      B21=nu/freq/sin(nu)-nu*nu/freq*cos(nu)/sin(nu)/sin(nu)
      T = lambda**3*(sin(lambda)*cosh(lambda)+cos(lambda)*sinh(lambda))/phi
      R = lambda**3*(sin(lambda)+sinh(lambda))/phi
      Q = lambda**2*(sin(lambda)*sinh(lambda))/phi
      H = lambda**2*(cosh(lambda)-cos(lambda))/phi
      S = lambda*(sin(lambda)*cosh(lambda)-cos(lambda)*sinh(lambda))/phi
      C = lambda*(sinh(lambda)-sin(lambda))/phi

      T1=One/Two/freq*lambda*(Three*T/lambda-T*S/lambda+Two*lambda**3*cos(lambda)*cosh(lambda)/phi)
      R1=One/Two/freq*lambda*(Three*R/lambda-R*S/lambda+lambda**3*(cos(lambda)+cosh(lambda))/phi)
      Q1=One/Two/freq*lambda*(Two*Q/lambda-Q*S/lambda+T/lambda)
      H1=One/Two/freq*lambda*(Two*H/lambda-H*S/lambda+R/lambda)
      S1=One/Two/freq*lambda*(S/lambda-S*S/lambda+Two*Q/lambda)
      C1=One/Two/freq*lambda*(C/lambda-C*S/lambda+H/lambda)

      EK(1,1) = B11*EA/Length
      EK(1,4) = -B21*EA/Length
      EK(2,2) = T1*EI/Length**3
      EK(2,3) = Q1*EI/Length**2
      EK(2,5) = -R1*EI/Length**3
      EK(2,6) = H1*EI/Length**2
      EK(3,3) = S1*EI/Length
      EK(3,5) = -H1*EI/Length**2
      EK(3,6) = C1*EI/Length
      EK(4,4) = B11*EA/Length
      Ek(5,5) = T1*EI/Length**3
      EK(5,6) = -Q1*EI/Length**2
      EK(6,6) = S1*EI/Length

      do j=1,6
         do i=j+1,6
            EK(i,j)=EK(j,i)
         end do
      end do
      return
   end subroutine DEStifMat

   !==================================================
   subroutine ElemDisp2(EDisp,Disp,Elem)
   !==================================================
        real (rkind),intent(out) :: EDisp(:)  
        real (rkind),intent(in) :: Disp(:)    
        type (typ_Element),intent(in) :: Elem    
        integer (ikind) :: i

        do i=1,6 
            if(Elem%GlbDOF(i)==0) then 
                EDisp(i)=Zero 
            else 
                EDisp(i)=Disp(Elem%GlbDOF(i)) 
            end if 
        end do 
        return 
    end subroutine ElemDisp2 

end module Vibration

!==============================================
program SM_90            ! main prog
!==============================================
   use NumKind           ! numeric precision module
   use TypeDef           ! data type defination module
   use DispMethod        ! displacement method module
   use Vibration         ! vibration module 
   implicit none
   type (typ_Element),allocatable :: Elem(:)
   type (typ_Joint),allocatable :: Joint(:)
   type (typ_JointLoad),allocatable :: JLoad(:)
   type (typ_ElemLoad),allocatable :: ELoad(:)
   integer (ikind) :: ProbType
   integer (ikind) :: NElem,NJoint,NGlbDOF,NJLoad,NELoad
   integer (ikind) :: NFreq,SFreq
   real (rkind) :: Tol
   real (rkind),pointer :: Disp(:)
   type (typ_Kcol),allocatable :: Kcol(:)

   call InputData()                                             ! internal sub in below
   call SetElemProp(Elem,Joint)                                 ! contained in TypeDef
   if(Probtype==1) call SolveDisp(Disp, Elem,Joint,JLoad,ELoad) ! contained in DispMethod
   call OutputResults()                                         ! internal sub in below
 
   stop

   contains

   !--------------------------
   subroutine InputData ()
   !--------------------------
      integer (ikind) :: i,ie
      open(3,file="SM90.IPT",status="OLD",position="REWIND",action="READ")
      read(3,*) ProbType
      if(ProbType==2) read(3,*) NFreq,SFreq,Tol
      read(3,*) NElem,NJoint,NGlbDOF,NJLoad,NELoad
      allocate (Joint(NJoint))
      allocate (Elem(NElem))
      allocate (JLoad(NJLoad))
      allocate (ELoad(NELoad))
      allocate (Disp(NGlbDOF))
      allocate (Kcol(NGlbDOF))
      Disp=Zero
      read(3,*) (Joint(i),i=1,NJoint)
      if(ProbType==1) read(3,*) (Elem(ie)%JointNo,Elem(ie)%EA,Elem(ie)%EI,ie=1,NElem)
      if(ProbType==2) read(3,*) (Elem(ie)%JointNo,Elem(ie)%EA,Elem(ie)%EI,Elem(ie)%mass,ie=1,NElem)
      if(NJLoad>0) read(3,*) (JLoad(i),i=1,NJLoad)
      if(NELoad>0) read(3,*) (ELoad(i),i=1,NELoad)
      close(3)
      return
   end subroutine InputData

   !--------------------------
   subroutine OutputResults ()
   !--------------------------
      integer(ikind) :: i,ie
      real (rkind) :: EDisp(6),EForce(6),EDispLoad(6),ET(6,6)
      open(8,file="SMCAI90.OUT",status="REPLACE",position="REWIND",action="WRITE")
      write(8,*) 10,0
      if(ProbType==1) then
          do ie=1,NElem
            call ElemDisp(EDisp,Disp,Elem(ie))
            EDispLoad=Zero
            call TransMatrix(ET,Elem(ie)%CosA,Elem(ie)%SinA)
            do i=1,NELoad 
                    if(ELoad(i)%ElemNo==ie)then 
                        if(ELoad(i)%Indx==3)then 
                            if(ELoad(i)%pos==0)then !ʼ���ƶ�
                                EDispLoad(1)=EDispLoad(1)+ELoad(i)%LodVal 
                            else !ĩ���ƶ�
                                EDispLoad(4)=EDispLoad(4)+ELoad(i)%LodVal 
                            end if 
                        else if(ELoad(i)%Indx==4)then 
                            if(ELoad(i)%pos==0)then 
                                EDispLoad(2)=EDispLoad(2)+ELoad(i)%LodVal 
                            else 
                                EDispLoad(5)=EDispLoad(5)+ELoad(i)%LodVal 
                            end if 
                        end if 
                    end if 
                end do 
                EDispLoad=matmul(transpose(ET),EDispLoad)
                EDisp=EDisp+EDispLoad
                write(8,*) Edisp(1),Edisp(2),Edisp(3),Edisp(4),Edisp(5),Edisp(6)
            end do
            do ie=1,NElem
                call ElemForce (EForce,ie,Disp,Elem,ELoad)
                write(8,*) -EForce(1),EForce(2),EForce(3),EForce(4),EForce(5),EForce(6)
            end do
      else if(ProbType==2) then
          call Mode(Elem,Kcol,SFreq,NFreq,Tol)
      end if
      deallocate(JLoad) 
      deallocate(Joint) 
      deallocate(Elem) 
      deallocate(ELoad) 
      deallocate(Disp) 
      deallocate(Kcol) 
      close(8)
      return
   end subroutine OutputResults

end program SM_90
