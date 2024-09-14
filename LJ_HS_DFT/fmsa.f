	SUBROUTINE CP(T,RHOB,DP,DPE,RC,MU,MUA0)
! Calculate the chemical potential of a bulk LJ fluid at T and RHO
	IMPLICIT NONE

	REAL*8 T,RHOB,DPE,PI,DP,RC
	REAL*8 AN0,AN1,AN2,AN3
	REAL*8 AHS0,AHS1,AHS2,AHS3
	REAL*8 Z(2),AK(2),AA(2),BB(2)
 	REAL*8 AL(2),AS(2),AQ(2),AQD(2)
	REAL*8 A1TEM(2),A2TEM1(2),A2TEM2(2),A2TEMP
	REAL*8 B1TEM(2),B2TEM1(2),B2TEM2(2),B2TEMP
	REAL*8 ABTEM1,ABTEM2
	REAL*8 MUID,MUHS,MUAT,MUCUT,MU,MUA0,MUS
	PARAMETER (PI=3.1415926)
	INTEGER J

! Ideal-gas term
      MUID=DLOG(RHOB)

! Hard-sphere term
	AN0=RHOB
	AN1=0.5*RHOB*DPE
	AN2=PI*RHOB*DPE**2
	AN3=PI*RHOB*DPE**3/6.0
 
 	AHS0=-DLOG(1.-AN3)
	AHS1=AN2/(1.-AN3)
	AHS2=AN1/(1.-AN3)+AN2**2*DLOG(1.-AN3)/(12.*PI*AN3**2)
     &   +AN2**2/(12.*PI*AN3*(1.-AN3)**2)
	AHS3=AN0/(1.-AN3)+AN1*AN2/(1.-AN3)**2
     &   -AN2**3*DLOG(1.-AN3)/(18.*PI*AN3**3)
     &   -AN2**3/(36.*PI*AN3**2*(1.-AN3))
     &   -AN2**3*(1.-3.*AN3)/(36.*PI*AN3**2*(1.-AN3)**3)

      MUHS=AHS0+AHS1*DPE*0.5+AHS2*PI*DPE**2+AHS3*PI*DPE**3/6.0

! van der Waals attraction
	Z(1)=2.9637/DP
 	Z(2)=14.0167/DP

	DO J=1,2
	AK(J)=2.1714*DP*DEXP(Z(J)*(DP-DPE))
	ENDDO

	DO J=1,2
	AL(J)=(1.+0.5*AN3)*Z(J)*DPE+1.+2.*AN3
	AS(J)=(1.-AN3)**2*(Z(J)*DPE)**3+6.*AN3*(1.-AN3)*(Z(J)*DPE)**2
     &     +18.*AN3**2*Z(J)*DPE-12.*AN3*(1.+2.*AN3)
	AQ(J)=(AS(J)+12.*AN3*AL(J)*DEXP(-Z(J)*DPE))
     &     /((1.-AN3)**2*(Z(J)*DPE)**3)
	AQD(J)=(6.*(1.-AN3)*(Z(J)*DPE)**2+36.*AN3*Z(J)*DPE-12.*(1.+5.*AN3)
     &      +12.*((1.+2.*AN3)*Z(J)*DPE+1.+5.*AN3)*DEXP(-Z(J)*DPE))
     &      /((1.-AN3)**3*(Z(J)*DPE)**3)
	ENDDO

	DO J=1,2
	A1TEM(J)=AK(J)*(AL(J)/(Z(J)**2*(1.-AN3)**2*AQ(J))
     &                -(1.+Z(J)*DPE)/Z(J)**2)
	ENDDO

	DO J=1,2
	A2TEM1(J)=AK(J)**2/2./Z(J)/AQ(J)**4
	A2TEM2(J)=AK(J)/DPE/AQ(J)**2
	ENDDO
	A2TEMP=2.*AK(1)*AK(2)/(Z(1)+Z(2))/AQ(1)**2/AQ(2)**2

	DO J=1,2
	B1TEM(J)=AK(J)*(((2.5+0.5*AN3)*Z(J)*DPE+4.+2.*AN3)
     &                /Z(J)**2/(1.-AN3)**3/AQ(J)
     &                -AL(J)*AQD(J)/Z(J)**2/(1.-AN3)**2/AQ(J)**2) 
      ENDDO
	
      DO J=1,2
	B2TEM1(J)=AK(J)**2*AQD(J)/Z(J)/AQ(J)**5
	B2TEM2(J)=AK(J)*AQD(J)/DPE/AQ(J)**3
	ENDDO
	B2TEMP=2.*AK(1)*AK(2)*(AQD(1)*AQ(2)+AQ(1)*AQD(2))
     &      /(Z(1)+Z(2))/AQ(1)**3/AQ(2)**3
										
     	ABTEM1=(DP/DPE)**12/9.-(DP/DPE)**6/3.
 	ABTEM2=(DP/DPE)**12/9.-(DP/DPE)**6/3.+2.*(DP/DPE)**3/9.

	AA(1)=-12.*AN3/(T*DPE**3)*(A1TEM(1)-A1TEM(2))
     &      +48.*AN3/T*ABTEM1
     &      -48.*AN3*(1.+0.5*AN3)/(T*(1.-AN3)**2)*ABTEM2
      AA(2)=-6.*AN3/(T**2*DPE**3)*(A2TEM1(1)+A2TEM1(2)-A2TEMP)
     &      -24.*AN3/T**2*(A2TEM2(1)-A2TEM2(2))*ABTEM2

	BB(1)=AA(1)-12.*AN3**2/(T*DPE**3)*(B1TEM(1)-B1TEM(2))
     &     -24.*AN3**2*(5.+AN3)/(T*(1.-AN3)**3)*ABTEM2
	BB(2)=AA(2)+12.*AN3**2/(T**2*DPE**3)*(B2TEM1(1)+B2TEM1(2)-B2TEMP)
     &     +48.*AN3**2/T**2*(B2TEM2(1)-B2TEM2(2))*ABTEM2
	
	MUAT=AA(1)+AA(2)+BB(1)+BB(2)

! Cut the tail at RC 	
	MUCUT=-16.*PI*RHOB*DP**3*((DP/RC)**9/9.-(DP/RC)**3/3.)/T
	MUA0=MUAT+MUCUT

! Cut and shifted potential (residual)
	MUS=-16.*PI*RHOB*DP**3*(4.*(DP/RC)**9/9.-2.*(DP/RC)**3/3.)/T

! Overall chemical potential
	MU=MUID+MUHS+MUAT+MUS

c	write(*,*) mu,MUID,MUHS,MUA0,MUS
c	write(*,*) mucut,MUID+MUHS+MUA0,MUID+MUHS+MUAT
	RETURN
	END


	SUBROUTINE PRESSURE(T,RHOB,DP,DPE,RC,P,PA0)
! Calculate the pressure of a bulk LJ fluid at T and RHO

	IMPLICIT NONE

	REAL*8 T,RHOB,DP,DPE,RC,PI
	REAL*8 RHOR,ETA

	REAL*8 Z(2),AK(2)
 	REAL*8 AA(2),BB(2)
 	REAL*8 AL(2),AS(2),AQ(2),AQD(2)
	REAL*8 A1TEM(2),A2TEM1(2),A2TEM2(2),A2TEMP
	REAL*8 B1TEM(2),B2TEM1(2),B2TEM2(2),B2TEMP
	REAL*8 ABTEM1,ABTEM2
	integer j
	REAL*8 PCS,PAT,PCUT,P,PA0,PS
	PARAMETER (PI=3.1415926)

	RHOR=RHOB*DP**3
	ETA=RHOB*DPE**3*PI/6.0

      PCS=(1.+ETA+ETA**2-ETA**3)*RHOR*T/(1.-ETA)**3
	 
	Z(1)=2.9637/DP
	Z(2)=14.0167/DP

	DO J=1,2
 	AK(J)=2.1714*DP*DEXP(Z(J)*(DP-DPE))
	ENDDO

	DO J=1,2
	AL(J)=(1.+0.5*ETA)*Z(J)*DPE+1.+2.*ETA
	AS(J)=(1.-ETA)**2*(Z(J)*DPE)**3+6.*ETA*(1.-ETA)*(Z(J)*DPE)**2
     &     +18.*ETA**2*Z(J)*DPE-12.*ETA*(1.+2.*ETA)
	AQ(J)=(AS(J)+12.*ETA*AL(J)*DEXP(-Z(J)*DPE))
     &     /((1.-ETA)**2*(Z(J)*DPE)**3)
	AQD(J)=(6.*(1.-ETA)*(Z(J)*DPE)**2+36.*ETA*Z(J)*DPE-12.*(1.+5.*ETA)
     &      +12.*((1.+2.*ETA)*Z(J)*DPE+1.+5.*ETA)*DEXP(-Z(J)*DPE))
     &      /((1.-ETA)**3*(Z(J)*DPE)**3)
	ENDDO

	DO J=1,2
	A1TEM(J)=AK(J)*(AL(J)/(Z(J)**2*(1.-ETA)**2*AQ(J))
     &                -(1.+Z(J)*DPE)/Z(J)**2)
	ENDDO

	DO J=1,2
	A2TEM1(J)=AK(J)**2/2./Z(J)/AQ(J)**4
	A2TEM2(J)=AK(J)/DPE/AQ(J)**2
	ENDDO
	A2TEMP=2.*AK(1)*AK(2)/(Z(1)+Z(2))/AQ(1)**2/AQ(2)**2

	DO J=1,2
	B1TEM(J)=AK(J)*(((2.5+0.5*ETA)*Z(J)*DPE+4.+2.*ETA)
     &                /Z(J)**2/(1.-ETA)**3/AQ(J)
     &                -AL(J)*AQD(J)/Z(J)**2/(1.-ETA)**2/AQ(J)**2) 
      ENDDO
	
      DO J=1,2
	B2TEM1(J)=AK(J)**2*AQD(J)/Z(J)/AQ(J)**5
	B2TEM2(J)=AK(J)*AQD(J)/DPE/AQ(J)**3
	ENDDO
	B2TEMP=2.*AK(1)*AK(2)*(AQD(1)*AQ(2)+AQ(1)*AQD(2))
     &      /(Z(1)+Z(2))/AQ(1)**3/AQ(2)**3
										
     	ABTEM1=(DP/DPE)**12/9.-(DP/DPE)**6/3.
 	ABTEM2=(DP/DPE)**12/9.-(DP/DPE)**6/3.+2.*(DP/DPE)**3/9.

	AA(1)=-12.*ETA/(T*DPE**3)*(A1TEM(1)-A1TEM(2))
     &      +48.*ETA/T*ABTEM1
     &      -48.*ETA*(1.+0.5*ETA)/(T*(1.-ETA)**2)*ABTEM2
      AA(2)=-6.*ETA/(T**2*DPE**3)*(A2TEM1(1)+A2TEM1(2)-A2TEMP)
     &      -24.*ETA/T**2*(A2TEM2(1)-A2TEM2(2))*ABTEM2

	BB(1)=AA(1)-12.*ETA**2/(T*DPE**3)*(B1TEM(1)-B1TEM(2))
     &     -24.*ETA**2*(5.+ETA)/(T*(1.-ETA)**3)*ABTEM2
	BB(2)=AA(2)+12.*ETA**2/(T**2*DPE**3)*(B2TEM1(1)+B2TEM1(2)-B2TEMP)
     &     +48.*ETA**2/T**2*(B2TEM2(1)-B2TEM2(2))*ABTEM2

	PAT=(BB(1)+BB(2))*RHOR*T

! Cut the tail at RC 	
	PCUT=-8.*PI*RHOR**2*((DP/RC)**9/9.-(DP/RC)**3/3.)
	PA0=PAT+PCUT

! Cut-and-shifted potential
	PS=-8.*PI*RHOR**2*(4./9.*(DP/RC)**9-2./3.*(DP/RC)**3)

! Overall
	P=PCS+PAT+PS

c	write(*,*) P, RHOR*T,PCS-RHOR*T,PA0,PS
c	write(*,*) Pcut,PCS+PA0, PCS+PAT+PS, PCS+PAT

	RETURN
	END
