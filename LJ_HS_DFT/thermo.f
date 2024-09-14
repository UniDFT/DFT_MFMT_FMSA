	SUBROUTINE THERMO(IR,NR,RHO,GAMMA,RGP)
c calculate the reduced EXCESS grand potential (RGP) 
C   and the reduced EXCESS Helmholtz energy functional (RF)

	IMPLICIT NONE

      REAL*8 RHOB,DC,DP,T,MIX_F,RB,RC,DR,TOLE,PI
      COMMON RHOB,DC,DP,T,MIX_F,RB,RC,DR,TOLE,PI

	REAL*8 DPE,Z(2),X(2,7)
	COMMON /FMSA/ DPE,Z,X

      INTEGER NP
      PARAMETER(NP=2**13)

	INTEGER I,K,KIN,IR,NR
	REAL*8 DRHO(NP),RHO(NP)
	REAL*8 N(3,NP),DDR(3,NP),DFATT
	REAL*8 N2,N3,NV2,DFS, RFEXS

	REAL*8 RFRD,DDA,MUB,MUA0,FAB,PA0,P
	REAL*8 RFID,RFEXR,RFEXA,RF,RGP,RG0

	REAL*8 R,GAMMA,VACC,RMAX,RMIN 	
	REAL*8 FSB, MuSB

 
! Integration boundary
	RMIN=(DC+DP)/2.
	RMAX=DC/2+FLOAT(NR-1)*DR

!	system volume
      VACC=4.*PI*(RMAX**3-RMIN**3)/3.

! THE FIRST NON-ZERO DENSITY POINT (FROM THE PARTICLE SURFACE)
	KIN=NINT(DP/2./DR)+1  

c Number of molecules
	GAMMA=0.D0

	DO 10 K=KIN,NR
	I=K-IR+1
	R=DC/2+FLOAT(K-1)*DR

	IF ((K.EQ.KIN).OR.(K.EQ.NR)) THEN
	GAMMA=GAMMA+2.*PI*R**2*RHO(I)*DR

	ELSE

	GAMMA=GAMMA+4.*PI*R**2*RHO(I)*DR

	ENDIF

10	CONTINUE


C Reduced excess Helmholtz free energy  (/kBT)

! ideal-gas term
	RFID=0.0
	DO K=KIN,NR
	I=K-IR+1
	R=DC/2+FLOAT(K-1)*DR
	IF ((K.EQ.KIN).OR.(K.EQ.NR)) THEN
	RFID=RFID+2.*PI*R**2*RHO(I)*(DLOG(RHO(I))-1.)*DR
	ELSE
	RFID=RFID+4.*PI*R**2*RHO(I)*(DLOG(RHO(I))-1.)*DR
	ENDIF
	ENDDO

! hard-sphere term
	RFEXR=0.0

! CALCULATE THE WEIGHTED DENSITIES
	CALL DREP(IR,NR,DPE,RHO,N,DDR)		

	DO K=KIN,NR
	I=K-IR+1
	R=DC/2+FLOAT(K-1)*DR
	N2=N(1,I)
	N3=N(2,I)
	NV2=N(3,I)

	IF (N3.GT.1.E-8) THEN
	RFRD=-N2*DLOG(1.-N3)/(PI*DPE**2)
     &     +(N2**2-NV2**2)/(1.-N3)/(2.*PI*DPE)
     &     +(DLOG(1.-N3)+N3/(1.-N3)**2)
     &     *(N2**3-3.*N2*NV2**2)/(36.*PI*N3**2)
	ELSE
	RFRD=0.D0
	ENDIF

	IF ((K.EQ.KIN).OR.(K.EQ.NR)) THEN
	RFEXR=RFEXR+2.*PI*R**2*RFRD*DR
	ELSE
	RFEXR=RFEXR+4.*PI*R**2*RFRD*DR
	ENDIF
	ENDDO

! van der Waals attraction
	DO I=1,NP
	DRHO(I)=RHO(I)-RHOB
	ENDDO

!  bulk chemical potential and pressure
	CALL CP(T,RHOB,DP,DPE,RC,MUB,MUA0)
	CALL PRESSURE(T,RHOB,DP,DPE,RC,P,PA0)

!  bulk Helmholtz energy density (attractive part, 0 to rc) 
	FAB=MUA0*RHOB-PA0/T

	RFEXA=0.0
 	DO K=KIN,NR
	I=K-IR+1
	R=DC/2+FLOAT(K-1)*DR
	CALL DFEXATT(K,IR,NR,DRHO,DFATT)
	DDA=MUA0+0.5*DFATT
	IF ((K.EQ.KIN).OR.(K.EQ.NR)) THEN
	RFEXA=RFEXA+2.*PI*R**2*(FAB+DRHO(I)*DDA)*DR
	ELSE
	RFEXA=RFEXA+4.*PI*R**2*(FAB+DRHO(I)*DDA)*DR
	ENDIF
	ENDDO

! Excess Helmholtz energy due to the shifted potential 
	RFEXS=0.0

! bulk shifted free energy density and chemical potential
	FSB=-8.*PI*RHOB**2*DP**3*((DP/RC)**9-(DP/RC)**3)/(3.*T)
	MuSB=-16.*PI*RHOB*DP**3*((DP/RC)**9-(DP/RC)**3)/(3.*T)

	DO K=KIN,NR
        I=K-IR+1
	R=DC/2+FLOAT(K-1)*DR
	CALL DFEXS(K,IR,NR,RHO,DFS)
	DDA=(MUSB+DFS)*0.5
	IF ((K.EQ.KIN).OR.(K.EQ.NR)) THEN
	RFEXS=RFEXS+2.*PI*R**2*(FSB+DDA*DRHO(I))*DR
	ELSE
	RFEXS=RFEXS+4.*PI*R**2*(FSB+DDA*DRHO(I))*DR
	ENDIF
	ENDDO

! overall Helmgolz energy/KBT
	RF=RFID+RFEXR+RFEXA+RFEXS
		
!  reduced one-body potential 
	RG0=0.0

	DO K=KIN,NR
	I=K-IR+1
	R=DC/2+FLOAT(K-1)*DR
	IF ((K.EQ.KIN).OR.(K.EQ.NR)) THEN
	RG0=RG0+2.*PI*R**2*RHO(I)*MUB*DR
	ELSE
	RG0=RG0+4.*PI*R**2*RHO(I)*MUB*DR
	ENDIF
	ENDDO

! reduced grand potential/kBT
	RGP=RF-RG0

c	write(*,*) 'grand potential density=',rgp/VACC
c	write(*,*) 'helmholtz energy density', RF/VACC
c	write(*,*) 'gibbs energy density', RG0/VACC
c	write(*,*) 'bulk pressure=', P/T
c	write(*,*) 'chemical potential=', MUB,RG0/gamma
c	write(*,*) 'average density', gamma/VACC
	RETURN
	END
