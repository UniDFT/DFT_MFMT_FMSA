	SUBROUTINE DFEXATT(I,IR,NR,RHO,DFATT)

!     THIS SUBROUTINE IS USED TO CALCULATE THE FUNCTIONAL DERIVITIVE OF HELMHOLTZ
!     ENERGY DUE TO ATTRACTION

	IMPLICIT NONE
      REAL*8 RHOB,DC,DP,T,MIX_F,RB,RC,DR,TOLE,PI
      COMMON RHOB,DC,DP,T,MIX_F,RB,RC,DR,TOLE,PI

      INTEGER NP
      PARAMETER(NP=2**13)

      REAL*8 RHO(NP)
	REAL*8 DPE,DFATT

	REAL*8 R,R1,lambda
	INTEGER KIN,KIP,I,K,J
	INTEGER IR,NR,IP

	REAL*8 Z(2),X(2,7)
	REAL*8 TEMP,TEMP1(3)

	COMMON /FMSA/DPE,Z,X

	R=DC/2+FLOAT(I-1)*DR

      DFATT=0.0

! INTEGRATION BOUNDARIES	 
      KIN=IR+NINT(DP/DR)		  !STARTING POINT OF NONE ZERO DENSITY
      KIP=NR+NINT(RC/DR)        !RANGE OF CALCUATION
			
	DO 10 K=KIN,KIP
	IP=K-IR+1	 !INDEX FOR THE DENSITY PROFILE
	R1=DC/2+FLOAT(K-1)*DR

	lambda=ABS(R-R1)/DPE

! SEVEN CASES OF INTEGRATION BOUNDARUES

	IF ((R+R1).GT.RC) THEN	 

! CASE I:|R-R1|<DPE, AND R+R1>RC
		IF (ABS(R-R1).LT.DPE) THEN

		DO J=1,2
		TEMP1(J)=-(X(J,1)+X(J,2))*(1.-EXP(Z(J)*(1.-lambda)))/Z(J)
     &         +X(J,3)*(1.-EXP(Z(J)*(lambda-1.0)))/Z(J)
     &         +X(J,4)*(1.-lambda**5)/5.0
     &         +X(J,5)*(1.-lambda**3)/3.0
     &         +X(J,6)*(1.-lambda**2)/2.0
     &         +X(J,7)*(1.-lambda)
		ENDDO

		TEMP1(3)=(0.4*DP**12/RC**10-DP**6/RC**4+0.6*DP**2)/T

		TEMP=TEMP1(1)-TEMP1(2)+TEMP1(3)

! CASE II:DPE<|R-R1|<DP, AND R+R1>RC
		ELSEIF ((ABS(R-R1).GE.DPE).AND.(ABS(R-R1).LT.DP)) THEN

		TEMP=(0.4*DP**12/RC**10-DP**6/RC**4+0.6*DP**2)/T

! CASE III: |R-R1|>DP, AND R+R1>RC
		ELSEIF ((ABS(R-R1).GE.DP).AND.(ABS(R-R1).LE.RC)) then

		TEMP=(0.4*(DP**12/RC**10-DP**12/(R-R1)**10)
     &     -(DP**6/RC**4-DP**6/(R-R1)**4))/T

! CASE IV: |R-R1|>RC, C(R)=0
 		ELSEIF (ABS(R-R1).GT.RC) THEN	 
		TEMP=0.0

		ENDIF  ! END OF R+R1>RC

	ELSE 
! CASE V: |R-R1|<DPE, AND R+R1<RC
		IF (ABS(R-R1).LT.DPE) THEN

		DO J=1,2
 		TEMP1(J)=-(X(J,1)+X(J,2))*(1.-EXP(Z(J)*(1.-lambda)))/Z(J)
     &         +X(J,3)*(1.-EXP(Z(J)*(lambda-1.0)))/Z(J)
     &         +X(J,4)*(1.-lambda**5)/5.0
     &         +X(J,5)*(1.-lambda**3)/3.0
     &         +X(J,6)*(1.-lambda**2)/2.0
     &         +X(J,7)*(1.-lambda)
		ENDDO

		TEMP1(3)=(0.4*DP**12/(R+R1)**10-DP**6/(R+R1)**4+0.6*DP**2)/T

		TEMP=TEMP1(1)-TEMP1(2)+TEMP1(3)

! CASE VI: DPE<|R-R1|<DP, AND R+R1<RC
		ELSEIF ((ABS(R-R1).GE.DPE).AND.(ABS(R-R1).LT.DP)) THEN
		TEMP=(0.4*DP**12/(R+R1)**10-DP**6/(R+R1)**4+0.6*DP**2)/T

! CASE VII: |R-R1|>DP, AND R+R1<RC
		ELSE
		TEMP=(0.4*(DP**12/(R+R1)**10-DP**12/(R-R1)**10)
     &     -(DP**6/(R+R1)**4-DP**6/(R-R1)**4))/T

		ENDIF
      ENDIF

	DFATT=DFATT+RHO(IP)*R1*TEMP

10    CONTINUE

	DFATT=-DFATT*DR*2.*PI/R

      RETURN 
	END
     

	SUBROUTINE CPM(AMU)

!     THIS SUBROUTINE IS USED TO CALCULATE THE BULK CHEMICAL POTENTIAL
!      FOR ID AND HS TERMS PLUS A RESIDUE TERM DUE TO THE ATTRACTION 
      REAL*8 RHOB,DC,DP,T,MIX_F,RB,RC,DR,TOLE,PI
      COMMON RHOB,DC,DP,T,MIX_F,RB,RC,DR,TOLE,PI

	REAL*8 DPE,ETA
	REAL*8 CK0(2),Z(2),X(2,7)
	REAL*8 MUID,MUHS,MUAT,AMU,MUCUT,MUS

	COMMON /FMSA/DPE,Z,X

	ETA=RHOB*PI*DPE**3/6.0

	MUID=DLOG(RHOB)

	MUHS=ETA*(8.-9.*ETA+3.*ETA**2)/(1.-ETA)**3

	DO I=1,2
	CK0(I)=(X(I,1)+X(I,2))*(EXP(Z(I))-Z(I)-1.0)/Z(I)**2
     &      +X(I,3)*(EXP(-Z(I))+Z(I)-1.0)/Z(I)**2
     &      +X(I,4)/6.+X(I,5)/4.+X(I,6)/3.+X(I,7)/2.
	ENDDO

	MUAT=-4.*PI*RHOB*DPE*(CK0(1)-CK0(2))

! CONTRIBUTION TO MU FROM 1 TO RC	
	MUCUT=-16.*PI*RHOB*(DP**12/RC**9-3.*DP**6/RC**3
     &      +2.*DP**3)/(9.*T)

! Shifted potential in the bulk
	MUS=-16.*PI*RHOB*DP**3*((DP/RC)**9-(DP/RC)**3)/(3.*T)

	AMU=MUID+MUHS+MUAT+MUCUT+MUS

c	write(*,*) RHOB, MUHS, MUAT+MUCUT,MUS,'bulk'

	RETURN
	END


	SUBROUTINE CRATT

!     THIS SUBROUTINE IS USED TO CALCULATE THE PARAMETERS FOR 
!     DIRECT CORRELATION FUNCTION DUE TO ATTRACTION

        REAL*8 RHOB,DC,DP,T,MIX_F,RB,RC,DR,TOLE,PI
        COMMON RHOB,DC,DP,T,MIX_F,RB,RC,DR,TOLE,PI

	REAL*8 DPE
	REAL*8 Z(2),X(2,7)

	REAL*8 AK0,Z0(2)
	REAL*8 BETA(2)
	REAL*8 ETA

	REAL*8 AQ,AL,AS

	REAL*8 TEMP(2)

	COMMON /FMSA/DPE,Z,X

! EFFECTIVE HS SPHERE DIAMETER FOR LJ PARTICLES
	DPE=DP*(1.+0.2977*T)/(1.+0.33163*T+0.00104771*T**2)

  	AK0=2.1714
      Z0(1)=2.9637
	Z0(2)=14.0167
      ETA=RHOB*DPE**3*PI/6.0

	DO I=1,2
	Z(I)=Z0(I)*DPE/DP
	BETA(I)=AK0*EXP(Z0(I)*(1.-DPE/DP))/(DPE/DP)/T
	ENDDO

	DO I=1,2
	AL=(1.+0.5*ETA)*Z(I)+1.+2.*ETA
	AS=(1.-ETA)**2*Z(I)**3+6.*ETA*(1.-ETA)*Z(I)**2
     &   +18.*ETA**2*Z(I)-12.*ETA*(1.+2.*ETA)
	AQ=(AS+12.*ETA*AL*EXP(-Z(I)))/((1.-ETA)**2*Z(I)**3)

	TEMP(I)=-DPE**2*BETA(I)/((1.-ETA)**4*Z(I)**6*AQ**2)

	X(I,1)=DPE**2*BETA(I)
	X(I,2)=TEMP(I)*AS**2
	X(I,3)=TEMP(I)*144.*ETA**2*AL**2
	X(I,4)=-TEMP(I)*12.*ETA**2
     &       *((1.+2.*ETA)**2*Z(I)**4+(1.-ETA)*(1.+2.*ETA)*Z(I)**5)
	X(I,5)=TEMP(I)*12.*ETA
     &       *(AS*AL*Z(I)**2-(1.-ETA)**2*(1.+0.5*ETA)*Z(I)**6)
	X(I,6)=-TEMP(I)*24.*ETA
     &       *((1.+2.*ETA)**2*Z(I)**4
     &       +(1.-ETA)*(1.+2.*ETA)*Z(I)**5)
	X(I,7)=TEMP(I)*24.*ETA*AS*AL

	ENDDO

	RETURN
	END

