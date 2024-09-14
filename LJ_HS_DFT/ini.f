	  SUBROUTINE READDAT (ICODE) 
	  IMPLICIT NONE
	  INTEGER ICODE
	   
        REAL*8 RHOB,DC,DP,T,MIX_F,RB,RC,DR,TOLE,PI
        COMMON RHOB,DC,DP,T,MIX_F,RB,RC,DR,TOLE,PI
	   
	  OPEN (66, FILE='input')
! ICODE=0 OR 1; 1 FOR RESUMING CALCULATION 
	  READ(66,*) 
	  READ(66,*) ICODE

! BULK DENSITY OF  THE LENNARD-JONES FLUID 
	  READ(66,*) 
	  READ(66,*) RHOB

! REDUCED TEMPERATURE BETWEEN LENNARD-JONES SPHERES              
	  READ(66,*) 
	  READ(66,*) T

! DIAMETER OF THE FIXED COLLOID PARTICLE            
	  READ(66,*) 
	  READ(66,*) DC

! DIAMETER OF THE LENNARD-JONES SPHERES              
	  READ(66,*) 
        READ(66,*) DP

! CUTOFF DISTANCE OF THE LENNARD-JONES POTENTIAL
	  READ(66,*) 
        READ(66,*) RC

! UPPER COMPUTATION BOUNDARY                  
	  READ(66,*) 
        READ(66,*) RB

! STEP LENGTH             
	  READ(66,*) 
        READ(66,*) DR

! MIXING PARAMETER FOR PICARD ITERATION       
	  READ(66,*) 
        READ(66,*) MIX_F

! TOLERANCE
	  READ(66,*) 
        READ(66,*) TOLE           
	
! PI              
	  READ(66,*) 
        READ(66,*) PI           

	  CLOSE(66)

	return
	END

