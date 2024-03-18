! ======================================================================
! User Subroutine UEL for Abaqus: two and three dimension for phase
! and displacement problem:
! Type 1: C2D6  displacement triangle element
! Type 2: C2D8  displacement rectangular element
! Type 3: C3D12 displacement tetrahedral element
! Type 4: C3D24 displacement brick element
!
! Type 5: C2P3 phase-field triangle element
! Type 6: C2P4 phase-field rectangular element
! Type 7: C3P4 phase-field tetrahedral element
! Type 8: C3P8 phase-field brick element
! ======================================================================
! Material properties to be given through the input file (*.inp):
!
! For Type 1 element (stress-strain):
! PROPS(1) = Young's modulus (E)
! PROPS(2) = Poisson's ratio (nu)
! PROPS(3) = AT1 (xi = 1) and AT2 (xi = 0)
! PROPS(4) = Length scale parameter (lc)
! PROPS(5) = Crack surface energy (gc)
!
! For Type 2 element (phase field):
! PROPS(1) = Young's modulus (E)
! PROPS(2) = AT1 (xi = 1) and AT2 (xi = 0)
! PROPS(3) = Length scale parameter (lc)
! PROPS(4) = Crack surface energy (gc)
!
!     ==================================================================
!     Comments on solution dependent variables
!     ==================================================================
!
!     Stress/strain element:
!                           SDV(1) - phase field
!                           SDV(2) - Elastic strain energy
!                           SDV(3) - von Mises stress
!
!     Phase-field element:
!                         SDV(1) - phase field
!                         SDV(2) - history energy
!
! ======================================================================
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     4     PERIOD)
!     ==================================================================
      INCLUDE 'ABA_PARAM.INC'
!     ==================================================================
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,TWO=2.D0,THREE=3.D0,
     1 TOLER=1.0D-8,FOUR=4.D0,RP75=0.375D0,HALF=0.5D0,FIVE=5.D0,
     2 SIX=6.D0,SEVEN=7.D0,EIGHT=8.D0,TEN=10.D0,DTMIN=1.0D-5,
     3 DPHMIN=1.0D-2,NSTVTO=2,NSTVTT=3,NSTV=9,N_ELEM=82406)
!     ==================================================================
!     Initialization for all the element types
!     ==================================================================
      DIMENSION RHS(MLVARX,1),AMATRX(NDOFEL,NDOFEL),
     1     SVARS(NSVARS),ENERGY(8),PROPS(NPROPS),COORDS(MCRD,NNODE),
     2     U(NDOFEL),DU(MLVARX,1),V(NDOFEL),A(NDOFEL),TIME(2),
     3     PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     4     DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
     5     JPROPS(*)
     
       INTEGER I,J,L,K,K1,K2,K3,K4,IX,IY,NALL

       REAL*8 AINTW(8),XII(8,MCRD),XI(MCRD),dNdxi(NNODE,MCRD),
     1 VJACOB(3,3),dNdx(NNODE,MCRD),VJABOBINV(3,3),BP(MCRD,NDOFEL),
     2 AN(8),DP(MCRD),SDV(NSTV),BB(6,NDOFEL),CMAT(6,6),EPS(6),STRESS(6),
     3 VNI(MCRD,NDOFEL),ULOC(MCRD),PHASENOD(NNODE),PNEWDTIP(8),
     4 EIGV(3),ALPHAI(3),VECTI(3),EPSC(6)

       REAL*8 DTM,HIST,CLPAR,GCPAR,EMOD,ENU,PARK,ENG,ENGDG,ENGD,
     1 REITER,PHASE0,DPHASE0,ELAMEG,ELAMEL,THCK,ALPHA,DALPHA,DDALPHA,
     2 FT,C0,OMEGA,DOMEGA,DDOMEGA,PHASE_SOURCE,DPHASE_SOURCE,NN,NP,
     3 ENGP,ATPAR,ENGMAX
	  
!
       COMMON/KUSER/USRVAR(N_ELEM,NSTV,8)
!
!     History variables
        ENG=ZERO
        ENGDG=ZERO
        NALL=NSTVTO+NSTVTT           				
!     ==================================================================
!     ******************************************************************
!     Constructing elemet TYPE 1 and 2 and 3 and 4 (stress/strain-displacement)
!     ******************************************************************
!     ==================================================================
!      
       IF ((JTYPE.EQ.ONE).OR.(JTYPE.EQ.TWO).OR.
     1    (JTYPE.EQ.THREE).OR.(JTYPE.EQ.FOUR)) THEN
!     ==================================================================
!     Time an iteration variables
!     ==================================================================
!
       IF (TIME(2).EQ.ZERO) THEN
        TIMEZ=-999.D0
        TIMEZLOC=-99.D0
        DO K2=1,NSTV
         DO K3=1,8
           USRVAR(JELEM,K2,K3)=ZERO
         END DO
        END DO        
       ELSE
        TIMEZ=USRVAR(JELEM,NALL+1,1)
        TIMEZLOC=TIME(2)-DTIME
       ENDIF
       DTZERO=USRVAR(JELEM,NALL+2,1)
       IF (TIMEZ.LT.TIMEZLOC) THEN
        USRVAR(JELEM,NALL+1,1)=TIMEZLOC
        USRVAR(JELEM,NALL+2,1)=DTIME
        USRVAR(JELEM,NALL+3,1)=ZERO
        USRVAR(JELEM,NALL+4,1)=ZERO
       ELSE
        IF (DTZERO.GT.DTIME*(ONE+TOLER)) THEN
         USRVAR(JELEM,NALL+2,1)=DTIME
         USRVAR(JELEM,NALL+3,1)=USRVAR(JELEM,NALL+3,1)+ONE
         USRVAR(JELEM,NALL+4,1)=ZERO
        ELSE
         USRVAR(JELEM,NALL+4,1)=USRVAR(JELEM,NALL+4,1)+ONE
        ENDIF
       ENDIF      
       REITER=USRVAR(JELEM,NALL+3,1)
       STEPITER=USRVAR(JELEM,NALL+4,1)       
!     ==================================================================
!     Material parameters
!     ==================================================================
       EMOD =PROPS(1)
       ENU  =PROPS(2)
       ATPAR=PROPS(3)
       CLPAR=PROPS(4)
       GCPAR=PROPS(5)
       PARK =1.D-6 	   
       THCK =ONE
       ELAMEG=EMOD/(TWO*(ONE+ENU))
       ELAMEL=ELAMEG*TWO*ENU/(ONE-TWO*ENU) 
!     ==================================================================
!     Crack geometry function 
!     ==================================================================	   
!        	  
!      The critical strength upon crack nucleation (not
!      necessarily the peak strength	  
       IF (NINT(ATPAR)==0) THEN ! AT2
         FT = ZERO
         C0 = HALF
       ELSEIF (NINT(ATPAR)==1) THEN ! AT1
         FT = SQRT(RP75*EMOD*GCPAR/CLPAR)
         C0 = TWO/THREE
       ELSE
         WRITE (*,*) '**ERROR: MODEL LAW NO. ', ATPAR, 
     1                 'DOES NOT EXIST!'
         CALL XIT
       ENDIF 	   	   
!	   
!     ==================================================================
!     Initial preparations
!     ==================================================================
       DO K1 = 1, NDOFEL                      
        DO KRHS = 1, NRHS
         RHS(K1,KRHS) = ZERO
        END DO
        DO K2 = 1, NDOFEL
         AMATRX(K2,K1) = ZERO
        END DO
       END DO
!     ==================================================================
!     Local coordinates and weights
!     ==================================================================
       IF (JTYPE.EQ.ONE) THEN
        XII(1,1) = ONE/THREE
        XII(1,2) = ONE/THREE
        INNODE = 1.D0
        AINTW(1) = HALF
       ELSEIF (JTYPE.EQ.TWO) THEN
        XII(1,1) = -ONE/THREE**HALF
        XII(1,2) = -ONE/THREE**HALF
        XII(2,1) = ONE/THREE**HALF
        XII(2,2) = -ONE/THREE**HALF
        XII(3,1) = ONE/THREE**HALF
        XII(3,2) = ONE/THREE**HALF
        XII(4,1) = -ONE/THREE**HALF
        XII(4,2) = ONE/THREE**HALF
        INNODE = 4.D0
        DO I=1,INNODE
         AINTW(I) = ONE
        END DO	   
       ELSEIF (JTYPE.EQ.THREE) THEN
        XII(1,1) = ONE/FOUR
        XII(1,2) = ONE/FOUR
        XII(1,3) = ONE/FOUR
        INNODE=1.D0
        DO I=1,INNODE
         AINTW(I) = ONE/SIX
        END DO
       ELSEIF (JTYPE.EQ.FOUR) THEN
        XII(1,1) = MONE/THREE**HALF
        XII(1,2) = MONE/THREE**HALF
        XII(1,3) = MONE/THREE**HALF
        XII(2,1) = ONE/THREE**HALF
        XII(2,2) = MONE/THREE**HALF
        XII(2,3) = MONE/THREE**HALF
        XII(3,1) = ONE/THREE**HALF
        XII(3,2) = ONE/THREE**HALF
        XII(3,3) = MONE/THREE**HALF
        XII(4,1) = MONE/THREE**HALF
        XII(4,2) = ONE/THREE**HALF
        XII(4,3) = MONE/THREE**HALF
        XII(5,1) = MONE/THREE**HALF
        XII(5,2) = MONE/THREE**HALF
        XII(5,3) = ONE/THREE**HALF
        XII(6,1) = ONE/THREE**HALF
        XII(6,2) = MONE/THREE**HALF
        XII(6,3) = ONE/THREE**HALF
        XII(7,1) = ONE/THREE**HALF
        XII(7,2) = ONE/THREE**HALF
        XII(7,3) = ONE/THREE**HALF
        XII(8,1) = MONE/THREE**HALF
        XII(8,2) = ONE/THREE**HALF
        XII(8,3) = ONE/THREE**HALF
        INNODE = 8.D0
        DO I=1,INNODE
         AINTW(I) = ONE
        END DO
       ENDIF
!
!     ==================================================================
!     Determining ration of new time increment
!     ==================================================================

        DO K1=1,8
         PNEWDTIP(K1)=TWO
        END DO
!
!     ==================================================================
!     Calculating properties at each integration point
!     ==================================================================
       DO INPT=1,INNODE
!     Initializing solution dependent variables (phase,history)
        DO I=1,NSTVTT
          SDV(I)=SVARS(NSTVTT*(INPT-1)+I)
        END DO
!
!     Local coordinates of the integration point
        DO I=1,MCRD
         XI(I) = XII(INPT,I)
        END DO
!     Shape functions and local derivatives
        IF (JTYPE.EQ.ONE) THEN
         CALL SHAPEFUNTRI(AN,dNdxi,XI)
        ELSEIF (JTYPE.EQ.TWO) THEN
         CALL SHAPEFUNQUAD(AN,dNdxi,XI)
        ELSEIF (JTYPE.EQ.THREE) THEN
         CALL SHAPEFUNTET(AN,dNdxi,XI)
        ELSEIF (JTYPE.EQ.FOUR) THEN
         CALL SHAPEFUNBRICK(AN,dNdxi,XI)
        ENDIF
!     Shape functions		
        IF (MCRD.EQ.2) THEN
! --------- 2D Elements ------------ 		
         IY=ZERO
         DO I = 1,NNODE
          IX=IY+1
          IY=IX+1
          VNI(1,IX)=AN(I)
          VNI(1,IY)=ZERO
          VNI(2,IX)=ZERO
          VNI(2,IY)=AN(I)
         END DO
        ELSE
! --------- 3D Elements ------------  
         IZ=ZERO
         DO I = 1,NNODE
          IX=IZ+ONE
          IY=IX+ONE
          IZ=IY+ONE
          VNI(1,IX)=AN(I)
          VNI(2,IX)=ZERO
          VNI(3,IX)=ZERO
          VNI(1,IY)=ZERO
          VNI(2,IY)=AN(I)
          VNI(3,IY)=ZERO
          VNI(1,IZ)=ZERO
          VNI(2,IZ)=ZERO
          VNI(3,IZ)=AN(I)
         END DO
        ENDIF		
		
!     Jacobian
        DO I = 1,MCRD
         DO J = 1,MCRD
          VJACOB(I,J) = ZERO
          DO K = 1,NNODE
           VJACOB(I,J) = VJACOB(I,J) + COORDS(I,K)*dNdxi(K,J)
          END DO
         END DO
        END DO
!        
        DTM = ZERO
        IF (MCRD.EQ.2) THEN
! --------- 2D Elements ------------        
         DTM = VJACOB(1,1)*VJACOB(2,2)-VJACOB(1,2)*VJACOB(2,1)
        ELSE
! --------- 3D Elements ------------        
         DTM = VJACOB(1,1)*VJACOB(2,2)*VJACOB(3,3)+VJACOB(1,2)*
     1   VJACOB(2,3)*VJACOB(3,1)+VJACOB(1,3)*VJACOB(2,1)*
     2   VJACOB(3,2)-VJACOB(3,1)*VJACOB(2,2)*VJACOB(1,3)-
     3   VJACOB(3,2)*VJACOB(2,3)*VJACOB(1,1)-VJACOB(3,3)*
     4   VJACOB(2,1)*VJACOB(1,2)        
        ENDIF	 	 	 	 
	 
        IF (DTM.LT.ZERO) THEN
         WRITE(7,*) 'Negative Jacobian',DTM
         CALL XIT
        ENDIF
		
!     Inverse of Jacobian
        IF (MCRD.EQ.2) THEN
! --------- 2D Elements ------------        
         VJABOBINV(1,1)= VJACOB(2,2)/DTM
         VJABOBINV(1,2)=-VJACOB(1,2)/DTM
         VJABOBINV(2,1)=-VJACOB(2,1)/DTM
         VJABOBINV(2,2)= VJACOB(1,1)/DTM
        ELSE
! --------- 3D Elements ------------        
         VJABOBINV(1,1)=(VJACOB(2,2)*VJACOB(3,3)-VJACOB(2,3)*
     1    VJACOB(3,2))/DTM
         VJABOBINV(1,2)=-(VJACOB(1,2)*VJACOB(3,3)-VJACOB(3,2)*
     1    VJACOB(1,3))/DTM
         VJABOBINV(1,3)=(VJACOB(1,2)*VJACOB(2,3)-VJACOB(1,3)*
     1    VJACOB(2,2))/DTM
         VJABOBINV(2,1)=-(VJACOB(2,1)*VJACOB(3,3)-VJACOB(2,3)*
     1    VJACOB(3,1))/DTM
         VJABOBINV(2,2)=(VJACOB(1,1)*VJACOB(3,3)-VJACOB(1,3)*
     1    VJACOB(3,1))/DTM
         VJABOBINV(2,3)=-(VJACOB(1,1)*VJACOB(2,3)-VJACOB(1,3)*
     1    VJACOB(2,1))/DTM
         VJABOBINV(3,1)=(VJACOB(2,1)*VJACOB(3,2)-VJACOB(2,2)*
     1    VJACOB(3,1))/DTM
         VJABOBINV(3,2)=-(VJACOB(1,1)*VJACOB(3,2)-VJACOB(1,2)*
     1    VJACOB(3,1))/DTM
         VJABOBINV(3,3)=(VJACOB(1,1)*VJACOB(2,2)-VJACOB(1,2)*
     1    VJACOB(2,1))/DTM      
        ENDIF
!        
!     Derivatives of shape functions respect to global ccordinates
        DO K = 1,NNODE
         DO I = 1,MCRD
          dNdx(K,I) = ZERO
          DO J = 1,MCRD
           dNdx(K,I) = dNdx(K,I) + dNdxi(K,J)*VJABOBINV(J,I)
          END DO
         END DO
        END DO
!
!     Calculating B matrix (B=LN)	   
        IF (MCRD.EQ.2) THEN
! --------- 2D Elements ------------ 
         IY=0
         DO INODE=1,NNODE
          IX=IY+1
          IY=IX+1
          BB(1,IX)= dNdx(INODE,1)
          BB(1,IY)= ZERO
          BB(2,IX)= ZERO
          BB(2,IY)= dNdx(INODE,2)
          BB(3,IX)= dNdx(INODE,2)
          BB(3,IY)= dNdx(INODE,1)
         END DO
!
        ELSE
! --------- 3D Elements ------------ 
         IZ=ZERO
         DO INODE=1,NNODE
          IX=IZ+ONE
          IY=IX+ONE
          IZ=IY+ONE
          BB(1,IX)= dNdx(INODE,1)
          BB(2,IX)= ZERO
          BB(3,IX)= ZERO
          BB(4,IX)= dNdx(INODE,2)
          BB(5,IX)= dNdx(INODE,3)
          BB(6,IX)= ZERO
          BB(1,IY)= ZERO
          BB(2,IY)= dNdx(INODE,2)
          BB(3,IY)= ZERO
          BB(4,IY)= dNdx(INODE,1)
          BB(5,IY)= ZERO
          BB(6,IY)= dNdx(INODE,3)
          BB(1,IZ)= ZERO
          BB(2,IZ)= ZERO
          BB(3,IZ)= dNdx(INODE,3)
          BB(4,IZ)= ZERO
          BB(5,IZ)= dNdx(INODE,1)
          BB(6,IZ)= dNdx(INODE,2)
         END DO
        ENDIF	   
!       
!     ==================================================================
!     Nodal displacements
!     ==================================================================
        DO J=1,MCRD
         ULOC(J)=ZERO
        END DO
        DO J=1,MCRD
         DO I=1,NDOFEL
          ULOC(J)=ULOC(J)+VNI(J,I)*U(I)
         END DO
        END DO  
!   
!     ==================================================================
!     Nodal phase-field
!     ==================================================================
        IF ((STEPITER.EQ.ZERO).AND.(REITER.EQ.ZERO)) THEN
         PHASE=USRVAR(JELEM,NSTVTT+1,INPT)	
        ELSE
         PHASE=USRVAR(JELEM,1,INPT)
        ENDIF
        IF (PHASE.GT.ONE) THEN
         PHASE=ONE
        ELSEIF (PHASE.LT.ZERO) THEN
         PHASE=ZERO
        ENDIF
!
!        SDV(1)=PHASE
        PHASE0=SDV(1) 
        IF (PHASE.LT.SDV(1)) THEN
         PHASE=SDV(1)
        ENDIF
        SDV(1)=PHASE
		
        DPHASE0=ZERO
        DPHASE0=PHASE-PHASE0
        IF (DPHASE0.LT.ZERO) THEN
         DPHASE0=ZERO
        ENDIF
		
!        write(7,*) 'PHASEE',PHASE	
!
!     ==================================================================
!     Calculating materials stiffness matrix
!     ==================================================================
        IF (MCRD.EQ.2) THEN
! --------- 2D Elements ------------ 
         DO I=1,3
          DO J=1,3
           CMAT(I,J)=ZERO
          END DO
         END DO
         CMAT(1,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
         CMAT(2,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
         CMAT(3,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
         CMAT(1,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
         CMAT(2,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
!		
        ELSE	   
! --------- 3D Elements ------------       
         DO I=1,6
          DO J=1,6
           CMAT(I,J)=ZERO
          END DO
         END DO
         CMAT(1,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
         CMAT(2,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
         CMAT(3,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(ONE-ENU)
         CMAT(1,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
         CMAT(2,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
         CMAT(1,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
         CMAT(3,1)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
         CMAT(2,3)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
         CMAT(3,2)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*ENU
         CMAT(4,4)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
         CMAT(5,5)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
         CMAT(6,6)=EMOD/((ONE+ENU)*(ONE-TWO*ENU))*(HALF-ENU)
        ENDIF
!		
!     ==================================================================
!     Calculating strain
!     ==================================================================
        IF (MCRD.EQ.2) THEN
! --------- 2D Elements ------------ 
         NP=3		
        ELSE	   
! --------- 3D Elements ------------       
         NP=6
        ENDIF
	   
        DO J=1,NP
         EPS(J)=ZERO
        END DO
        DO I=1,NP
         DO J=1,NDOFEL
          EPS(I)=EPS(I)+BB(I,J)*U(J)    
         END DO
        END DO
!
        DO K1=1,6
         EPSC(K1)=ZERO
        END DO		
        IF (MCRD.EQ.2) THEN
! --------- 2D Elements ------------ 
         EPSC(1)=EPS(1)
         EPSC(2)=EPS(2)
         EPSC(4)=EPS(3) 		 
        ELSE	   
! --------- 3D Elements ------------       
         EPSC=EPS
        ENDIF		
!		
!     ==================================================================
!     Calculating stresses
!     ==================================================================		
        DO K1=1,NP
         STRESS(K1)=ZERO
        END DO
        DO K1=1,NP
         DO K2=1,NP
          STRESS(K2)=STRESS(K2)+CMAT(K2,K1)*EPS(K1)
         END DO
        END DO
		
!      Calculating von Mises Stress		
        IF (MCRD.EQ.2) THEN
! --------- 2D Elements ------------ 
        VM_STRESS=(HALF*((STRESS(1)-STRESS(2))**TWO+STRESS(1)**TWO+ 
     1  STRESS(2)**TWO+SIX*(STRESS(3)**TWO)))**HALF	 
!		
        ELSE	   
! --------- 3D Elements ------------  
        VM_STRESS=(HALF*((STRESS(1)-STRESS(2))**TWO+
     1  (STRESS(2)-STRESS(3))**TWO+(STRESS(3)-STRESS(1))**TWO+ 
     2  SIX*(STRESS(4)**TWO+STRESS(5)**TWO+STRESS(6)**TWO)))**HALF	 
        ENDIF
		
        SDV(3)=VM_STRESS								
!
!     ==================================================================
!     Energy degradation function
!     ==================================================================
!
        OMEGA   =  (ONE-PHASE)**TWO+PARK        
        DOMEGA  = -TWO*(ONE-PHASE)
        DDOMEGA =  TWO
!
!     ==================================================================
!     Calculating elastic ENERGY
!     ==================================================================
!        ENG=ZERO
!        DO I=1,NP
!         ENG=ENG+STRESS(I)*EPS(I)*HALF
!        END DO

        ENG0=HALF*FT**TWO/EMOD

        CALL EIGOWN(EPSC,EIGV,ALPHAE,ALPHAI,VECTI)
!       
        ENG=(ELAMEL*(ALPHAE*(EIGV(1)+EIGV(2)+EIGV(3)))**TWO)/
     1       TWO+ELAMEG*((EIGV(1)*ALPHAI(1))**TWO+(EIGV(2)*
     2       ALPHAI(2))**TWO+(EIGV(3)*ALPHAI(3))**TWO)
	 
        ENGMAX=max(ENG,ENG0)	 

        ENERGY(2)=ENG	
        
!        write(7,*) 'ENG',ENG	 
!
!     ==================================================================
!     Adaptive time step control
!     ==================================================================
!
        DENG=ENGMAX-SDV(2)
!        write(7,*) 'DENG',DENG
!        write(7,*) 'ENGMAX',ENGMAX	
!        write(7,*) 'SDV(2)',SDV(2)			
!
        DENGMAXV=GCPAR		
        IF ((DPHASE0.GT.DPHMIN).AND.(REITER.LT.4)) THEN 
          IF ((DENG.GT.DENGMAXV)) THEN
            IF (DTIME.GT.DTMIN*(ONE+TOLER)) THEN
              PNEWDTIP(INPT)=DENGMAXV/DENG/TEN ! PNEWDTIP(INPT)=DENGMAXV/DENG
              DTIMENEXT=PNEWDTIP(INPT)*DTIME
              IF (DTIMENEXT.LT.(DTMIN*(ONE+TOLER))) THEN
                PNEWDTIP(INPT)=ONE 
              ENDIF
            ELSE
            PNEWDTIP(INPT)=ONE 
            ENDIF
          ENDIF
        ENDIF
!		
        SDV(2)=max(ENG,ENG0)
!
!     ==================================================================
!     Calculating element stiffness matrix
!     ==================================================================
!
        DO K=1,NDOFEL
         DO L=1,NDOFEL
          DO I=1,NP
           DO J=1,NP
            AMATRX(K,L)=AMATRX(K,L)+AINTW(INPT)*BB(I,K)*CMAT(I,J)*
     1       BB(J,L)*DTM*OMEGA*THCK
           END DO
          END DO
         END DO
        END DO
!       
!     ==================================================================
!     Internal forces (residual vector)
!     ==================================================================
        DO K1=1,NDOFEL
         DO K4=1,NP
           RHS(K1,1)=RHS(K1,1)-AINTW(INPT)*BB(K4,K1)*STRESS(K4)*DTM*
     1      OMEGA*THCK
         END DO
        END DO
!       
!     ==================================================================
!     Uploading solution dep. variables
!     ==================================================================
        DO I=1,NSTVTT
         SVARS(NSTVTT*(INPT-1)+I)=SDV(I)
         USRVAR(JELEM,I,INPT)=SVARS(NSTVTT*(INPT-1)+I)
        END DO
       END DO
!        
!     New time increment
       PNEWDTE=MINVAL(PNEWDTIP)
       IF (PNEWDTE.LT.(ONE+TOLER)) THEN
        PNEWDT=PNEWDTE
       ENDIF
!      
!     ==================================================================
!     ******************************************************************
!     Constructing elemet TYPE 5 and 6 and 7 and 8 (damage phase-field)
!     ******************************************************************
!     ==================================================================
      ELSEIF ((JTYPE.EQ.FIVE).OR.(JTYPE.EQ.SIX).OR.
     1       (JTYPE.EQ.SEVEN).OR.(JTYPE.EQ.EIGHT)) THEN
       REITER=USRVAR(JELEM-N_ELEM,NALL+3,1)
       STEPITER=USRVAR(JELEM-N_ELEM,NALL+4,1)
!     ==================================================================
!     Material parameters
!     ==================================================================
       EMOD =PROPS(1)
       ATPAR=PROPS(2)
       CLPAR=PROPS(3)
       GCPAR=PROPS(4)
       PARK =1.D-6 	
       THCK =ONE
!     ==================================================================
!     Crack geometry function 
!     ==================================================================	   
!        	  
!      The critical strength upon crack nucleation (not
!      necessarily the peak strength	  
       IF (NINT(ATPAR)==0) THEN ! AT2
         FT = ZERO
         C0 = HALF
       ELSEIF (NINT(ATPAR)==1) THEN ! AT1
         FT = SQRT(RP75*EMOD*GCPAR/CLPAR)
         C0 = TWO/THREE
       ELSE
         WRITE (*,*) '**ERROR: MODEL LAW NO. ', ATPAR, 
     1                 'DOES NOT EXIST!'
         CALL XIT
       ENDIF 
!
!     ==================================================================
!     Initial preparations
!     ==================================================================
       DO K1 = 1, NDOFEL                      
        DO KRHS = 1, NRHS
         RHS(K1,KRHS) = ZERO
        END DO
        DO K2 = 1, NDOFEL
         AMATRX(K2,K1) = ZERO
        END DO
       END DO
!     ==================================================================
!     Local coordinates and weights
!     ==================================================================
       IF (JTYPE.EQ.FIVE) THEN
        XII(1,1) = ONE/THREE
        XII(1,2) = ONE/THREE
        INNODE = 1.D0
        AINTW(1) = HALF
       ELSEIF (JTYPE.EQ.SIX) THEN
        XII(1,1) = -ONE/THREE**HALF
        XII(1,2) = -ONE/THREE**HALF
        XII(2,1) = ONE/THREE**HALF
        XII(2,2) = -ONE/THREE**HALF
        XII(3,1) = ONE/THREE**HALF
        XII(3,2) = ONE/THREE**HALF
        XII(4,1) = -ONE/THREE**HALF
        XII(4,2) = ONE/THREE**HALF
        INNODE = 4.D0
        DO I=1,INNODE
         AINTW(I) = ONE
        END DO	   
       ELSEIF (JTYPE.EQ.SEVEN) THEN
        XII(1,1) = ONE/FOUR
        XII(1,2) = ONE/FOUR
        XII(1,3) = ONE/FOUR
        INNODE=1.D0
        DO I=1,INNODE
         AINTW(I) = ONE/SIX
        END DO
       ELSEIF (JTYPE.EQ.EIGHT) THEN
        XII(1,1) = MONE/THREE**HALF
        XII(1,2) = MONE/THREE**HALF
        XII(1,3) = MONE/THREE**HALF
        XII(2,1) = ONE/THREE**HALF
        XII(2,2) = MONE/THREE**HALF
        XII(2,3) = MONE/THREE**HALF
        XII(3,1) = ONE/THREE**HALF
        XII(3,2) = ONE/THREE**HALF
        XII(3,3) = MONE/THREE**HALF
        XII(4,1) = MONE/THREE**HALF
        XII(4,2) = ONE/THREE**HALF
        XII(4,3) = MONE/THREE**HALF
        XII(5,1) = MONE/THREE**HALF
        XII(5,2) = MONE/THREE**HALF
        XII(5,3) = ONE/THREE**HALF
        XII(6,1) = ONE/THREE**HALF
        XII(6,2) = MONE/THREE**HALF
        XII(6,3) = ONE/THREE**HALF
        XII(7,1) = ONE/THREE**HALF
        XII(7,2) = ONE/THREE**HALF
        XII(7,3) = ONE/THREE**HALF
        XII(8,1) = MONE/THREE**HALF
        XII(8,2) = ONE/THREE**HALF
        XII(8,3) = ONE/THREE**HALF
        INNODE = 8.D0
        DO I=1,INNODE
         AINTW(I) = ONE
        END DO
       ENDIF
!
!     ==================================================================
!     Calculating properties at each integration point
!     ==================================================================
       DO INPT=1,INNODE
!     Initializing solution dependent variables (phase,history)
        DO I=1,NSTVTO
          SDV(I)=SVARS(NSTVTO*(INPT-1)+I)
        END DO
!
!     Local coordinates of the integration point 
        DO I=1,MCRD
         XI(I) = XII(INPT,I)
        END DO		
!     Shape functions and local derivatives
        IF (JTYPE.EQ.FIVE) THEN
         CALL SHAPEFUNTRI(AN,dNdxi,XI)
        ELSEIF (JTYPE.EQ.SIX) THEN
         CALL SHAPEFUNQUAD(AN,dNdxi,XI)
        ELSEIF (JTYPE.EQ.SEVEN) THEN
         CALL SHAPEFUNTET(AN,dNdxi,XI)
        ELSEIF (JTYPE.EQ.EIGHT) THEN
         CALL SHAPEFUNBRICK(AN,dNdxi,XI)
        ENDIF
		
!     Jacobian
        DO I = 1,MCRD
         DO J = 1,MCRD
          VJACOB(I,J) = ZERO
          DO K = 1,NNODE
           VJACOB(I,J) = VJACOB(I,J) + COORDS(I,K)*dNdxi(K,J)
          END DO
         END DO
        END DO
!        
        DTM = ZERO
        IF (MCRD.EQ.2) THEN
! --------- 2D Elements ------------        
         DTM = VJACOB(1,1)*VJACOB(2,2)-VJACOB(1,2)*VJACOB(2,1)
        ELSE
! --------- 3D Elements ------------        
         DTM = VJACOB(1,1)*VJACOB(2,2)*VJACOB(3,3)+VJACOB(1,2)*
     1   VJACOB(2,3)*VJACOB(3,1)+VJACOB(1,3)*VJACOB(2,1)*
     2   VJACOB(3,2)-VJACOB(3,1)*VJACOB(2,2)*VJACOB(1,3)-
     3   VJACOB(3,2)*VJACOB(2,3)*VJACOB(1,1)-VJACOB(3,3)*
     4   VJACOB(2,1)*VJACOB(1,2)        
        ENDIF	 	 
!     
        IF (DTM.LT.ZERO) THEN
         WRITE(7,*) 'Negative Jacobian',DTM
         CALL XIT	
        ENDIF
		
!     Inverse of Jacobian
        IF (MCRD.EQ.2) THEN
! --------- 2D Elements ------------        
         VJABOBINV(1,1)= VJACOB(2,2)/DTM
         VJABOBINV(1,2)=-VJACOB(1,2)/DTM
         VJABOBINV(2,1)=-VJACOB(2,1)/DTM
         VJABOBINV(2,2)= VJACOB(1,1)/DTM
        ELSE
! --------- 3D Elements ------------        
         VJABOBINV(1,1)=(VJACOB(2,2)*VJACOB(3,3)-VJACOB(2,3)*
     1    VJACOB(3,2))/DTM
         VJABOBINV(1,2)=-(VJACOB(1,2)*VJACOB(3,3)-VJACOB(3,2)*
     1    VJACOB(1,3))/DTM
         VJABOBINV(1,3)=(VJACOB(1,2)*VJACOB(2,3)-VJACOB(1,3)*
     1    VJACOB(2,2))/DTM
         VJABOBINV(2,1)=-(VJACOB(2,1)*VJACOB(3,3)-VJACOB(2,3)*
     1    VJACOB(3,1))/DTM
         VJABOBINV(2,2)=(VJACOB(1,1)*VJACOB(3,3)-VJACOB(1,3)*
     1    VJACOB(3,1))/DTM
         VJABOBINV(2,3)=-(VJACOB(1,1)*VJACOB(2,3)-VJACOB(1,3)*
     1    VJACOB(2,1))/DTM
         VJABOBINV(3,1)=(VJACOB(2,1)*VJACOB(3,2)-VJACOB(2,2)*
     1    VJACOB(3,1))/DTM
         VJABOBINV(3,2)=-(VJACOB(1,1)*VJACOB(3,2)-VJACOB(1,2)*
     1    VJACOB(3,1))/DTM
         VJABOBINV(3,3)=(VJACOB(1,1)*VJACOB(2,2)-VJACOB(1,2)*
     1    VJACOB(2,1))/DTM      
        ENDIF 
!        
!     Derivatives of shape functions respect to global ccordinates
        DO K = 1,NNODE
         DO I = 1,MCRD
          dNdx(K,I) = ZERO
          DO J = 1,MCRD
           dNdx(K,I) = dNdx(K,I) + dNdxi(K,J)*VJABOBINV(J,I)
          END DO
         END DO
        END DO
!
!     Calculating B matrix (B=LN)
       DO INODE=1,NNODE
        DO K1=1,MCRD
         BP(K1,INODE)=dNdx(INODE,K1)
        END DO
       END DO
!
!     ==================================================================
!     Nodal phase-field
!     ==================================================================
        PHASE=ZERO
        DPHASE=ZERO
        DO I=1,NDOFEL
         PHASE=PHASE+AN(I)*U(I)
        END DO
        DO I=1,NDOFEL
         DPHASE=DPHASE+AN(I)*DU(I,1)
        END DO
        SDV(1)=PHASE
!
!     Gradient
        DO I=1,MCRD
         DP(I)=ZERO
        END DO
        DO I=1,MCRD
         DO J=1,NNODE
          DP(I)=DP(I)+BP(I,J)*U(J)
         END DO
        END DO
!
!     ==================================================================
!     Energy degradation function
!     ==================================================================
!
        OMEGA   =  (ONE-PHASE)**TWO+PARK        
        DOMEGA  = -TWO*(ONE-PHASE)
        DDOMEGA =  TWO
!
!     ==================================================================
!     Calculating elastic energy history
!     ==================================================================
        IF ((STEPITER.EQ.ZERO).AND.(REITER.EQ.ZERO)) THEN
         ENGN=USRVAR(JELEM-N_ELEM,2,INPT)
        ELSE
         ENGN=USRVAR(JELEM-N_ELEM,NSTVTT+2,INPT)
        ENDIF
!       
        HISTN=USRVAR(JELEM-N_ELEM,NSTVTT+2,INPT)
        IF (ENGN.GT.HISTN) THEN
         HIST=ENGN
        ELSE
         HIST=HISTN
        ENDIF	
        SDV(2)=HIST	
!        write(7,*) 'PHASEH',PHASE	
!        write(7,*) 'HIST',HIST	
!        
!     ==================================================================
!     Calculating fracture energy for history output
!     ==================================================================
!
        ALPHA  =ATPAR*PHASE+(ONE-ATPAR)*PHASE**TWO  
        DALPHA =ATPAR+TWO*(ONE-ATPAR)*PHASE
        DDALPHA=TWO*(ONE-ATPAR)
!
        ENGD=ZERO
        ENGD=ALPHA/(FOUR*C0*CLPAR)*DTM*GCPAR
!
        DO J=1,MCRD
         ENGD=ENGD+DP(J)*DP(J)*CLPAR*DTM/(FOUR*C0)*GCPAR
        END DO
        ENGDG=ENGDG+ENGD
!
!     ==================================================================
!     Calculating element stiffness matrix (K_dd)
!     ==================================================================

        PHASE_SOURCE=DOMEGA*HIST+GCPAR/(FOUR*C0*CLPAR)*DALPHA
        DPHASE_SOURCE=DDOMEGA*HIST+GCPAR/(FOUR*C0*CLPAR)*DDALPHA

        DO I=1,NNODE
         DO K=1,NNODE
          DO J=1,MCRD
           AMATRX(I,K)=AMATRX(I,K)+BP(J,I)*BP(J,K)*DTM*
     1      GCPAR*CLPAR/(TWO*C0)*AINTW(INPT)*THCK
          END DO
          AMATRX(I,K)=AMATRX(I,K)+AN(I)*AN(K)*DTM*
     1     AINTW(INPT)*(DPHASE_SOURCE)*THCK
         END DO
        END DO
!        
!     ==================================================================
!     Internal forces (residual vector) (r_d)
!     ==================================================================
        DO I=1,NDOFEL
         DO J=1,MCRD
           RHS(I,1)=RHS(I,1)-BP(J,I)*DP(J)*GCPAR*CLPAR/(TWO*C0)*
     1      AINTW(INPT)*DTM*THCK
         END DO
         RHS(I,1)=RHS(I,1)-AN(I)*AINTW(INPT)*DTM*
     1    (PHASE_SOURCE)*THCK 
        END DO
!
!     ==================================================================
!     Uploading solution dep. variables
!     ==================================================================
        DO I=1,NSTVTO
         SVARS(NSTVTO*(INPT-1)+I)=SDV(I)
         USRVAR(JELEM-N_ELEM,I+NSTVTT,INPT)=SVARS(NSTVTO*(INPT-1)+I)
        END DO
       END DO
!       
      ENDIF
        ENERGY(7)=ENGDG
      RETURN
      END	  
!
! ==============================================================
! Eigenstrains from Voigt notation
! ==============================================================
!
      SUBROUTINE EIGOWN(EPS,EIGV,ALPHAE,ALPHAI,VECTI)
      INCLUDE 'ABA_PARAM.INC'
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,TWO=2.D0,
     1  TS=27.D0,THREE=3.D0,HALF=0.5D0,TOLER=1.0D-12,FOUR=4.D0,
     2  CNTN=100,TOLERE=1.0D-12,PI=3.1415926535897932384626433832D0)
      INTEGER I, J, K
      REAL*8 EPS(6), EIGV(3), ALPHAI(3), VECTI(3)
      REAL*8 PC, QC, ALPHAE, DISC, PI, CNT
!
!   Scaling the strain vector
       VMAXE=MAXVAL(ABS(EPS))
       IF (VMAXE.GT.TOLERE) THEN
        DO K1=1,6
         EPS(K1)=EPS(K1)/VMAXE
        END DO
       ENDIF
!    
!   Calculating eigenvalues
       VECTI(1)=EPS(1)+EPS(2)+EPS(3)
       VECTI(2)=EPS(1)*EPS(2)+EPS(2)*EPS(3)+EPS(1)*EPS(3)-
     1 EPS(4)**TWO/FOUR-EPS(5)**TWO/FOUR-EPS(6)**TWO/FOUR
       VECTI(3)=EPS(1)*EPS(2)*EPS(3)+EPS(4)*EPS(5)*EPS(6)/
     1 FOUR-EPS(1)*EPS(6)**TWO/FOUR-EPS(2)*EPS(5)**TWO/FOUR-
     2 EPS(3)*EPS(4)**TWO/FOUR
!
!   Depressed coefficients    
       PC=VECTI(2)-VECTI(1)**TWO/THREE
       QC=VECTI(1)*VECTI(2)/THREE-TWO*VECTI(1)**THREE/TS-VECTI(3)
       DISC=MONE*(FOUR*PC**THREE+TS*QC**TWO)
!
       DO I=1,3
        EIGV(I)=ZERO
       END DO
       CNT=ZERO
       IF (ABS(DISC).LT.TOLER) THEN
        IF ((ABS(QC).LT.TOLER).AND.(ABS(PC).LT.TOLER)) THEN
         EIGV(1)=VECTI(1)/THREE
         EIGV(2)=VECTI(1)/THREE
         EIGV(3)=VECTI(1)/THREE
        ELSE
         EIGV(1)=-THREE*QC/TWO/PC+VECTI(1)/THREE
         EIGV(2)=-THREE*QC/TWO/PC+VECTI(1)/THREE
         EIGV(3)=THREE*QC/PC+VECTI(1)/THREE
         IF (EIGV(1).GT.EIGV(3)) THEN
          EONE=EIGV(1)
          EIGV(1)=EIGV(3)
          EIGV(3)=EONE
         ENDIF
        ENDIF
       ELSE
        DO I=1,3
         EIGV(I)=VECTI(1)/THREE+TWO*(MONE*PC/THREE)**HALF*
     1   COS(ONE/THREE*ACOS(MONE*QC/TWO*(TS/(MONE*PC**THREE))**
     2   HALF)+TWO*I*PI/THREE)
        END DO
       ENDIF
!       
       ALPHAE=ZERO
       IF ((EIGV(1)+EIGV(2)+EIGV(3)).GT.TOLER) THEN
        ALPHAE=ONE
       ENDIF
       DO K1=1,3
        ALPHAI(K1)=ZERO
        IF (EIGV(K1).GT.TOLER) THEN
         ALPHAI(K1)=ONE
        ENDIF
       END DO
!
!    Rescaling eigenvalues       
       IF (VMAXE.GT.TOLERE) THEN
        DO K1=1,6
         EPS(K1)=EPS(K1)*VMAXE
        END DO
        DO K1=1,3
         EIGV(K1)=EIGV(K1)*VMAXE
        END DO
         VECTI(1)=EIGV(1)+EIGV(2)+EIGV(3)
         VECTI(2)=EIGV(1)*EIGV(2)+EIGV(1)*EIGV(3)+EIGV(3)*EIGV(2)
         VECTI(3)=EIGV(1)*EIGV(2)*EIGV(3)
       ENDIF
!   
       RETURN
       END
!
! ======================================================================
! 3-node linear elements
! 4-node bilinear elements
! ======================================================================
!	 	          
! ----------------------------------------------------------------------      
! Shape functions for 3-node linear elements
! ----------------------------------------------------------------------      
      SUBROUTINE SHAPEFUNTRI(AN,dNdxi,XI)
      INCLUDE 'ABA_PARAM.INC'
      Real*8 AN(4),dNdxi(3,2)
      Real*8 XI(2)
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,FOUR=4.D0)

!     Values of shape functions as a function of local coord.
      AN(1) = XI(1)
      AN(2) = XI(2)
      AN(3) = ONE-XI(1)-XI(2)
!
!     Derivatives of shape functions respect to local ccordinates
      DO I=1,4
        DO J=1,2
            dNdxi(I,J) =  ZERO
        END DO
      END DO
      dNdxi(1,1) =  ONE
      dNdxi(1,2) =  ZERO
      dNdxi(2,1) =  ZERO
      dNdxi(2,2) =  ONE
      dNdxi(3,1) =  MONE
      dNdxi(3,2) =  MONE
      RETURN
      END
!
! -----------------------------------------------------------------------      
! Shape functions for 4-node bilinear elements
! -----------------------------------------------------------------------      
      SUBROUTINE SHAPEFUNQUAD(AN,dNdxi,xi)
      INCLUDE 'ABA_PARAM.INC'
      Real*8 AN(4),dNdxi(4,2)
      Real*8 XI(2)
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,FOUR=4.D0)
!
!     Values of shape functions as a function of local coord.
      AN(1) = ONE/FOUR*(ONE-XI(1))*(ONE-XI(2))
      AN(2) = ONE/FOUR*(ONE+XI(1))*(ONE-XI(2))
      AN(3) = ONE/FOUR*(ONE+XI(1))*(ONE+XI(2))
      AN(4) = ONE/FOUR*(ONE-XI(1))*(ONE+XI(2))
!
!     Derivatives of shape functions respect to local coordinates
      DO I=1,4
        DO J=1,2
            dNdxi(I,J) =  ZERO
        END DO
      END DO
      dNdxi(1,1) =  MONE/FOUR*(ONE-XI(2))
      dNdxi(1,2) =  MONE/FOUR*(ONE-XI(1))
      dNdxi(2,1) =  ONE/FOUR*(ONE-XI(2))
      dNdxi(2,2) =  MONE/FOUR*(ONE+XI(1))
      dNdxi(3,1) =  ONE/FOUR*(ONE+XI(2))
      dNdxi(3,2) =  ONE/FOUR*(ONE+XI(1))
      dNdxi(4,1) =  MONE/FOUR*(ONE+XI(2))
      dNdxi(4,2) =  ONE/FOUR*(ONE-XI(1))
      RETURN
      END
!
! ======================================================================
! 4-node linear tetrahedron elements
! 8-node linear brick elements
! ======================================================================
! 
! ----------------------------------------------------------------------      
! Shape functions for 4-node linear tetrahedron elements
! ----------------------------------------------------------------------  
      SUBROUTINE SHAPEFUNTET(AN,dNdxi,XI)
      INCLUDE 'ABA_PARAM.INC'
      Real*8 AN(4),dNdxi(4,3)
      Real*8 XI(3)
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,FOUR=4.D0)

!     Values of shape functions as a function of local coord.
      AN(1) = ONE-XI(1)-XI(2)-XI(3)
      AN(2) = XI(1)
      AN(3) = XI(2)
      AN(4) = XI(3)
!
!     Derivatives of shape functions respect to local coordinates
      DO I=1,4
        DO J=1,3
            dNdxi(I,J) =  ZERO
        END DO
      END DO
      dNdxi(1,1) =  MONE
      dNdxi(1,2) =  MONE
      dNdxi(1,3) =  MONE
!
      dNdxi(2,1) =  ONE
      dNdxi(2,2) =  ZERO
      dNdxi(2,3) =  ZERO
!
      dNdxi(3,1) =  ZERO
      dNdxi(3,2) =  ONE
      dNdxi(3,3) =  ZERO
!
      dNdxi(4,1) =  ZERO
      dNdxi(4,2) =  ZERO
      dNdxi(4,3) =  ONE
!
      RETURN
      END
!
! --------------------------------------------------------------------      
! Shape functions for 8-node linear brick elements
! --------------------------------------------------------------------     
      SUBROUTINE SHAPEFUNBRICK(AN,dNdxi,XI)
      INCLUDE 'ABA_PARAM.INC'
      REAL*8 AN(8),dNdxi(8,3)
      REAL*8 XI(3)
      PARAMETER(ZERO=0.D0,ONE=1.D0,MONE=-1.D0,FOUR=4.D0,EIGHT=8.D0)

!     Values of shape functions as a function of local coord.
      AN(1) = ONE/EIGHT*(ONE-XI(1))*(ONE-XI(2))*(ONE-XI(3))
      AN(2) = ONE/EIGHT*(ONE+XI(1))*(ONE-XI(2))*(ONE-XI(3))
      AN(3) = ONE/EIGHT*(ONE+XI(1))*(ONE+XI(2))*(ONE-XI(3))
      AN(4) = ONE/EIGHT*(ONE-XI(1))*(ONE+XI(2))*(ONE-XI(3))
      AN(5) = ONE/EIGHT*(ONE-XI(1))*(ONE-XI(2))*(ONE+XI(3))
      AN(6) = ONE/EIGHT*(ONE+XI(1))*(ONE-XI(2))*(ONE+XI(3))
      AN(7) = ONE/EIGHT*(ONE+XI(1))*(ONE+XI(2))*(ONE+XI(3))
      AN(8) = ONE/EIGHT*(ONE-XI(1))*(ONE+XI(2))*(ONE+XI(3))
      
!     Derivatives of shape functions respect to local coordinates
      DO I=1,8
        DO J=1,3
            dNdxi(I,J) =  ZERO
        END DO
      END DO
      dNdxi(1,1) =  MONE/EIGHT*(ONE-XI(2))*(ONE-XI(3))
      dNdxi(1,2) =  MONE/EIGHT*(ONE-XI(1))*(ONE-XI(3))
      dNdxi(1,3) =  MONE/EIGHT*(ONE-XI(1))*(ONE-XI(2))
      dNdxi(2,1) =  ONE/EIGHT*(ONE-XI(2))*(ONE-XI(3))
      dNdxi(2,2) =  MONE/EIGHT*(ONE+XI(1))*(ONE-XI(3))
      dNdxi(2,3) =  MONE/EIGHT*(ONE+XI(1))*(ONE-XI(2))
      dNdxi(3,1) =  ONE/EIGHT*(ONE+XI(2))*(ONE-XI(3))
      dNdxi(3,2) =  ONE/EIGHT*(ONE+XI(1))*(ONE-XI(3))
      dNdxi(3,3) =  MONE/EIGHT*(ONE+XI(1))*(ONE+XI(2))
      dNdxi(4,1) =  MONE/EIGHT*(ONE+XI(2))*(ONE-XI(3))
      dNdxi(4,2) =  ONE/EIGHT*(ONE-XI(1))*(ONE-XI(3))
      dNdxi(4,3) =  MONE/EIGHT*(ONE-XI(1))*(ONE+XI(2))
      dNdxi(5,1) =  MONE/EIGHT*(ONE-XI(2))*(ONE+XI(3))
      dNdxi(5,2) =  MONE/EIGHT*(ONE-XI(1))*(ONE+XI(3))
      dNdxi(5,3) =  ONE/EIGHT*(ONE-XI(1))*(ONE-XI(2))
      dNdxi(6,1) =  ONE/EIGHT*(ONE-XI(2))*(ONE+XI(3))
      dNdxi(6,2) =  MONE/EIGHT*(ONE+XI(1))*(ONE+XI(3))
      dNdxi(6,3) =  ONE/EIGHT*(ONE+XI(1))*(ONE-XI(2))
      dNdxi(7,1) =  ONE/EIGHT*(ONE+XI(2))*(ONE+XI(3))
      dNdxi(7,2) =  ONE/EIGHT*(ONE+XI(1))*(ONE+XI(3))
      dNdxi(7,3) =  ONE/EIGHT*(ONE+XI(1))*(ONE+XI(2))
      dNdxi(8,1) =  MONE/EIGHT*(ONE+XI(2))*(ONE+XI(3))
      dNdxi(8,2) =  ONE/EIGHT*(ONE-XI(1))*(ONE+XI(3))
      dNdxi(8,3) =  ONE/EIGHT*(ONE-XI(1))*(ONE+XI(2))
!      
      RETURN
      END
!
!      
! Subroutine UMAT  : 
! Dummy material
!
! ==============================================================
! !!! NOTE: N_ELEM has to be changed according to the UEL !!!!!
! ==============================================================
!
       SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,STRAN,DSTRAN,
     2 TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,MATERL,NDI,NSHR,NTENS,
     3 NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,CELENT,
     4 DFGRD0,DFGRD1,NOEL,NPT,KSLAY,KSPT,KSTEP,KINC)
!
      INCLUDE 'ABA_PARAM.INC'
!
       CHARACTER*80 CMNAME
       DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),
     2 DDSDDT(NTENS),DRPLDE(NTENS),
     3 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     4 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
! 
       PARAMETER(zero=0.D0,ONE=1.D0,TWO=2.D0,THREE=3.D0,SIX=6.D0, 
     1 HALF=0.5D0,N_ELEM=82406,NSTV=9) 
!	 
       DATA NEWTON,TOLER/40,1.D-6/ 
!       
       COMMON/KUSER/USRVAR(N_ELEM,NSTV,8)
! 
! ----------------------------------------------------------- 
!
!   Stiffness tensor
!
       DDSDDE=0.D0
!
       NELEMAN=NOEL-TWO*N_ELEM
       IF (NPT.EQ.3) THEN
        NPT=4
       ELSEIF (NPT.EQ.4) THEN
        NPT=3
       ENDIF
       IF (NPT.EQ.7) THEN
        NPT=8
       ELSEIF (NPT.EQ.8) THEN
        NPT=7
       ENDIF
       
       DO I=1,NSTATV
        STATEV(I)=USRVAR(NELEMAN,I,NPT)
       END DO
       
       RETURN
       END      
!
! ------------------------------------------------------------

      