!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Adrian Soto
!  24-03-2016
!  Stony Brook University
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This file is part of the code QEdark v1.1.1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
SUBROUTINE qedark_f2_3d( restartmode, &
     nksf, numval, numcond, &
     vearth_SI, vesc_SI, v0_SI, deltav_SI, &
     Er_bin_type, num_er_bins, ermax_NU, er_binsize, &
     numqbins, dq, &
     do_scissor_correction, scissorgap)
  !
  ! Main driver to evaluate and print the spherically form factor
  ! squared as a function of |q| and E.
  ! 
  ! 
  !
  !
  USE constants, ONLY: rytoev
  USE wavefunctions_module,           ONLY: evc     ! For collinear, evc(npwx, nbnd) [look at allocate_wfc.f90]
  USE kinds,                          ONLY: DP
  USE wvfct,                          ONLY: igk, nbnd, npwx, et,g2kin ,ecutwfc !, btype
  USE klist,                          ONLY: nks, ngk, wk, xk, nelec
  USE lsda_mod,                       ONLY: nspin
  USE io_files,                       ONLY: nwordwfc, iunwfc, iunigk
  USE buffers,                        ONLY: get_buffer
  USE gvect,                          ONLY: g
  USE cell_base,                      ONLY: bg, tpiba, tpiba2, omega
  USE noncollin_module,               ONLY: noncolin

  use omp_lib

  
  IMPLICIT NONE

  INTEGER, EXTERNAL :: find_bin
  REAL(DP), EXTERNAL :: vmin, bzvol, eucl_norm_fast
  

  ! Variables specific to int_formfact(..)
  REAL(DP), PARAMETER :: Ry2eV = 13.60569253_DP                   ! ==2/alpha 
  REAL(DP), PARAMETER :: twobyalpha=274.07199814110116            ! ==2/alpha, which is the conversion factor for velocities from N.U. to Rydberg A.U. 
  REAL(DP), PARAMETER :: speedoflight=3.0E08_DP                   ! Speed of light in SI units 
  REAL(DP), PARAMETER :: RAU2kmps=(speedoflight/1000.0_DP)/twobyalpha 
  REAL(DP), PARAMETER :: pi = 3.141592653589793238462643383_DP
  
  LOGICAL :: restartmode
  LOGICAL :: f1exists, f2exists, f3exists
  INTEGER :: ik1init, ik2init
  
  LOGICAL :: do_scissor_correction != .false.                      ! 
  REAL(DP) :: scissorgap                                          ! Band gap (in eV) for scissor operator corrected band energies 


  REAL(DP) :: vearth_SI                                           ! Earth's average velocity around the galactic center in SI units 
  REAL(DP) :: vearth_RAU !=(twobyalpha/speedoflight)*vearth_SI      ! Earth's average velocity in RAU
  REAL(DP) :: vearth_kmps !=vearth_SI/1000.0                        ! Earth's average velocity in km/s

  REAL(DP) :: vesc_SI                                             ! Earth's escape velocity in SI units 
  REAL(DP) :: vesc_kmps != vesc_SI/1000.0                            ! Earth's escape velocity in km/s
  REAL(DP) :: vesc_RAU !=(twobyalpha/speedoflight)*vesc_SI          ! Earth's escape velocity in RAU
  INTEGER :: iparastart
  REAL(DP) :: v0_SI                                                 ! DM typical velocity in SI
  REAL(DP) :: v0_kmps !=v0_SI/1000.0                                ! DM typical velocity in km/s
  REAL(DP) :: v0_RAU !=(twobyalpha/speedoflight)*v0_SI              ! DM typical velocity in RAU
  REAL(DP) :: RAU2NU


  INTEGER :: idmff, imonth
  INTEGER :: nmonths = 3  
  REAL(DP) :: deltav_SI(3)
  REAL(DP) :: deltav_kmps(3) !=v0_SI/1000.0                                
  REAL(DP) :: deltav_RAU(3) !=(twobyalpha/speedoflight)*v0_SI              
  
  
  REAL(DP) :: me_NU = 0.510998910E6_DP                             ! Electron mass in NU 
  

  INTEGER, PARAMETER :: max_num_mx=99
  INTEGER :: imx, num_mx
  REAL(DP) :: mx_NU(max_num_mx) != 100.0E6_DP                     ! Dark matter particle mass in NU 
  REAL(DP) :: mx_RAU(max_num_mx) != 0.5_DP*(mx_NU/me_NU)          ! Dark matter particle mass in RAU (me=0.5)


  REAL(DP) :: vmin_aux(max_num_mx)                                ! auxiliary vmin for looping

  REAL(DP) :: bb(6)                                               ! Dot products of reciprocal lattice basis vectors
  REAL(DP) :: bzv                                                 ! 1BZ volume

  REAL(DP) :: fcross
  REAL(DP) :: q(3)                                                ! Coordinates of q vector
  REAL(DP) :: qnorm                                               ! Norm of q vector
  REAL(DP) :: wq                                                  ! weight of q-point for q-space integration
  REAL(DP) :: deltaEmax  
  
  INTEGER :: er_bin_type                                          ! Type of bins for integrated form factors
  INTEGER :: num_er_bins                                          ! Number of bins for integrated form factors
  REAL(DP) :: er_binsize                                          ! Bin size (Ry) for linear bins with user-provided size 

  INTEGER :: iE ! for Energy bin  
  REAL(DP) :: deltaE
  REAL(DP) :: ermax_NU !=0.0_DP                                   ! Max recoil energy cutoff in eV
  REAL(DP) :: ermax_RAU ! = ercut_NU / Ry2eV                      ! Max recoil energy cutoff in Ry
  REAL(DP), ALLOCATABLE :: binedgesE(:)
  
  CHARACTER(2) :: numbins                                         ! String containing number of bins. It can't be larger than 99
  CHARACTER(26) :: FMT                                            ! Format specifier
  
  
  REAL(DP) :: aux(max_num_mx)
  REAL(DP) :: aux_dmff(max_num_mx)
  REAL(DP) :: kgvec(3)  
  COMPLEX(DP) :: f1(4)                                               ! Current |f|  used for iteration, 4 possible spin pairs
  COMPLEX(DP) :: f2                                                  ! Current |f|^2  used for iteration
  !REAL(DP), ALLOCATABLE :: ctot(:,:), cbinned(:,:,:)              ! Integrated form factor C^{i-->i'}
  REAL(DP) :: nc                                               ! Norm check, nc(1)=sum_ig1 conjg(evc)*evc and nc(2)=sum_ig1 conjg(evcouter)*evcouter
  COMPLEX(DP) :: fvec(4,3)                                        ! Matrix with 4 values for each component of the f-vector, analogous to f(4)
  COMPLEX(DP) :: fvecs(3)                                         ! f-vector summed over spins
  COMPLEX(DP) :: f2vec(3)                                             
  COMPLEX(DP) :: f1temp(4)                                          

  REAL(DP), ALLOCATABLE :: ctot(:,:,:,:)                          ! Integrated form factor. No band indices
  REAL(DP), ALLOCATABLE :: ctotA1(:,:,:,:)                          ! Integrated form factor. No band indices
  REAL(DP), ALLOCATABLE :: ctotA2(:,:,:,:)                          ! Integrated form factor. No band indices
  REAL(DP), ALLOCATABLE :: ctotA3(:,:,:,:)                          ! Integrated form factor. No band indices
  REAL(DP), ALLOCATABLE :: ctotA1c(:,:,:,:)                          ! Integrated form factor. No band indices
  REAL(DP), ALLOCATABLE :: ctotA2c(:,:,:,:)                          ! Integrated form factor. No band indices
  REAL(DP), ALLOCATABLE :: ctotA3c(:,:,:,:)                          ! Integrated form factor. No band indices
  REAL(DP), ALLOCATABLE :: ctotf1(:,:,:,:)                          ! Integrated form factor. No band indices
  REAL(DP), ALLOCATABLE :: ctotf2(:,:,:,:)                          ! Integrated form factor. No band indices
  REAL(DP), ALLOCATABLE :: ctotf3(:,:,:,:)                          ! Integrated form factor. No band indices


  CHARACTER(20) :: outfile = "C.dat"                               ! Output file name 
  INTEGER :: if1, if2     
  INTEGER :: ik1, ik2                                             ! Indices for k-point loops 
  INTEGER :: ig1, ig2                                             ! Indices for G-vector loops
  INTEGER :: iband1, iband2                                       ! Indices for band loops 
  INTEGER :: ipara, ik3
  COMPLEX(DP), ALLOCATABLE :: evcouter(:,:)                       ! Auxiliar wavefunction array for looping. Will be used in the OUTERmost loop

  INTEGER, ALLOCATABLE :: alligk(:,:)                             ! Here we load from file the igk(:) for all k-points, i.e. alligk(ig, ik) = igk(ig)
  REAL(DP) :: dk(3,nks,nks)                                       ! Table storing all k1-k2 values  TODO: this table is antisymmetric and can be reduced
  INTEGER :: gsi(npwx, npwx)                                      ! Sum index table for band2band

  INTEGER :: nksf                                                 ! Number of k-points for formfactor calculation
  
  INTEGER :: numval, numcond                                      ! Number of occupied and unoccupied Kohn-Sham orbitals. TO BE SET BY USER!
  INTEGER :: numvaltot, numcondtot                                ! Total # bands in DFT run. numvaltot==nelec/2 and numcondtot=nbnd-numvaltot 
  INTEGER :: ivalbottom, ivaltop                                  ! Minimum and maximum values for valence band index
  INTEGER :: icondbottom, icondtop                                ! Minimum and maximum values for conduction index
  
  INTEGER :: ierr                                                 ! Error index

  REAL(DP) :: tol = 1.0E-6                                        ! Tolerance for G-vector distances


  INTEGER :: numqbins
  INTEGER :: iqx, iqy, iqz                                        ! for q bin
  REAL(DP) :: dq                                                  ! size of |q| bin
  REAL(DP), ALLOCATABLE :: binedgesq(:)  
  
  !DEBUGGING VARIABLES
  LOGICAL :: runff                                                ! For debugging purposes: set to false to skip the form factor evaluation                       

  runff = .true. ! Set to false to not skip ff calculation --for debugging purposes


  CALL start_clock( ' qedark_f2_3d ')

  print *, "           -------             "
  print *, " calculation_mode == f2_3d"
  print *, " Calculating formfactor squared binned in E and q"
  print *, "           -------             "



!  IF (nspin .ne. 1) THEN
!     CALL errore ('qedark_f2_3d', 'Form factor calculation works only for spin-unpolarized systems!', 1)
!  ENDIF


  IF ( nksf > nks .or. nksf < 0 ) &
       CALL errore( 'qedark_f2_3d ',' nksf has a non-allowed value. Check input. ', ABS(ierr) )


  ! Band indices
  IF ( noncolin .eqv. .false.) THEN
     numvaltot = nelec/2
  ELSE
     numvaltot = nelec
  ENDIF
  numcondtot= nbnd-numvaltot
  ivalbottom = numvaltot-numval+1
  ivaltop = numvaltot
  icondbottom = numvaltot+1
  icondtop = numvaltot + numcond



  IF( numval>numvaltot .or. numcond>numcondtot .or. numval<1 .or. numcond<1 ) &
       CALL errore( 'qedark_f2_3d ',' Check numval and numcond values in input ', ABS(ierr) )

  IF (numval /= numvaltot) THEN
     PRINT *, " "
     PRINT *, "  WARNING: this calculation will NOT iterate over all valence bands"
     PRINT *, " ivalbottom, ivaltop = ", ivalbottom, ivaltop
     PRINT *, " "
  ENDIF


  IF (numcond /= numcondtot) THEN
     PRINT *, " "
     PRINT *, "  WARNING: this calculation will NOT iterate over all conduction bands"
     PRINT *, " icondbottom, icondtop = ", icondbottom, icondtop
     PRINT *, " "
  ENDIF

    

  !CALL volume (tpiba, bg(:,1), bg(:,2), bg(:,3), bzv) ! could also do 
  bzv = bzvol(bg) 


  print *, "  npwx      ==", npwx
  print *, "  nbnd      ==", nbnd
  print *, "  numval    ==", numval
  print *, "  numvaltot ==", numvaltot
  print *, "  numcond   ==", numcond
  print *, "  numcondtot==", numcondtot
  PRINT *, "  cell vol  == ", omega, "bohr^(3)"
  PRINT *, "  1BZ vol   == ", bzv, "bohr^(-3)"


  !CALL qspace(bzv, .false.)



  IF (nspin .ne. 1) THEN
     CALL errore ('qedark_f2_3d', 'Form factor calculation works only for spin-unpolarized systems!', 1)
  ENDIF

  IF (noncolin .eqv. .false.) THEN
     ALLOCATE ( evcouter(npwx, nbnd) , STAT=ierr )
     IF( ierr /= 0 ) &
          CALL errore( 'qedark_f2_3d',' error allocating evcouter ', ABS(ierr) )

  ELSE
     ALLOCATE ( evcouter(2*npwx, nbnd) , STAT=ierr )
     IF( ierr /= 0 ) &
          CALL errore( 'qedark_f2_3d',' error allocating evcouter ', ABS(ierr) )
  ENDIF





  ALLOCATE ( alligk(npwx, nks) , STAT=ierr )
  IF( ierr /= 0 ) &
       CALL errore( 'qedark_f2_3d',' error allocating alligk ', ABS(ierr) )
  alligk(:,:)=0.0_DP ! when alligk=0, ig doesn't correspond to a G-vector in the set and should be disregarded 

  



  WRITE(*,*), " "


  IF (num_er_bins > 999) CALL errore( 'qedark_f2_3d','Number of energy recoil bins cannot exceed 999 ', ABS(ierr) )
  IF (numqbins > 999) CALL errore( 'qedark_f2_3d','Number of momentum transfer bins cannot exceed 999 ', ABS(ierr) )
  

  vearth_RAU = (twobyalpha/speedoflight)*vearth_SI
  vearth_kmps = vearth_SI/1000.0
  vesc_kmps = vesc_SI/1000.0
  vesc_RAU = (twobyalpha/speedoflight)*vesc_SI
  v0_kmps = v0_SI/1000.00                      
  v0_RAU = (twobyalpha/speedoflight)*v0_SI    
  deltav_kmps(:) = deltav_SI(:)/1000.0
  deltav_RAU(:) = (twobyalpha/speedoflight)*deltav_SI(:)
  mx_RAU(:) = 0.5_DP*(mx_NU(:)/me_NU)
  ermax_RAU = ermax_NU / Ry2eV 
  !RAU2NU = dsqrt(2*Ry2eV*me_NU)
  RAU2NU = me_NU/137.0
  !Calculating dq from E_cut and the number of q-bins according to eq. 4.3 of essig et al. The factor of 1.1 is added to overshoot E_cut since eq 4.3 is just an approximation.
  !dq in dm.in does not do anything and gets redefined here. The factor of 2 in front of the squareroot is because q has to cover both positive and negative values.
  dq=1.1*2*dsqrt(2*ecutwfc*Ry2eV*me_NU)/numqbins
!!!!!!!!!  CALL print_DM_data(vesc_kmps, vearth_kmps, v0_kmps, mx_NU, ermax_NU)


  print *, " "
  IF (do_scissor_correction) THEN
     CALL scissor(numvaltot, numcondtot, scissorgap, et)
  ENDIF


  ! Stuff for binning integrated form factors. 
  ! Needs to go after ermax_RAU is calculated
  IF (num_er_bins > 0) THEN
     
     ALLOCATE (binedgese(num_er_bins+1), STAT=ierr )
     IF( ierr /= 0 ) &
          CALL errore( 'qedark_f2_3d',' error allocating binedgesE ', ABS(ierr) )
     
     CALL create_bins(er_bin_type, 0.0_DP, ermax_RAU, &
          num_er_bins, er_binsize, binedgese)
     WRITE (*,*), " "
     WRITE (*,*) "Creating E bins for formfactor sum ..."
     WRITE (*,*), "bintype=", er_bin_type
     WRITE (*,*), "Energy bin size=", Ry2eV*er_binsize, "eV"
  ENDIF



  ! Stuff for binning integrated form factors. 
  ! Needs to go after ermax_RAU is calculated
  IF (num_er_bins > 0) THEN
     ALLOCATE (binedgesq(num_er_bins+1), STAT=ierr )
     IF( ierr /= 0 ) &
          CALL errore( 'qedark_f2_3d ',' error allocating binedgesq ', ABS(ierr) )

          
     CALL create_bins(er_bin_type, -dq*numqbins/2.0, dq*numqbins/2.0, &
          numqbins, dq, binedgesq)
     WRITE (*,*), " "
     WRITE (*,*) "Creating q bins for formfactor sum ..."
     WRITE (*,*), "bintype=", er_bin_type
     WRITE (*,*), "Momentum transfer bin size (tpiba)=", dq
  ENDIF



  ALLOCATE( ctot(numqbins+1, numqbins+1, numqbins+1, num_Er_bins+1) , STAT=ierr )
  IF( ierr /= 0 ) &
       CALL errore( 'qedark_f2_3d ',' cannot allocate ctot ', ABS(ierr) )
  ctot(:,:,:,:) = 0.0_DP
  
  ALLOCATE( ctotA1(numqbins+1, numqbins+1, numqbins+1, num_Er_bins+1) , STAT=ierr )
  IF( ierr /= 0 ) &
       CALL errore( 'qedark_f2_3d ',' cannot allocate ctot ', ABS(ierr) )
  ctotA1(:,:,:,:) = 0.0_DP

  ALLOCATE( ctotA2(numqbins+1, numqbins+1, numqbins+1, num_Er_bins+1) , STAT=ierr )
  IF( ierr /= 0 ) &
       CALL errore( 'qedark_f2_3d ',' cannot allocate ctot ', ABS(ierr) )
  ctotA2(:,:,:,:) = 0.0_DP

  ALLOCATE( ctotA3(numqbins+1, numqbins+1, numqbins+1, num_Er_bins+1) , STAT=ierr )
  IF( ierr /= 0 ) &
       CALL errore( 'qedark_f2_3d ',' cannot allocate ctot ', ABS(ierr) )
  ctotA3(:,:,:,:) = 0.0_DP

  ALLOCATE( ctotA1c(numqbins+1, numqbins+1, numqbins+1, num_Er_bins+1) , STAT=ierr )
  IF( ierr /= 0 ) &
       CALL errore( 'qedark_f2_3d ',' cannot allocate ctot ', ABS(ierr) )
  ctotA1c(:,:,:,:) = 0.0_DP

  ALLOCATE( ctotA2c(numqbins+1, numqbins+1, numqbins+1, num_Er_bins+1) , STAT=ierr )
  IF( ierr /= 0 ) &
       CALL errore( 'qedark_f2_3d ',' cannot allocate ctot ', ABS(ierr) )
  ctotA2c(:,:,:,:) = 0.0_DP

  ALLOCATE( ctotA3c(numqbins+1, numqbins+1, numqbins+1, num_Er_bins+1) , STAT=ierr )
  IF( ierr /= 0 ) &
       CALL errore( 'qedark_f2_3d ',' cannot allocate ctot ', ABS(ierr) )
  ctotA3c(:,:,:,:) = 0.0_DP

  ALLOCATE( ctotf1(numqbins+1, numqbins+1, numqbins+1, num_Er_bins+1) , STAT=ierr )
  IF( ierr /= 0 ) &
       CALL errore( 'qedark_f2_3d ',' cannot allocate ctot ', ABS(ierr) )
  ctotf1(:,:,:,:) = 0.0_DP

  ALLOCATE( ctotf2(numqbins+1, numqbins+1, numqbins+1, num_Er_bins+1) , STAT=ierr )
  IF( ierr /= 0 ) &
       CALL errore( 'qedark_f2_3d ',' cannot allocate ctot ', ABS(ierr) )
  ctotf2(:,:,:,:) = 0.0_DP

  ALLOCATE( ctotf3(numqbins+1, numqbins+1, numqbins+1, num_Er_bins+1) , STAT=ierr )
  IF( ierr /= 0 ) &
       CALL errore( 'qedark_f2_3d ',' cannot allocate ctot ', ABS(ierr) )
  ctotf3(:,:,:,:) = 0.0_DP
  er_binsize=ermax_NU/FLOAT(num_er_bins)
  print *, er_binsize,Ry2eV*(binedgesE(2)-binedgesE(1))
  deltaEmax=0.0
  !Turning bzv into the propper prefactor
  bzv=bzv*RAU2NU**3/(2*pi*dq**3*er_binsize)
  IF (restartmode) THEN  

     open(59,file='W_3d.dat')
     DO iqx=1, numqbins
            DO iqy=1, numqbins
                DO iqz=1, numqbins
                        DO iE=1, num_er_bins
READ(59,*), ctot(iqx, iqy, iqz,iE),ctotA1(iqx, iqy, iqz, iE),ctotA2(iqx, iqy, iqz, iE),ctotA3(iqx, iqy, iqz,iE)&
,ctotA1c(iqx, iqy, iqz,iE),ctotA2c(iqx, iqy, iqz, iE),ctotA3c(iqx, iqy, iqz, iE)&
,ctotf1(iqx, iqy, iqz,iE),ctotf2(iqx, iqy, iqz, iE),ctotf3(iqx, iqy, iqz, iE)
                        ENDDO
                ENDDO
        ENDDO
     ENDDO
     close(59)
     DO iqx=1, numqbins
            DO iqy=1, numqbins
                DO iqz=1, numqbins
qnorm=dsqrt((binedgesQ(iqx)+binedgesQ(iqx+1))**2+(binedgesQ(iqy)+binedgesQ(iqy+1))**2+(binedgesQ(iqz)+binedgesQ(iqz+1))**2)/2.0
                        DO iE=1, num_er_bins
ctot(iqx, iqy, iqz, iE) = ctot(iqx, iqy,iqz, iE) /(bzv*qnorm)
ctotA1(iqx, iqy, iqz, iE) = ctotA1(iqx,iqy, iqz, iE) /(bzv*qnorm)
ctotA2(iqx, iqy, iqz, iE) = ctotA2(iqx,iqy, iqz, iE) /(bzv*qnorm)
ctotA3(iqx, iqy, iqz, iE) = ctotA3(iqx,iqy, iqz, iE) /(bzv*qnorm)
ctotA1c(iqx, iqy, iqz, iE) =ctotA1c(iqx, iqy, iqz, iE) /(bzv*qnorm)
ctotA2c(iqx, iqy, iqz, iE) =ctotA2c(iqx, iqy, iqz, iE) /(bzv*qnorm)
ctotA3c(iqx, iqy, iqz, iE) =ctotA3c(iqx, iqy, iqz, iE) /(bzv*qnorm)
ctotf1(iqx, iqy, iqz, iE) = ctotf1(iqx,iqy, iqz, iE)/(bzv*qnorm)
ctotf2(iqx, iqy, iqz, iE) = ctotf2(iqx,iqy, iqz, iE)/(bzv*qnorm)
ctotf3(iqx, iqy, iqz, iE) = ctotf3(iqx,iqy, iqz, iE)/(bzv*qnorm)
                        ENDDO
                ENDDO
        ENDDO
     ENDDO

  ENDIF
  iparastart=3
  IF (restartmode .eqv. .false.) iparastart=1
  ! Initialize tables
  CALL bdotb(bb)
  CALL all_igk(nks, npwx, alligk)
  CALL create_dk_table(nks, xk, dk)
  !Computing the energy binsize. Below will be the first time in the code it is
  !used.
  print *, " "
  print *, " "
  print *, " "
  print *, SIZE(ctot)
  IF (runff) THEN
 
     !  WRITE(18,*) "            ik          iq         i          i'          |f_\{i-->i'}(ik,iq)|^2"     


     ik1init=1
     ik2init=1
     DO ipara=iparastart,CEILING(FLOAT(nksf)/omp_get_max_threads())
     ! Loop over k-vector index of outer wavefunction evcouter
     !$omp parallel &
     !$omp private(ig1, ig2, iband1, iband2, q, qnorm, deltaE, iE, &
     !$omp iqx, iqy, iqz, &
     !$omp ik1,ik2,ik3,evcouter, evc, gsi, f1, f2, fvecs, f1temp, &
     !$omp if1, if2, f2vec, kgvec, nc, fvec,fcross)
     !$omp do
     DO ik3=ik2init, omp_get_max_threads()
        ik2=ik3+(ipara-1)*omp_get_max_threads()
        IF (ik2>nksf) CYCLE
        print *, "Iterating... @ ik2=", ik2, "from thread", omp_get_thread_num() 
        !print *, "Iterating... @ ik2=", ik2, "from thread"
	FLUSH(6)
        ! Load wavefunctions from file         
        CALL get_buffer (evcouter, nwordwfc, iunwfc, ik2)                  
	!print *, "Evcouter loaded", SIZE(evcouter)
        !FLUSH(6)
        DO ik1=ik1init, nksf
           print *, "Iterating... @ ik1=", ik1, "from thread",omp_get_thread_num() 
           FLUSH(6)
           ! Load wavefunctions from file
           CALL get_buffer (evc, nwordwfc, iunwfc, ik1)  
           !print *, "Evc loaded", SIZE(evc)
	  ! FLUSH(6)

           ! Update G-vector table for current (ik1, ik2)
           CALL which_sums(ik1, ik2, alligk, tol, gsi)
           
           ! weight of current point in q-vector space
           wq = 0.5*(wk(ik2) + wk(ik1))                   
           ! Loop over band index of outer wavefunction
           DO iband2=icondbottom, icondtop

           ! Loop over band index of inner wavefunction
              DO iband1=ivalbottom, ivaltop
                
                 deltaE = et(iband2, ik2) - et(iband1, ik1)
           !      IF (deltaE*Ry2eV>70) THEN
           !             print *,deltaE*Ry2eV
           !             FLUSH(6)
           !      ENDIF 
                 iE = find_bin(num_er_bins, binedgesE, deltaE)
		 !Find the maximum deltaE to be able to see how it compares to Ermax
		! IF (deltaE>deltaEmax) THEN
		!	deltaEmax=deltaE
		! ENDIF
                 ! Loop over G-vector index of outer wavefunction
                 DO ig2=1, ngk(ik2) 

                    ! Make sure that G-vector exists
                    IF (alligk(ig2,ik2) < 1) CYCLE                          
                    
                    ! q in cartesian coordinates and convert to RAU multiplying by tpiba
                    q(:) = tpiba * RAU2NU * (xk(:,ik2) - xk(:,ik1) + g(:, alligk(ig2,ik2)) )
                    
                    ! Calculate |q|. This is not needed in this case as the binning happens for each component of q.
                    !qnorm = dsqrt(sum(q(:)**2))
                    

                    iqx = find_bin(numqbins, binedgesQ, q(1))
                    iqy = find_bin(numqbins, binedgesQ, q(2))
                    iqz = find_bin(numqbins, binedgesQ, q(3))
    
                    
                    ! Initialize f1, fvec, and nc for this (ik1, ik2, ig2)
                    f1(:)=0.0
                    fvec(:,:)=0.0
                    !nc=0.0
                    ! Sum over G to calculate formfactor
                    DO ig1=1, ngk(ik1)
                       ! Make sure that G-vector exists
                       IF (alligk(ig1,ik1) < 1) CYCLE
                       !nc=nc+CONJG(evc(ig1,iband1))*evc(ig1,iband1)
                       IF (gsi(ig2,ig1) < 1) CYCLE

                       !Add together the norms

                       IF (noncolin .eqv. .false.) THEN
                          f1temp =  CONJG( evcouter(gsi(ig2,ig1) , iband2) ) * evc(ig1, iband1)

                       ELSE
                          f1temp(1) =  CONJG( evcouter(gsi(ig2,ig1), iband2) ) * evc(ig1, iband1)            ! u-->u
                          f1temp(2) =  CONJG( evcouter(gsi(ig2,ig1)+npwx, iband2) ) * evc(ig1, iband1)       ! u-->d
                          f1temp(3) =  CONJG( evcouter(gsi(ig2,ig1), iband2) ) * evc(ig1+npwx, iband1)       ! d-->u
                          f1temp(4) =  CONJG( evcouter(gsi(ig2,ig1)+npwx, iband2) ) * evc(ig1+npwx, iband1)  ! d-->d
                       ENDIF
                       ! The vector appearing in the expression for f-vector.
                       f1(:)=f1(:)+f1temp(:)
                       kgvec=(xk(:,ik1) + g(:,alligk(ig1,ik1)))

                       DO if1=1, 4
                          DO if2=1, 3
                             fvec(if1,if2)=fvec(if1,if2)+f1temp(if1)* kgvec(if2)
                          ENDDO

                       ENDDO

                    ENDDO !G-vector sum for formfactor
                    ! Take norm squared and sum over spins (when non-collinear)
                    IF(noncolin .eqv. .false.) THEN
                       f2= CONJG(f1(1))*f1(1)
                       fvecs(:)=fvec(1,:)
                    ELSE
                       f1(1)= SUM(f1(:))
                       DO if2=1, 3
                          fvecs(if2)=SUM(fvec(:,if2))
                       ENDDO
                       f2   = DBLE( CONJG(f1(1))*f1(1) )
                    ENDIF
		    !Multiply with suitable prefactors
                    fvecs(:)=-fvecs(:)*tpiba*RAU2NU/me_NU
                    f2vec(:)=CONJG(fvecs(:))*f1(1)
		    !print *, "Updating ctot"
		    !FLUSH(6)
                    fcross=q(1)*2*(REAL(fvecs(2))*AIMAG(fvecs(3))-REAL(fvecs(3))*AIMAG(fvecs(2)))&
+ q(2)*2*(REAL(fvecs(3))*AIMAG(fvecs(1))-REAL(fvecs(1))*AIMAG(fvecs(3))) &
+ q(3)*2*(REAL(fvecs(1))*AIMAG(fvecs(2))-REAL(fvecs(2))*AIMAG(fvecs(1)))

                    ! Add contribution to corresponding bin
		    !ctot is the first crystal response, similar to the one found by Essig et al.
		    !$OMP ATOMIC UPDATE
		    ctot(iqx, iqy, iqz, iE) = ctot(iqx, iqy, iqz, iE) + f2 * wk(ik1) * wk(ik2)
		    !$OMP END ATOMIC
		    !ctotA1-ctotA3 and ctotA1c-ctotA3c are the real parts and immaginary parts of respectively of the second crystal response.
		    !$OMP ATOMIC  UPDATE
		    ctotA1(iqx, iqy, iqz, iE) = ctotA1(iqx, iqy, iqz, iE) + REAL(f2vec(1)) * wk(ik1) * wk(ik2)
		    !$OMP END ATOMIC
		    !$OMP ATOMIC UPDATE
		    ctotA2(iqx, iqy, iqz, iE) = ctotA2(iqx, iqy, iqz, iE) + REAL(f2vec(2)) * wk(ik1) * wk(ik2)
		    !$OMP END ATOMIC
		    !$OMP ATOMIC UPDATE
	            ctotA3(iqx, iqy, iqz, iE) = ctotA3(iqx, iqy, iqz, iE) + REAL(f2vec(3)) * wk(ik1) * wk(ik2)
		    !$OMP END ATOMIC
		    !$OMP ATOMIC UPDATE
		    ctotA1c(iqx, iqy, iqz, iE) = ctotA1c(iqx, iqy, iqz, iE) + AIMAG(f2vec(1)) * wk(ik1) * wk(ik2)
		    !$OMP END ATOMIC
		    !$OMP ATOMIC UPDATE
		    ctotA2c(iqx, iqy, iqz, iE) = ctotA2c(iqx, iqy, iqz, iE) + AIMAG(f2vec(2)) * wk(ik1) * wk(ik2)
		    !$OMP END ATOMIC
		    !$OMP ATOMIC UPDATE
		    ctotA3c(iqx, iqy, iqz, iE) = ctotA3c(iqx, iqy, iqz, iE) + AIMAG(f2vec(3)) * wk(ik1) * wk(ik2)
		    !$OMP END ATOMIC
		    !ctotf1 is the 3. crystal response
		    !$OMP ATOMIC UPDATE
		    ctotf1(iqx, iqy, iqz, iE) = ctotf1(iqx, iqy, iqz, iE) +(REAL(fvecs(1))**2 + AIMAG(fvecs(1))**2 + REAL(fvecs(2))**2 &
+ AIMAG(fvecs(2))**2 + REAL(fvecs(3))**2 + AIMAG(fvecs(3))**2) * wk(ik1) * wk(ik2)
		    !$OMP END ATOMIC
		    !ctotf2 is the 4. crystal response
		    !$OMP ATOMIC UPDATE
		    ctotf2(iqx, iqy, iqz, iE) = ctotf2(iqx, iqy, iqz, iE) + ((REAL(fvecs(1))*q(1)+REAL(fvecs(2))*q(2)+REAL(fvecs(3))*q(3))**2 &
+ (AIMAG(fvecs(1))*q(1)+AIMAG(fvecs(2))*q(2)+AIMAG(fvecs(3))*q(3))**2 ) * wk(ik1) * wk(ik2)/me_NU**2
		    !$OMP END ATOMIC
		    !ctotf3 is the 5. and last crystal response.
                    !$OMP ATOMIC UPDATE
		    ctotf3(iqx, iqy, iqz, iE) = ctotf3(iqx, iqy, iqz, iE)+fcross* wk(ik1) * wk(ik2)/me_NU 
                    !$OMP END ATOMIC
		 ENDDO
                !print *, "G2-loop done"
	        !FLUSH(6) 
              ENDDO
              
           ENDDO !iband2 loop
	   !print *,"Done with first parallel-run"
	   !FLUSH(6) 
        ENDDO
     ENDDO 
     !$omp end do
     !$omp end parallel
!     print *, "Maximum deltaE=",deltaEmax*Ry2eV, "eV."   
          
     ! Multiply by prefactors
     !ctot(:,:) =  ctot(:,:)/(er_binsize*dq)
     !Computing the energy binsize. Below will be the first time in the code it is used.
 !    er_binsize=ermax_NU/FLOAT(num_er_bins)
  !   print *, er_binsize,Ry2eV*(binedgesE(2)-binedgesE(1))
  !   FLUSH(6)
     ! Print to file
     print *, "Writing to file for ipara=", ipara, "out of ",CEILING(FLOAT(nksf)/omp_get_max_threads())
     FLUSH(6)
     open(58,file='q_3d.dat',status="replace")
     DO iqx=1, numqbins
 	    DO iqy=1, numqbins
        	DO iqz=1, numqbins
           		DO iE=1, num_er_bins
!				IF (ABS(ctot(iqx, iqy, iqz, iE))+ABS(ctotA1(iqx, iqy, iqz, iE))+ABS(ctotA2(iqx, iqy, iqz, iE))&
!					+ABS(ctotA3(iqx, iqy, iqz, iE))+ABS(ctotA1c(iqx, iqy, iqz, iE))+ABS(ctotA2c(iqx, iqy, iqz, iE))&
!					+ABS(ctotA3c(iqx, iqy, iqz, iE))+ABS(ctotf1(iqx, iqy, iqz, iE))+ABS(ctotf2(iqx, iqy, iqz, iE))&
!					+ABS(ctotf3(iqx, iqy, iqz, iE))>0.0) THEN
              				WRITE(58,*), (binedgesQ(iqx)+binedgesQ(iqx+1))/2.0, (binedgesQ(iqy)+binedgesQ(iqy+1))/2.0 &
						,(binedgesQ(iqz)+binedgesQ(iqz+1))/2.0,(binedgesE(iE)+binedgesE(iE+1))*Ry2eV/2.0
!           			ENDIF
			ENDDO
        	ENDDO
     	ENDDO
     ENDDO
     close(58)
     !Turning bzv into the propper prefactor
     !bzv=bzv*RAU2NU**3/(2*pi*dq**3*er_binsize)
     open(59,file='W_3d.dat',status="replace")
     DO iqx=1, numqbins
            DO iqy=1, numqbins
                DO iqz=1, numqbins
			qnorm=dsqrt((binedgesQ(iqx)+binedgesQ(iqx+1))**2+(binedgesQ(iqy)+binedgesQ(iqy+1))**2+(binedgesQ(iqz)+binedgesQ(iqz+1))**2)/2.0
                        DO iE=1, num_er_bins
					WRITE(59,*), ctot(iqx, iqy, iqz, iE)*bzv*qnorm,ctotA1(iqx, iqy, iqz, iE)*bzv*qnorm,&
ctotA2(iqx, iqy, iqz, iE)*bzv*qnorm,ctotA3(iqx, iqy, iqz, iE)*bzv*qnorm&
						,ctotA1c(iqx, iqy, iqz, iE)*bzv*qnorm,ctotA2c(iqx, iqy, iqz, iE)*bzv*qnorm,&
ctotA3c(iqx, iqy, iqz, iE)*bzv*qnorm,ctotf1(iqx, iqy, iqz, iE)*bzv*qnorm&
,ctotf2(iqx, iqy, iqz, iE)*bzv*qnorm,ctotf3(iqx, iqy, iqz, iE)*bzv*qnorm
                        ENDDO
                ENDDO
        ENDDO
     ENDDO
     close(59)
  print *, "Done writing to file"
  FLUSH(6)
  
  ENDDO
  ENDIF  ! runff
  
  DEALLOCATE (evcouter)     
  
  
  print *, " " 
  
  
  CALL stop_clock( ' qedark_f2_3d ' )
  
END SUBROUTINE qedark_f2_3d
