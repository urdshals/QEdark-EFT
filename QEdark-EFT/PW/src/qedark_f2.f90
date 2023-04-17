!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Adrian Soto
!  06-10-2015
!  Stony Brook University
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Modified by Einar Urdshals
!  2022
!  Chalmers University
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This file is part of the code QEdark v1.1.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!
SUBROUTINE qedark_f2( restartmode, &
     nksf, numval, numcond, &
     vearth_SI, vesc_SI, v0_SI, deltav_SI, &
     Er_bin_type, num_er_bins, ermax_NU, er_binsize, &
     numqbins, num_of_nodes, node, dq, &
     do_scissor_correction, scissorgap)
  !
  ! Main driver to evaluate and print the spherically form factor
  ! squared as a function of |q| and E.
  ! 
  ! 
  !
  !
  USE constants, ONLY: rytoev
  USE wavefunctions,           ONLY: evc     ! For collinear, evc(npwx, nbnd) [look at allocate_wfc.f90]
  USE kinds,                          ONLY: DP
  USE wvfct,                          ONLY: nbnd, npwx, et,g2kin  !, btype
  USE gvecw,                          ONLY: ecutwfc 
  USE klist,                          ONLY: nks, ngk, wk, xk, nelec
  USE lsda_mod,                       ONLY: nspin
  USE io_files,                       ONLY: nwordwfc, iunwfc
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
  INTEGER :: ik1init
  LOGICAL :: do_scissor_correction != .false.                      ! 
  REAL(DP) :: scissorgap                                          ! Band gap (in eV) for scissor operator corrected band energies 

  REAL(DP) :: RAU2NU

  REAL(DP) :: vearth_SI, vesc_SI, v0_SI, deltav_SI                !Placeholders
  
  REAL(DP) :: me_NU = 0.510998910E6_DP                             ! Electron mass in NU 
  
  REAL(DP) :: bb(6)                                               ! Dot products of reciprocal lattice basis vectors
  REAL(DP) :: bzv                                                 ! 1BZ volume
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
  REAL(DP) :: Erange(2)

  REAL(DP) :: kgvec(3)


  REAL(DP), ALLOCATABLE :: fp(:,:,:,:)
  REAL(DP), ALLOCATABLE :: f(:,:,:,:)  
  REAL(DP)  :: ftemp

  INTEGER :: ik1                                            ! Indices for k-point loops 
  INTEGER :: ig1                                            ! Indices for G-vector loops
  INTEGER :: iband1                                      ! Indices for band loops 
  INTEGER, ALLOCATABLE :: alligk(:,:)                             ! Here we load from file the igk(:) for all k-points, i.e. alligk(ig, ik) = igk(ig)
  REAL(DP) :: dk(3,nks,nks)                                       ! Table storing all k1-k2 values  TODO: this table is antisymmetric and can be reduced
  INTEGER :: nksf                                                 ! Number of k-points for formfactor calculation

  INTEGER :: numval, numcond                                      ! Number of occupied and unoccupied Kohn-Sham orbitals. TO BE SET BY USER!
  INTEGER :: numvaltot, numcondtot                                ! Total # bands in DFT run. numvaltot==nelec/2 and numcondtot=nbnd-numvaltot 
  INTEGER :: ivalbottom, ivaltop                                  ! Minimum and maximum values for valence band index
  INTEGER :: icondbottom, icondtop                                ! Minimum and maximum values for conduction index
  COMPLEX(DP) :: evctemp
  INTEGER :: ierr                                                 ! Error index

  REAL(DP) :: tol = 1.0E-6                                        ! Tolerance for G-vector distances


  INTEGER :: numqbins
  INTEGER :: num_of_nodes
  INTEGER :: node
  INTEGER :: ikGx                                       
  INTEGER :: ikGy
  INTEGER :: ikGz                                      
  REAL(DP) :: dkG                                                  ! size of k+G bin
  REAL(DP), ALLOCATABLE :: binedgeskG(:)
  REAL(DP) :: dq
  !DEBUGGING VARIABLES
  LOGICAL :: runff                                                ! For debugging purposes: set to false to skip the form factor evaluation                       


  runff = .true. ! Set to false to not skip ff calculation --for debugging purposes


  CALL start_clock( ' qedark_f2 ')

  print *, "           -------             "
  print *, " calculation_mode == f2"
  PRINT *, " num_of_nodes, node = ", num_of_nodes, node
  print *, " Calculating formfactor squared binned in E and q"
  print *, "           -------             "
  print *, "This version contains f1temp"

!  IF (nspin .ne. 1) THEN
!     CALL errore ('qedark_f2', 'Form factor calculation works only for spin-unpolarized systems!', 1)
!  ENDIF


  IF ( nksf > nks .or. nksf < 0 ) &
       CALL errore( 'qedark_f2 ',' nksf has a non-allowed value. Check input. ', ABS(ierr) )


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
  print *, 'ivalbottom: ', ivalbottom, 'ivaltop: ',ivaltop
  IF( numval>numvaltot .or. numcond>numcondtot .or. numval<1 .or. numcond<1 ) &
       CALL errore( 'qedark_f2 ',' Check numval and numcond values in input ', ABS(ierr) )

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





  ALLOCATE ( alligk(npwx, nks) , STAT=ierr )
  IF( ierr /= 0 ) &
       CALL errore( 'qedark_f2',' error allocating alligk ', ABS(ierr) )
  alligk(:,:)=0.0_DP ! when alligk=0, ig doesn't correspond to a G-vector in the set and should be disregarded 

  



  WRITE(*,*), " "


  IF (num_er_bins > 999) CALL errore( 'qedark_f2','Number of energy recoil bins cannot exceed 999 ', ABS(ierr) )
  IF (numqbins > 999) CALL errore( 'qedark_f2','Number of momentum transfer bins cannot exceed 999 ', ABS(ierr) )
  

  ermax_RAU = ermax_NU / Ry2eV 
  RAU2NU = me_NU/137.0
  dkG=2*dsqrt(2*ecutwfc*Ry2eV*me_NU)/numqbins
!!!!!!!!!  CALL print_DM_data(vesc_kmps, vearth_kmps, v0_kmps, mx_NU, ermax_NU)


  print *, " "
  IF (do_scissor_correction) THEN
     CALL scissor(numvaltot, numcondtot, scissorgap, et)
  ENDIF
  FLUSH(6)

!Loops to find minimum and maximum valence band energy needed for the binning
  Erange(:)=et(ivalbottom,1)
    DO iband1=ivalbottom, ivaltop
       DO ik1=ik1init, nksf
          IF (et(iband1,ik1)<Erange(1)) THEN
             Erange(1)=et(iband1,ik1)
          ENDIF
          IF (et(iband1,ik1)>Erange(2)) THEN
             Erange(2)=et(iband1,ik1)
         ENDIF
      ENDDO
   ENDDO
   Erange(2)=Erange(2)+1/Ry2eV
   Erange(1)=Erange(1)-1/Ry2eV
   print *, "Erange found ",Erange
   FLUSH(6)
   er_binsize=Ry2eV*(Erange(2)-Erange(1))/num_er_bins
   !Stuff for binning integrated form factors. 
   !Needs to go after ermax_RAU is calculated
  IF (num_er_bins > 0) THEN
     
     ALLOCATE (binedgese(num_er_bins+2), STAT=ierr )
     IF( ierr /= 0 ) &
          CALL errore( 'qedark_f2',' error allocating binedgesE ', ABS(ierr) )
     
     CALL create_bins(1, Erange(1), Erange(2), &
          num_er_bins, (Erange(2)-Erange(1))/num_er_bins, binedgese)
     WRITE (*,*), " "
     WRITE (*,*) "Creating E bins for formfactor sum ..."
     WRITE (*,*), "bintype=", er_bin_type
     WRITE (*,*), "Energy bin size=", er_binsize, "eV"
  ENDIF

FLUSH(6)

  IF (num_er_bins > 0) THEN
     ALLOCATE (binedgeskG(numqbins+2), STAT=ierr )
     IF( ierr /= 0 ) &
          CALL errore( 'qedark_f2 ',' error allocating binedgesq ', ABS(ierr) )

          
     CALL create_bins(er_bin_type, -dkG*numqbins/2.0, dkG*numqbins/2.0, &
          numqbins, dkG, binedgeskG)

     WRITE (*,*), " "
     WRITE (*,*) "Creating sintheta bins for formfactor sum ..."
     WRITE (*,*), "bintype=", er_bin_type
     WRITE (*,*), "Momentum transfer bin size (tpiba)=", 1.0_DP/numqbins
  ENDIF
FLUSH(6)


  ALLOCATE( f(numqbins,numqbins,CEILING(numqbins/2.0), num_Er_bins) , STAT=ierr )
  IF( ierr /= 0 ) &
       CALL errore( 'qedark_f2 ',' cannot allocate f ', ABS(ierr) )
  f(:,:,:,:) = 0.0_DP
  ALLOCATE( fp(numqbins,numqbins,CEILING(numqbins/2.0), num_Er_bins) , STAT=ierr )
  IF( ierr /= 0 ) &
       CALL errore( 'qedark_f2 ',' cannot allocate fp ', ABS(ierr) )
  !Computing the energy binsize. Below will be the first time in the code it is
  !used.
  ! Initialize tables
  CALL bdotb(bb)
  CALL all_igk(nks, npwx, alligk)
  CALL create_dk_table(nks, xk, dk)
  print *, " "
  print *, " "
  print *, num_er_bins,er_binsize
  print *, " "
  IF (runff) THEN

     !  WRITE(18,*) "            ik          iq         i          i'          |f_\{i-->i'}(ik,iq)|^2"     

     ik1init=1

        !$omp parallel &
        !$omp private(ig1, iband1, deltaE, iE, &
        !$omp ikGx,ikGy,ikGz, &
        !$omp ik1, evc, evctemp, &
        !$omp kgvec,fp,ftemp) shared(f)
        fp(:,:,:,:) = 0.0_DP
        !$omp do
        DO ik1=ik1init, nksf
           print *, "Iterating... @ ik1=", ik1, "from thread",omp_get_thread_num()
           FLUSH(6)
           ! Load wavefunctions from file
           CALL get_buffer(evc, nwordwfc, iunwfc, ik1)  
           ! Loop over band index of inner wavefunction
              DO iband1=ivalbottom, ivaltop
                 deltaE = et(iband1, ik1)
                 iE = find_bin(num_er_bins, binedgesE, deltaE)
                 ! Loop over G-vector
                 DO ig1=1, ngk(ik1) 
                    ! Make sure that G-vector exists
                    IF (alligk(ig1,ik1) < 1) CYCLE
                    kgvec(:) = tpiba * RAU2NU * (xk(:,ik1) + g(:, alligk(ig1,ik1)))
                    ikGx = find_bin(numqbins, binedgeskG, (kgvec(1)))
                    ikGy = find_bin(numqbins, binedgeskG, (kgvec(2)))
                    ikGz = find_bin(numqbins, binedgeskG, -1*ABS(kgvec(3)))
                    evctemp=evc(ig1, iband1)
                    ftemp=(REAL(evctemp)**2+AIMAG(evctemp)**2)*wk(ik1)
                    fp(ikGx,ikGy,ikGz, iE) = fp(ikGx,ikGy,ikGz, iE) + ftemp
                 ENDDO
              ENDDO
        print *, ik1, " done"
        FLUSH(6)
        ENDDO
     !$omp end do
     !$OMP CRITICAL(Add_to_f)
     f=f+fp
     !$OMP END CRITICAL(Add_to_f)
     !$omp end parallel
! Print to file
print *, "fsum: ", SUM(f)
print *, "Writing to file..."
FLUSH(6)
     open(58,file='kG.dat',status="replace")
     DO ikGx=1, numqbins
WRITE(58,*),(binedgeskG(ikGx)+binedgeskG(ikGx+1))/2.0
     ENDDO
     close(58)

open(59,file='E.dat',status="replace")
        DO iE=1, num_er_bins
WRITE(59,*),(binedgesE(iE)+binedgesE(iE+1))*Ry2eV/2.0
     ENDDO
     close(59)
     
     open(60,file='W.dat',status="replace")
     DO ikGx=1, numqbins
     DO ikGy=1, numqbins
     DO ikGz=1, CEILING(numqbins/2.0)
        DO iE=1, num_er_bins
WRITE(60,*), f(ikGx,ikGy,ikGz, iE)/(4*(binedgeskG(ikGx+1)-binedgeskG(ikGx))* &
(binedgeskG(ikGy+1)-binedgeskG(ikGy))*(binedgeskG(ikGz+1)-binedgeskG(ikGz))* &
(binedgesE(iE+1)-binedgesE(iE))*Ry2eV) !The factor 4 comes from wk summing to 2 and dkG_z being 2 times dkG (or equivalently that the sampling in each bin in z happens twice; once for positive z and once for negative z).
        ENDDO
     ENDDO
     ENDDO
     ENDDO
     close(60)
  print *, "Done writing to file"
  FLUSH(6)



  ENDIF  ! runff
  
  
  
  print *, " " 
  
  CALL stop_clock( ' qedark_f2 ' )
  
END SUBROUTINE qedark_f2
