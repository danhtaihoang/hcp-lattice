!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!     HOANG Danh Tai - Laboratoire de Physique Théorique et Modélitab_spin_ation 
!     UMR 8089 CNRS-Université de Cergy-Pontoise
!     2, Avenue Adolphe Chauvin, 95032 Cergy-Pontoise Cedex, France
!----------------------------------------------------------------------------------------------------!
!     PROGRAMME: MONTE CARLO TRANSPORT OF SPINS - HCP
!     Bat dau tu ngay 8/9/2010: Tinh E, M, 
!     04.12.2011: bo sung subroutine cacul time relaxation
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!

      PROGRAM main
      IMPLICIT NONE
!*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*==*==*==   
      ! Déclarations 
!*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*==*==*== 

      CHARACTER (256)       Ligne11,Ligne121,Ligne122,Ligne13,Ligne14
      CHARACTER (LEN=150):: CONF,RELOAD,USE_RELAX
      CHARACTER (LEN=15) :: tmp
      CHARACTER (LEN=50) :: name
      CHARACTER (LEN=3)  :: spinAt

      INTEGER (KIND=4),PARAMETER :: n_layer=1
      INTEGER (KIND=4),PARAMETER :: n_atom_p_m=2
      INTEGER (KIND=4),PARAMETER :: n_value_Mz=5

      REAL    (KIND=8),PARAMETER :: nul=0.  

      INTEGER (KIND=4) :: tab_spin_aveconf,tab_spin_avedata,bugvar,case_histo,case_histo_O,case_histo_P
      INTEGER (KIND=4) :: natx,naty,natz,natx_p,naty_p,natz_p,n_atom,n_maille
      INTEGER (KIND=4) :: n_e_itin,n_T
      INTEGER (KIND=4) :: i_equi_reseau_1,n_equi_reseau_1,i_equi_reseau_2,n_equi_reseau_2
      INTEGER (KIND=4) :: i_value_Mz,i_average_thermal,n_average_thermal
      INTEGER (KIND=4) :: alp,bet,gam,gam_new,gam_old
      INTEGER (KIND=4) :: number,bcx,bcy,bcz,bcnm,pix,piy,piz,dim_tab_maille_vois
      INTEGER (KIND=4) :: i,j,k,l,ip,im,jp,jm,kp,km,i_T,itme,del,compt,n_motif
      INTEGER (KIND=4) :: i_a,i_e,i_m,i_e_tmp,dprime
      INTEGER (KIND=4) :: test_pauli,i_gauss,n_gauss
      INTEGER (KIND=4) :: i_equi_elec_1,n_equi_elec_1,i_equi_elec_2,n_equi_elec_2
      INTEGER (KIND=4) :: i_elec_transport,n_transport_elec,i_average_transport,n_average_transport
      INTEGER (KIND=4) :: n_traverse

      REAL    (KIND=8) :: aaa,aaa2,aaa3,x,y,z,x0,y0,z0
      REAL    (KIND=8) :: vx,vy,vz,costheta,cosphi,sintheta,sinphi,fieldE,midmaille,denom_gauss,tab_tmp
      REAL    (KIND=8) :: px,py,pz,rdm_vit1,fieldB,N0,Vsph,Vbox
      REAL    (KIND=8) :: kev,kJ,Na,pi,q,vitesse,I0,K0,costhetatmp
      REAL    (KIND=8) :: temp_min,temp_max,EN_tmp1,EN_tmp2,EN_tmp
      REAL    (KIND=8) :: nr1_tot,Rho1,nr1

      REAL    (KIND=8) :: phi,phitmp
      REAL    (KIND=8) :: T,dT,spin_tmp,energy,energy_moy,energy2,energy2_moy,Cv,Ksi

      REAL    (KIND=8) :: prx,pry,prz,d1,d2,d,energy_new,energy_old,rdm_mtrp,tmp_s_lattice
      REAL    (KIND=8) :: energy_e_e_new,energy_e_a_new,dmod_e_a_new,dmod_e_e_new,pvx,pvy,cptpos
      REAL    (KIND=8) :: energy_e_e_old,energy_e_a_old,dmod_e_a_old,dmod_e_e_old,nr2_tot,Rho2,nr2
      REAL    (KIND=8) :: density_new,density_old,Rho,cpt_vit_moy,n_e_layer1,nr3_tot,Rho3,nr3
      REAL    (KIND=8) :: vit_moy_x,vit_moy_y,vit_moy_z,vit_moy,n_e_layer2,comptmc,Rt1,Rt2,R_drift
      REAL    (KIND=8) :: energy_new_tot,energy_old_tot,vx_tmp,vy_tmp,vz_tmp,min_energy,wl
      

      REAL    (KIND=8) :: J11_stmp,J12_stmp,J21_stmp,J22_stmp,J2_stmp,J1v,J2v,J1s,J2s
      REAL    (KIND=8) :: J3_stmp,J31_stmp,J32_stmp,J3v,J3s
      REAL (KIND=8)    :: Tc,time_relax

      INTEGER (KIND=8),DIMENSION(3)  :: clock 

      INTEGER (KIND=4),DIMENSION(:,:)    ,ALLOCATABLE :: tab_maille_e  
      INTEGER (KIND=4),DIMENSION(:,:)    ,ALLOCATABLE :: tab_maille_a
      INTEGER (KIND=4),DIMENSION(:)      ,ALLOCATABLE :: tab_maille_vois_new
      INTEGER (KIND=4),DIMENSION(:)      ,ALLOCATABLE :: tab_maille_vois_old
  
      REAL    (KIND=8),DIMENSION(:,:,:,:),ALLOCATABLE :: tab_spin_a
      REAL    (KIND=8),DIMENSION(:)      ,ALLOCATABLE :: Mz
      REAL    (KIND=8),DIMENSION(:)      ,ALLOCATABLE :: Mz_2
      REAL    (KIND=8),DIMENSION(:)      ,ALLOCATABLE :: Mz_moy
      REAL    (KIND=8),DIMENSION(:)      ,ALLOCATABLE :: Mz_2_moy

      REAL    (KIND=8),DIMENSION(:,:)    ,ALLOCATABLE :: tab_pos_a
      REAL    (KIND=8),DIMENSION(:,:)    ,ALLOCATABLE :: tab_pos_e

      REAL    (KIND=8),DIMENSION(4)                   :: tab_pos_e_tmp
            
      REAL    (KIND=8),DIMENSION(100)                 :: tab_dist_ee_histogram
      REAL    (KIND=8),DIMENSION(-630:630)            :: tab_theta_histogram
      REAL    (KIND=8),DIMENSION(-630:630)            :: tab_phi_histogram
      REAL    (KIND=8),DIMENSION(-510:510,0:630)      :: tab_thetaphi_histogram
      REAL    (KIND=8),DIMENSION(-500:500)            :: tab_vx_histogram
      REAL    (KIND=8),DIMENSION(-500:500)            :: tab_vy_histogram
      REAL    (KIND=8),DIMENSION(-500:500)            :: tab_vz_histogram

   

!------------------------------------------ END DECLARATION ------------------------------------------!

!Création d'un dossier de stockage config_out dans le repertoire courant
      CALL system('rm -r config_out')
      CALL system('mkdir config_out')
      CALL system('rm -r bank_config_out')
      CALL system('mkdir bank_config_out')
      CALL system('rm *.dat*')
      CALL system('rm *.pdb*')

!Ouverture du fichier de sortie average
      OPEN(unit=121,file='average_thermal.dat')
      OPEN(unit=122,file='M_average_thermal.dat')
      OPEN(unit=13,file='average_transport.dat')
      !OPEN(unit=14,file='M_layer_by_layer.dat')
      !OPEN(unit=15,file='n_traverse.dat')
      
      OPEN(unit=200,file='para_relax.dat')
      
      
!Appel de la routine d initialitab_spin_ation du generateur aleatoire
      CALL init_rdm_number()

!Appel routine de lecture des parametres d entres de parameter.in
      CALL read_input_parameter_file()

!Initialitab_spin_ation variable
      natx_p=natx+1
      naty_p=naty+1
      natz_p=natz+1
      dim_tab_maille_vois=1

      dprime=ceiling(max(d1,d2))   
      del=dprime
      dim_tab_maille_vois=((2*del+1)**3)+1

!Rapport de bug
      IF((d1>real(natx)/2.).or.(d1>real(naty)/2.).or.(d2>real(natx)/2.).or.(d2>real(naty)/2.))THEN
                WRITE(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
                WRITE(*,*) 'BUG : Cutoff d1 and d2 need to be <= at an half' 
                WRITE(*,*) 'lenght box'
                WRITE(*,*) 'RUN WAS STOP'
                WRITE(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
                STOP
      ENDIF
        
!Allocation dynamique des tableaux           
      ALLOCATE(tab_spin_a(0:natx_p,0:naty_p,0:natz_p,2))
      ALLOCATE(tab_maille_e(n_maille,n_e_itin+1))
      ALLOCATE(tab_maille_a(n_atom,n_atom_p_m))
      ALLOCATE(tab_maille_vois_new(dim_tab_maille_vois))
      ALLOCATE(tab_maille_vois_old(dim_tab_maille_vois))
      ALLOCATE(tab_pos_e(n_e_itin,4))
      ALLOCATE(tab_pos_a(n_atom,4))
      ALLOCATE(Mz(n_value_Mz))
      ALLOCATE(Mz_2(n_value_Mz))
      ALLOCATE(Mz_moy(n_value_Mz))
      ALLOCATE(Mz_2_moy(n_value_Mz))


!Calcul du pas en temperature
      IF (n_T==1) THEN 
            dT=0.
      ELSE
            dT=(temp_max-temp_min)/real(n_T-1)
      END IF


      CALL load_config_thermal()
!Boucle en temperature
      DO i_T=1,n_T                                            !1::DEBUT i_T::!
            T=temp_min+dT*real(i_T-1)
            WRITE(*,*)'i_T=',i_T

            Mz_moy=0.
            energy_moy=0.
            Mz_2_moy=0.
            energy2_moy=0.

            !IF(RELOAD=='YES')THEN
                  !CALL read_conf3D_thermal()
            !ELSE IF(RELOAD=='NO')THEN
                  !CALL load_config_thermal()
            !ENDIF

            IF (USE_RELAX=='YES') THEN   
                  CALL cacul_time_relax()
            ENDIF

!*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*==*==*==*!      
      !THERMALItab_spin_aTION ET MOYENNE
!*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*==*==*==*!
            CALL equi_reseau_1()

            IF(RELOAD=='NO')THEN
                  CALL average_thermal()
                  CALL tab_spin_ave_conf3D_thermal()
            ENDIF

            CALL write_conf3D_thermal() 

!*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*==*==*==*!      
      ! DEPLACEMENT DES ELECTRONS
!*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*==*==*==*!
            IF(n_e_itin > 0)THEN
                  bugvar=0
                  CALL load_config_transport_random()
                  !CALL test_load_config_transport_random()
                  bugvar=1
                  CALL write_conf3D_transport_ini()
                  CALL average_transport()
                  CALL write_conf3D_transport()
                  !CALL magnetization_layer_by_layer()
            ENDIF
            
            WRITE(*,*)time_relax,n_equi_reseau_2,n_equi_elec_2,n_transport_elec,n_average_transport
      ENDDO  
                                                     !1::FIN i_T::!
!Detab_spin_allocation dynamique des tableaux et fermeture des fichiers de sorties       
      DEALLOCATE(tab_spin_a)
      DEALLOCATE(tab_maille_e)
      DEALLOCATE(tab_maille_a)
      DEALLOCATE(tab_maille_vois_new)
      DEALLOCATE(tab_maille_vois_old)
      DEALLOCATE(tab_pos_e)
      DEALLOCATE(tab_pos_a)
      CLOSE(121)
      CLOSE(122)
      CLOSE(13)
      !CLOSE(14)
      !CLOSE(15)
      CLOSE(200)

!=====================================================================================================!
!*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*==*==*==*! 
!=====================================================================================================!
!=====================================================================================================!
!*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*==*==*==*! 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
 ! SUBROUTINES     ! SUBROUTINES     ! SUBROUTINES     ! SUBROUTINES    ! SUBROUTINES    ! SUBROUTINES 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*==*==*==*!
!=====================================================================================================!
!=====================================================================================================!
!*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*==*==*==*! 
!=====================================================================================================!
      CONTAINS
!=====================================================================================================
!__ SUBROUTINE read_input-parameter_file : 
!              Lecture des parametres d entre de la simulation. ______________________________________
!_____________________________________________________________________________________________________
      SUBROUTINE read_input_parameter_file()
      IMPLICIT NONE

      CHARACTER (LEN=150) :: tamp

      OPEN(11,file='parameter.in')
      
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(A10))')   tamp, USE_RELAX
      READ(11, '(A30,(A10))')   tamp, CONF
      READ(11, '(A30,(A10))')   tamp, RELOAD
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(I8))')    tamp, tab_spin_aveconf
      READ(11, '(A30,(I8))')    tamp, tab_spin_avedata
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(3I5))')   tamp, natx, naty, natz
      READ(11, '(A30,(F7.4))')  tamp, aaa
      READ(11, '(A30,(I8))')    tamp, n_e_itin
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(F7.4))')  tamp, d1
      READ(11, '(A30,(F7.4))')  tamp, d2
      READ(11, '(A30,(F7.4))')  tamp, d
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(F7.4))')  tamp, J1v
      READ(11, '(A30,(F7.4))')  tamp, J1s
      READ(11, '(A30,(F7.4))')  tamp, J2v
      READ(11, '(A30,(F7.4))')  tamp, J2s
      READ(11, '(A30,(F7.4))')  tamp, J3v
      READ(11, '(A30,(F7.4))')  tamp, J3s
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(F7.4))')  tamp, I0
      READ(11, '(A30,(F7.4))')  tamp, K0
      READ(11, '(A30,(F7.4))')  tamp, fieldB
      READ(11, '(A30,(F7.4))')  tamp, fieldE
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(I8))')    tamp, n_T
      READ(11, '(A30,(F7.4))')  tamp, temp_min
      READ(11, '(A30,(F7.4))')  tamp, temp_max
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(I8))')    tamp, n_equi_reseau_1
      READ(11, '(A30,(I8))')    tamp, n_equi_reseau_2
      READ(11, '(A30,(I8))')    tamp, n_average_thermal
      READ(11, '(A50)')         tamp
      READ(11, '(A30,(I8))')    tamp, n_equi_elec_1
      READ(11, '(A30,(I8))')    tamp, n_equi_elec_2
      READ(11, '(A30,(I8))')    tamp, n_transport_elec
      READ(11, '(A30,(I8))')    tamp, n_average_transport
      READ(11, '(A50)')         tamp
      READ(11, '(A50)')         tamp

      CLOSE(11) 

      !Constantes
      kev=(1.3806505/1.60217653)*10.**(-4) !ev/K
      q=1.60217653E-19                     !C
      kJ=1.38065812E-23                    !kg.m**2.s**-2.K**-1
      Na=6.02214179E23                     !mol**-1
      pi=3.141592653589793238462643d0      !pi

      IF (CONF=='HCP') THEN
            n_atom=n_atom_p_m*natx*naty*natz
            n_motif=1
      ENDIF


      
      aaa2=aaa*sqrt(3.)/2.
      aaa3=aaa*sqrt(8./3.)
      n_maille=natx*naty*natz
      Vbox=real(n_maille)
      Vsph=((4./3.)*pi*(d2)**3.)/((4.*pi/3.)*real(n_motif))
      N0=real(n_e_itin)/real(n_maille*n_motif)
      

      END SUBROUTINE read_input_parameter_file
      
!!!=================================================================================================
      SUBROUTINE cacul_time_relax()
      IMPLICIT NONE
  
      Tc=2.265
      
      !WRITE(*,*)'Use effec of time relaxation'
      
      time_relax=1./((abs(1.-T/Tc))**1.28876)
      
      n_transport_elec=nint(time_relax)                               ! n_transport_elec
      n_equi_elec_2=n_transport_elec/20+3                             ! n_equi_elec2
      n_average_transport=50000/(n_equi_elec_2+n_transport_elec)      ! n_average_transport      
      n_equi_reseau_2=50000/n_average_transport+3                       ! n_equi_reseau_2 (3)
          
      WRITE(200,*)T,n_equi_reseau_2,n_equi_elec_2,n_transport_elec,n_average_transport
      
      END SUBROUTINE cacul_time_relax
      
!=====================================================================================================
!__ SUBROUTINE init_rdm_number : 
!              Initialitab_spin_ation du generateur de nombre aletaoire en fonction de
!              l horloge de la machine. _____________________________________________________________
!____________________________________________________________________________________________________
      SUBROUTINE init_rdm_number()

      IMPLICIT NONE

      INTEGER (KIND=4) :: i_time,i

      CALL ITIME(clock)
      i_time=(clock(1)+1)*(clock(2)+1)*(clock(3)+1)
        
      DO i=1,i_time
            CALL RANDOM_NUMBER(tab_tmp)
      ENDDO 

      END SUBROUTINE init_rdm_number
!=====================================================================================================
!__ SUBROUTINE load_config : 
!              Genere la configuration initiale de la structure cristalline a
!              simuler. ______________________________________________________________________________
!_____________________________________________________________________________________________________
      SUBROUTINE load_config_thermal()

      IMPLICIT NONE
      INTEGER (KIND=4) :: i,j,k,l
       
      IF (1==2) THEN
      WRITE (*,*)'config ini spin: GS2'       
      DO k=1,natz
            DO j=1,naty
                  DO i=1,natx
                        DO l=1,n_atom_p_m

                             IF (mod(j,2)==0) THEN
                                   tab_spin_a(i,j,k,1)=1.
                                   tab_spin_a(i,j,k,2)=-1.
                             ELSE
                                  tab_spin_a(i,j,k,1)=-1.
                                  tab_spin_a(i,j,k,2)=1.
                             END IF
                               
                             !tab_spin_a(i,j,k,l)=1.

                              IF (k==1) THEN
                                    tab_spin_a(i,j,0,l)=0.        
                              END IF
                        
                              IF (k==natz) THEN 
                                    tab_spin_a(i,j,natz_p,l)=0.
                              END IF   
                        ENDDO
                  ENDDO
            ENDDO    
      ENDDO

      ELSE
      WRITE (*,*)'config ini spin: GS1' 
      DO k=1,natz
            DO j=1,naty
                  DO i=1,natx
                        DO l=1,n_atom_p_m

                             tab_spin_a(i,j,k,1)=1.
                             tab_spin_a(i,j,k,2)=-1.

                              IF (k==1) THEN
                                     tab_spin_a(i,j,0,l)=0.        
                              END IF
                        
                              IF (k==natz) THEN 
                                    tab_spin_a(i,j,natz_p,l)=0.
                              END IF  

                        ENDDO
                  ENDDO
            ENDDO    
      ENDDO

      END IF

            
      END SUBROUTINE load_config_thermal

!=====================================================================================================
!__ SUBROUTINE sweep loop thermal : 
!              N pastab_spin_age sur tous les atomes de la structure. ________________________________________
!_____________________________________________________________________________________________________
      SUBROUTINE equi_reseau_1()

      IMPLICIT NONE

      DO i_equi_reseau_1=1,n_equi_reseau_1
            CALL equi_atom()
      ENDDO

      END SUBROUTINE equi_reseau_1

!=====================================================================================================
!__ SUBROUTINE sweep loop thermal : 
!              N pastab_spin_age sur tous les atomes de la structure. ________________________________________
!_____________________________________________________________________________________________________
      SUBROUTINE equi_reseau_2()

      IMPLICIT NONE

      DO i_equi_reseau_2=1,0
            CALL equi_atom()
      ENDDO

      END SUBROUTINE equi_reseau_2

!=====================================================================================================
!__ SUBROUTINE Rethermalitab_spin_ation pour transport des electrons ________________________________________
!_____________________________________________________________________________________________________
      SUBROUTINE equi_reseau_3()

      IMPLICIT NONE

      DO i=1,n_equi_reseau_2
            CALL equi_atom()
      ENDDO

      END SUBROUTINE equi_reseau_3

!=====================================================================================================
!__ SUBROUTINE equi_atom : Mz
!              Pastab_spin_age sur tous les atomes de la structure et critere 
!              Metropolis. ___________________________________________________________________________
!_____________________________________________________________________________________________________
      SUBROUTINE equi_atom()

      IMPLICIT NONE    
      INTEGER (KIND=4) :: i,j,k,l
    
      DO k=1,natz
            km=k-1
            kp=k+1

            J11_stmp=J1v
            J12_stmp=J1v
            J2_stmp=J2v
            J21_stmp=J2v
            J22_stmp=J2v
            J3_stmp=J3v
            J31_stmp=J3v
            J32_stmp=J3v
            
            IF (k==1) THEN
                  J11_stmp=J1s
                  J21_stmp=J2s
            ELSE IF (k==natz) THEN
                  J11_stmp=J1s
            END IF

            DO j=1,naty
                  jm=j-1+(1/j)*naty
                  jp=j+1-(j/naty)*naty
                  DO i=1,natx
                        im=i-1+(1/i)*natx
                        ip=i+1-(i/natx)*natx
!11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111
!11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111


                        DO l=1,n_atom_p_m
                                       
                        SELECT CASE(l)
                        
                        CASE(1)

                        EN_tmp1=-tab_spin_a(i,j,k,1) &
                                    *(J11_stmp*(tab_spin_a(ip,j,k,1)+tab_spin_a(i,jp,k,1)&
                                           +tab_spin_a(im,jp,k,1)+tab_spin_a(im,j,k,1)&
                                           +tab_spin_a(i,jm,k,1)+tab_spin_a(ip,jm,k,1))&
                                           
                                    +J22_stmp*(tab_spin_a(i,j,k,2)+tab_spin_a(im,j,k,2)&
                                             +tab_spin_a(i,jm,k,2)&
                                             +tab_spin_a(i,j,km,2)+tab_spin_a(im,j,km,2)&
                                             +tab_spin_a(i,jm,km,2)))


                        CASE(2)
                        EN_tmp1=-tab_spin_a(i,j,k,2) &
                                    *(J11_stmp*(tab_spin_a(ip,j,k,2)+tab_spin_a(i,jp,k,2)&
                                           +tab_spin_a(im,jp,k,2)+tab_spin_a(im,j,k,2)&
                                           +tab_spin_a(i,jm,k,2)+tab_spin_a(ip,jm,k,2))&
                                           
                                    +J22_stmp*(tab_spin_a(i,j,k,1)+tab_spin_a(ip,j,k,1)&
                                             +tab_spin_a(i,jp,k,1)&
                                             +tab_spin_a(i,j,kp,1)+tab_spin_a(ip,j,kp,1)&
                                             +tab_spin_a(i,jp,kp,1)))

                        END SELECT

!11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111
!11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111

                        spin_tmp= -tab_spin_a(i,j,k,l)
                        EN_tmp2 = -EN_tmp1

                        CALL random_number(rdm_mtrp) 
                        IF ((EN_tmp2<=EN_tmp1).OR.(exp(-(EN_tmp2-EN_tmp1)/T) > rdm_mtrp))THEN
                            tab_spin_a(i,j,k,l)=spin_tmp      
                        END IF
                        
                        ENDDO
                  ENDDO
            ENDDO
      ENDDO

      END SUBROUTINE equi_atom


!=====================================================================================================
!__ SUBROUTINE value thermal : 
!              Calcul de energie et aimantation total de la structure. _______________________________
!_____________________________________________________________________________________________________

      SUBROUTINE value_thermal()

      IMPLICIT NONE

      energy=0.
      Mz(:)=0.
 
      DO k=1,natz
            km=k-1
            kp=k+1
      
            J11_stmp=J1v
            J12_stmp=J1v
            J2_stmp=J2v
            J21_stmp=J2v
            J22_stmp=J2v
            J3_stmp=J3v
            J31_stmp=J3v
            J32_stmp=J3v
            
            IF (k==1) THEN
                  J11_stmp=J1s
                  J21_stmp=J2s
            ELSE IF (k==natz) THEN
                  J11_stmp=J1s
            END IF

            DO j=1,naty
                  jm=j-1+(1/j)*naty
                  jp=j+1-(j/naty)*naty
                  DO i=1,natx
                        im=i-1+(1/i)*natx
                        ip=i+1-(i/natx)*natx

!11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111
!11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111


                        DO l=1,n_atom_p_m
                        !Mz(1)=Mz(1)+tab_spin_a(i,j,k,l)
                        
                        SELECT CASE(l)
                        
                        CASE(1)

                        EN_tmp1=-tab_spin_a(i,j,k,1) &
                                    *(J11_stmp*(tab_spin_a(ip,j,k,1)+tab_spin_a(i,jp,k,1)&
                                           +tab_spin_a(im,jp,k,1)+tab_spin_a(im,j,k,1)&
                                           +tab_spin_a(i,jm,k,1)+tab_spin_a(ip,jm,k,1))&
                                           
                                    +J22_stmp*(tab_spin_a(i,j,k,2)+tab_spin_a(im,j,k,2)&
                                             +tab_spin_a(i,jm,k,2)&
                                             +tab_spin_a(i,j,km,2)+tab_spin_a(im,j,km,2)&
                                             +tab_spin_a(i,jm,km,2)))


                        CASE(2)
                        EN_tmp1=-tab_spin_a(i,j,k,2) &
                                    *(J11_stmp*(tab_spin_a(ip,j,k,2)+tab_spin_a(i,jp,k,2)&
                                           +tab_spin_a(im,jp,k,2)+tab_spin_a(im,j,k,2)&
                                           +tab_spin_a(i,jm,k,2)+tab_spin_a(ip,jm,k,2))&
                                           
                                    +J22_stmp*(tab_spin_a(i,j,k,1)+tab_spin_a(ip,j,k,1)&
                                             +tab_spin_a(i,jp,k,1)&
                                             +tab_spin_a(i,j,kp,1)+tab_spin_a(ip,j,kp,1)&
                                             +tab_spin_a(i,jp,kp,1)))

                        END SELECT

!11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111
!11111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111

                        EN_tmp2 = -EN_tmp1
                        spin_tmp= -tab_spin_a(i,j,k,l)

                        CALL random_number(rdm_mtrp)         
                        IF ((EN_tmp2<=EN_tmp1).OR.(exp(-(EN_tmp2-EN_tmp1)/T) > rdm_mtrp))THEN
                            EN_tmp=EN_tmp2
                            tab_spin_a(i,j,k,l)=spin_tmp      
                        ELSE     
                            EN_tmp=EN_tmp1
                        END IF

                        energy=energy+EN_tmp
                                                
                        
                        !Mz(4)=Mz(4)+tab_spin_a(i,j,k,l)

                        ENDDO

                        Mz(1)=Mz(1)+((-1.)**mod(j,2))*tab_spin_a(i,j,k,1) & 
                                   +((-1.)**mod(j+1,2))*tab_spin_a(i,j,k,2)

                        Mz(2)=Mz(2)+((-1.)**mod(i,2))*tab_spin_a(i,j,k,1) & 
                                   +((-1.)**mod(i+1,2))*tab_spin_a(i,j,k,2)

                        Mz(3)=Mz(3)+((-1.)**mod(i+j,2))*tab_spin_a(i,j,k,1) & 
                                   +((-1.)**mod(i+j,2))*tab_spin_a(i,j,k,2)
                        
                        Mz(4)=Mz(4)+((-1.)**mod(j,2))*tab_spin_a(i,j,k,1) & 
                                   +((-1.)**mod(j+1,2))*tab_spin_a(i,j,k,2) &

                                   +((-1.)**mod(i,2))*tab_spin_a(i,j,k,1) & 
                                   +((-1.)**mod(i+1,2))*tab_spin_a(i,j,k,2) &

                                   +((-1.)**mod(i+j,2))*tab_spin_a(i,j,k,1) & 
                                   +((-1.)**mod(i+j,2))*tab_spin_a(i,j,k,2)

                        Mz(5)=Mz(5)+tab_spin_a(i,j,k,1)+tab_spin_a(i,j,k,2)

                  ENDDO
            ENDDO
      ENDDO      
      
      DO i_value_Mz=1,n_value_Mz
            Mz(i_value_Mz)=abs(Mz(i_value_Mz)/real(n_atom))            
      END DO
      !Mz(1)=abs(Mz(1)/real(n_atom))
      
      !WRITE(*,*),'Mz1=', Mz(1)
      energy=energy/(2.*real(n_atom))
      energy2=energy**2.
      Mz_2(1)=Mz(1)**2.
      
      END SUBROUTINE value_thermal

!=====================================================================================================
!__ SUBROUTINE average thermal : 
!              Calcul des moyennes des grandeurs physique. ___________________________________________
!_____________________________________________________________________________________________________
      SUBROUTINE average_thermal()

      IMPLICIT NONE
      
      energy_moy=0.
      energy2_moy=0.
      Mz_moy(:)=0.
      Mz_2_moy(:)=0.

      DO i_average_thermal=1,n_average_thermal
            CALL equi_reseau_2()
            CALL value_thermal()           
            DO i_value_Mz=1,n_value_Mz
                  Mz_moy(i_value_Mz)=Mz_moy(i_value_Mz)+Mz(i_value_Mz)                
            END DO
            energy_moy=energy_moy+energy
            energy2_moy=energy2_moy+energy2
            Mz_2_moy(1)=Mz_2_moy(1)+Mz_2(1)
      ENDDO
      
      DO i_value_Mz=1,n_value_Mz  
            Mz_moy(i_value_Mz)=Mz_moy(i_value_Mz)/real(n_average_thermal)
      END DO
      energy_moy=energy_moy/real(n_average_thermal)
      energy2_moy=energy2_moy/real(n_average_thermal)
      Mz_2_moy(1)=Mz_2_moy(1)/real(n_average_thermal)

      Cv =real(n_atom)*(energy2_moy-energy_moy**2.)/(T**2.)
      Ksi=real(n_atom)*(Mz_2_moy(1)-Mz_moy(1)**2.)/T

      WRITE(Ligne121,*) T,energy_moy,Mz_moy(1),Cv,Ksi
      
      WRITE(Ligne122,*) T,Mz_moy(1),Mz_moy(2),Mz_moy(3),Mz_moy(4),Mz_moy(5)

      WRITE(121,'(a)') trim(Ligne121)  
      WRITE(122,'(a)') trim(Ligne122) 

      END SUBROUTINE average_thermal

!=====================================================================================================
!__ SUBROUTINE lecture conf3D thermal : 
!              Lecture du fichier de structure atome 3D _______________________________
!_____________________________________________________________________________________________________
      SUBROUTINE read_conf3D_thermal()

      IMPLICIT NONE

!------------------------------------------------------------------------------
!Variables locales
      INTEGER (KIND=4) :: i,j,k,l
!------------------------------------------------------------------------------

!Incrémentation du nom de fichier d entre
      number=i_T+10000000
      WRITE(tmp,'(I8)') number

      IF (CONF=='HCP') THEN
            name='_HCP_'//TRIM(tmp)
      ENDIF

!Lecture du fichier
      OPEN(unit=11,file='bank_config/struct_3D_thermalitab_spin_ation'//trim(name)//'.pdb')

      DO k=0,natz_p
            DO j=0,naty_p
                  DO i=0,natx_p
                        DO l=1,n_atom_p_m
                              READ(11,*) tmp_s_lattice
                              tab_spin_a(i,j,k,l)=tmp_s_lattice
                        ENDDO
                  ENDDO
            ENDDO
      ENDDO

      CLOSE(11) 

      END SUBROUTINE read_conf3D_thermal
!=====================================================================================================
!__ SUBROUTINE lecture conf3D thermal: 
!              Lecture du fichier de structure atome 3D _______________________________
!_____________________________________________________________________________________________________
      SUBROUTINE tab_spin_ave_conf3D_thermal()

      IMPLICIT NONE

      INTEGER (KIND=4) :: i,j,k,l

!Incrémentation du nom de fichier de sortie
      number=i_T+10000000
      WRITE(tmp,'(I8)') number

      IF (CONF=='HCP') THEN
            name='_HCP_'//TRIM(tmp)
      ENDIF

!Ecriture du fichier config_ incrémenté
      OPEN(unit=11,file='bank_config_out/struct_3D_thermalitab_spin_ation'//trim(name)//'.pdb')

      DO k=0,natz_p
            DO j=0,naty_p
                  DO i=0,natx_p
                        DO l=1,n_atom_p_m
                              WRITE(Ligne11,*) tab_spin_a(i,j,k,l)
                              WRITE(11,'(a)') trim(Ligne11)
                        ENDDO
                  ENDDO
            ENDDO
      ENDDO

      CLOSE(11) 

      END SUBROUTINE tab_spin_ave_conf3D_thermal

!=====================================================================================================
!__ SUBROUTINE write conf thermal : 
!              Ecriture du fichier de structure 3D. __________________________________________________
!_____________________________________________________________________________________________________
      SUBROUTINE write_conf3D_thermal()

      IMPLICIT NONE

      INTEGER (KIND=4) :: i,j,k
!------------------------------------------------------------------------------

!Incrémentation du nom de fichier de sortie
      number=i_T+10000000
      WRITE(tmp,'(I8)') number

      IF (CONF=='HCP') THEN
            name='_HCP_'//TRIM(tmp)
      ENDIF

!Ecriture du fichier config_ incrémenté
      OPEN(unit=11,file='config_out/conf_3D_thermalitab_spin_ation'//trim(name)//'.pdb')
   
      DO k=1,natz
            DO j=1,naty
                  DO i=1,natx

                        x0=real(i-1)*aaa
                        y0=real(j-1)*aaa2
                        z0=real(k-1)*aaa3
!!!!------------------------------------------------------------------------------------------------      
                        IF (INT(mod(j,2))==1) THEN

                        DO l=1,n_atom_p_m

                              SELECT CASE(l)

                              CASE(1)
                              x=x0 ; y=y0 ; z=z0

                              CASE(2)
                              x=x0+aaa/2. ; y=y0+aaa*sqrt(3.)/6. ; z=z0+aaa*sqrt(2./3.) 
                              END SELECT

                              IF (int(tab_spin_a(i,j,k,l))==1) THEN 
                                    spinAt='Au'
                              ELSE
                                    spinAt='Cu'
                              END IF

                              WRITE(11,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                                    spinAt,x,y,z,nul  

                        ENDDO
!!!!------------------------------------------------------------------------------------------------  
                        ELSE
                        DO l=1,n_atom_p_m                            

                              SELECT CASE(l)

                              CASE(1)
                              x=x0+aaa/2. ; y=y0 ; z=z0                                          

                              CASE(2)
                              x=x0+aaa   ; y=y0+aaa*sqrt(3.)/6. ; z=z0+aaa*sqrt(2./3.) 
                                   
                              END SELECT

                              IF (int(tab_spin_a(i,j,k,l))==1) THEN 
                                    spinAt='Au'
                              ELSE
                                    spinAt='Cu'
                              END IF

                              WRITE(11,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                                    spinAt,x,y,z,nul  

                        ENDDO

                        END IF                       
!!!!------------------------------------------------------------------------------------------------  

                  ENDDO
            ENDDO
      ENDDO

      CLOSE(11) 

      END SUBROUTINE write_conf3D_thermal


!=====================================================================================================!
!*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*==*==*==*! 
!=====================================================================================================!
!=====================================================================================================!
!*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*==*==*==*! 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
      ! champitre TRANSPORT             ! champitre TRANSPORT              ! champitre TRANSPORT      
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
!*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*==*==*==*!
!=====================================================================================================!
!=====================================================================================================!
!*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*===*==*==*==*==*==*==*==*==*==*! 
!=====================================================================================================!

!=====================================================================================================
!__ SUBROUTINE load config transport : 
!              Chargement du tableau atom pour conditions periodique optimise. _______________________
! Xac dinh toa do cua cac atom        
! Goi ten maille, ten atom, tung electron vao
!_____________________________________________________________________________________________________
      SUBROUTINE load_config_transport_random()

      IMPLICIT NONE

      INTEGER (KIND=4) :: i,j,k,l
      !---------------------------------------

      !initialitab_spin_ation du tableau electron dans num_maille de la structure
      tab_maille_e(:,:)=0
      tab_maille_a(:,:)=0

      i_a=0
      bcnm=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      DO k=1,natz
            DO j=1,naty
                  DO i=1,natx
                        bcnm=bcnm+1

                        x=real(i-1)*aaa
                        y=real(j-1)*aaa2
                        z=real(k-1)*aaa3
      
                        IF (INT(mod(j,2))==1) THEN

                        DO l=1,n_atom_p_m
                              i_a=i_a+1
                              tab_maille_a(bcnm,l)=i_a
                              tab_pos_a(i_a,4)=tab_spin_a(i,j,k,l)

                              SELECT CASE(l)

                              CASE(1)
                              tab_pos_a(i_a,1)=x
                              tab_pos_a(i_a,2)=y
                              tab_pos_a(i_a,3)=z                                          

                              CASE(2)
                              tab_pos_a(i_a,1)=x+aaa/2.   
                              tab_pos_a(i_a,2)=y+aaa*sqrt(3.)/6. 
                              tab_pos_a(i_a,3)=z+aaa*sqrt(2./3.) 
                                   
                              END SELECT
                        ENDDO

                        ELSE
                        DO l=1,n_atom_p_m
                              i_a=i_a+1
                              tab_maille_a(bcnm,l)=i_a
                              tab_pos_a(i_a,4)=tab_spin_a(i,j,k,l)

                              SELECT CASE(l)

                              CASE(1)
                              tab_pos_a(i_a,1)=x+aaa/2.
                              tab_pos_a(i_a,2)=y
                              tab_pos_a(i_a,3)=z                                          

                              CASE(2)
                              tab_pos_a(i_a,1)=x+aaa   
                              tab_pos_a(i_a,2)=y+aaa*sqrt(3.)/6. 
                              tab_pos_a(i_a,3)=z+aaa*sqrt(2./3.) 
                                   
                              END SELECT
                        ENDDO

                        END IF

                  ENDDO
            ENDDO
      ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!On place les electrons itinerants aleatoirement dans la structure en respectant 
!le critere d exclusion de Pauli
      DO i_e=1,n_e_itin   

!Test du critere de Pauli
            test_pauli=1
            DO WHILE(test_pauli/=0)

!genere aleatoirement les positions
                  CALL random_number(px)
                  CALL random_number(py)
                  CALL random_number(pz)
                  px=px*aaa*real(natx)
                  py=py*aaa2*real(naty)
                  pz=pz*aaa3*real(natz-0.5)

!Remplistab_spin_age tableau position et spin des electron itinerants
                  tab_pos_e(i_e,1)=px 
                  tab_pos_e(i_e,2)=py 
                  tab_pos_e(i_e,3)=pz 
                  tab_pos_e(i_e,4)=1.

!Appel du remplistab_spin_age tableau maille voisines
                  gam_new=-1
                  CALL bound_cond_transport_old()
                  test_pauli=0

                  !---------------------------------------------------------
                  !Test Pauli e-e  
                  i_m=1
                  DO WHILE((tab_maille_vois_old(i_m)/=0).and.(test_pauli==0))

                        bcnm=tab_maille_vois_old(i_m)

                        itme=1

                        DO WHILE((tab_maille_e(bcnm,itme)/=0).and.(test_pauli==0))

                              i_e_tmp=tab_maille_e(bcnm,itme)

                              IF((tab_maille_e(bcnm,itme) > 0).and.(i_e/=i_e_tmp))THEN
                                    CALL dmod_elec_elec_old()
                                    CALL test_pauli_elec_old()
                              ENDIF

                              itme=itme+1
                        ENDDO
                        i_m=i_m+1     
                  ENDDO

                  !---------------------------------------------------------
                  !Test Pauli e-a
                  i_m=1
                  DO WHILE((tab_maille_vois_old(i_m)/=0).and.(test_pauli==0))
                        bcnm=tab_maille_vois_old(i_m)
                        l=1
                        DO WHILE ((l <= n_atom_p_m).and.(test_pauli==0))
                              i_a=tab_maille_a(bcnm,l)                                         
                              CALL dmod_elec_atom_old()
                              CALL test_pauli_atom_old()
                              l=l+1
                        ENDDO         
                        i_m=i_m+1  
                  ENDDO
!---------------------------------------------------------------------------------
            ENDDO

!!Positionnement de l electron dans le tableau
            alp=int(px/aaa)+1     
            bet=alp+natx*int(py/aaa2)
            bcnm=bet+natx*naty*int(pz/aaa3)

            itme=1
            DO WHILE(tab_maille_e(bcnm,itme)/=0)
                  itme=itme+1                        
            ENDDO
            tab_maille_e(bcnm,itme)=i_e

      ENDDO

      END SUBROUTINE load_config_transport_random


!=====================================================================================================
!=====================================================================================================
!    SUBROUTINE test_load_config_transport_random(): 
!
!              Kiem tra xem viec xac dinh toa do cua cac atom da dung chua?
!              CT nay chi su dung de chay tab_spin_au khi sua lai load_config_transport_random()
!_____________________________________________________________________________________________________
      SUBROUTINE test_load_config_transport_random()

     IMPLICIT NONE

!------------------------------------------------------------------------------
!Variables locales
      INTEGER (KIND=4) :: i_a
!------------------------------------------------------------------------------

      OPEN(112,file='test_conf_3D_transport_ini.pdb')

!Ecriture des atomes
      DO i_a=1,n_atom
            IF (INT(mod(i_a,2))==1) THEN
                  spinAt='Au'
            ELSE
                  spinAt='Cu'
            END IF

            WRITE(112,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
            spinAt,tab_pos_a(i_a,1),tab_pos_a(i_a,2),tab_pos_a(i_a,3),nul
            
      ENDDO

      CLOSE(112) 

      END SUBROUTINE test_load_config_transport_random

!=====================================================================================================
!__ SUBROUTINE reload lattice : 
!              Rechargement du reseau (uniquement) apres rethermalitab_spin_ation de la
!              structure. ____________________________________________________________________________
!_____________________________________________________________________________________________________
      SUBROUTINE reload_lattice()

      IMPLICIT NONE

!------------------------------------------------------------------------------
!Variables locales
      INTEGER (KIND=4) :: i,j,k,l
!------------------------------------------------------------------------------

      bcnm=0
      i_a=0

      DO k=1,natz
            DO j=1,naty
                  DO i=1,natx
                        bcnm=bcnm+1

                        x=real(i-1)*aaa
                        y=real(j-1)*aaa2
                        z=real(k-1)*aaa3
      
                        IF (INT(mod(j,2))==1) THEN

                        DO l=1,n_atom_p_m
                              i_a=i_a+1
                              tab_maille_a(bcnm,l)=i_a
                              tab_pos_a(i_a,4)=tab_spin_a(i,j,k,l)

                              SELECT CASE(l)

                              CASE(1)
                              tab_pos_a(i_a,1)=x
                              tab_pos_a(i_a,2)=y
                              tab_pos_a(i_a,3)=z                                          

                              CASE(2)
                              tab_pos_a(i_a,1)=x+aaa/2.   
                              tab_pos_a(i_a,2)=y+aaa*sqrt(3.)/6. 
                              tab_pos_a(i_a,3)=z+aaa*sqrt(2./3.) 
                                   
                              END SELECT
                        ENDDO

                        ELSE
                        DO l=1,n_atom_p_m
                              i_a=i_a+1
                              tab_maille_a(bcnm,l)=i_a
                              tab_pos_a(i_a,4)=tab_spin_a(i,j,k,l)

                              SELECT CASE(l)

                              CASE(1)
                              tab_pos_a(i_a,1)=x+aaa/2.
                              tab_pos_a(i_a,2)=y
                              tab_pos_a(i_a,3)=z                                          

                              CASE(2)
                              tab_pos_a(i_a,1)=x+aaa   
                              tab_pos_a(i_a,2)=y+aaa*sqrt(3.)/6. 
                              tab_pos_a(i_a,3)=z+aaa*sqrt(2./3.) 
                                   
                              END SELECT
                        ENDDO

                        END IF

                  ENDDO
            ENDDO
      ENDDO

      END SUBROUTINE reload_lattice

!=====================================================================================================
!__ SUBROUTINE write conf3D transport : 
!              Ecriture du fichier de structure atome 3D + electron 3D. ______________________________
!_____________________________________________________________________________________________________
      SUBROUTINE write_conf3D_transport_ini()

      IMPLICIT NONE

!------------------------------------------------------------------------------
!Variables locales
      INTEGER (KIND=4) :: i_a,i_e
!------------------------------------------------------------------------------

!Incrémentation du nom de fichier de sortie
      number=i_T+10000000
      WRITE(tmp,'(I8)') number

      IF (CONF=='HCP') THEN
            name='_HCP_ini_'//TRIM(tmp)
      ENDIF

!Ecriture du fichier config_ incrémenté
      OPEN(unit=11,file='config_out/conf_3D_transport'//trim(name)//'.pdb')

!Ecriture des atomes
      DO i_a=1,n_atom
            IF (tab_pos_a(i_a,4)==1) THEN
                  spinAt='Au'
                  WRITE(11,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                  spinAt,tab_pos_a(i_a,1),tab_pos_a(i_a,2),tab_pos_a(i_a,3),nul
            END IF
            
            IF (tab_pos_a(i_a,4)==-1) THEN
                  spinAt='Cu'
                  WRITE(11,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                  spinAt,tab_pos_a(i_a,1),tab_pos_a(i_a,2),tab_pos_a(i_a,3),nul
            END IF
            
      ENDDO
!Ecriture des electrons
      DO i_e=1,n_e_itin
            spinAt='H'
            WRITE(11,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
            spinAt,tab_pos_e(i_e,1),tab_pos_e(i_e,2),tab_pos_e(i_e,3),nul
      ENDDO

      CLOSE(11) 

      END SUBROUTINE write_conf3D_transport_ini


!=====================================================================================================
!__ SUBROUTINE write conf3D transport : 
!              Ecriture du fichier de structure atome 3D + electron 3D. ______________________________
!_____________________________________________________________________________________________________
      SUBROUTINE write_conf3D_transport()

      IMPLICIT NONE

!------------------------------------------------------------------------------
!Variables locales
      INTEGER (KIND=4) :: i_a,i_e
!------------------------------------------------------------------------------

!Incrémentation du nom de fichier de sortie
      number=i_T+10000000
      WRITE(tmp,'(I8)') number

      IF (CONF=='HCP') THEN
            name='_HCP_'//TRIM(tmp)
      ENDIF

!Ecriture du fichier config_ incrémenté
      OPEN(unit=11,file='config_out/conf_3D_transport'//trim(name)//'.pdb')

!Ecriture des atomes
      DO i_a=1,n_atom
            IF(tab_pos_a(i_a,4)==1) THEN
                  spinAt='Au'
                  WRITE(11,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                  spinAt,tab_pos_a(i_a,1),tab_pos_a(i_a,2),tab_pos_a(i_a,3),nul
            ENDIF

            IF(tab_pos_a(i_a,4)==-1) THEN
                  spinAt='Cu'
                  WRITE(11,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
                  spinAt,tab_pos_a(i_a,1),tab_pos_a(i_a,2),tab_pos_a(i_a,3),nul
            ENDIF
      ENDDO

!Ecriture des electrons
      DO i_e=1,n_e_itin
            spinAt='H'
            WRITE(11,'("ATOM",9X,A3,14X,3F8.3,1X,"1.00",1X,F8.3)') &
            spinAt,tab_pos_e(i_e,1),tab_pos_e(i_e,2),tab_pos_e(i_e,3),nul
      ENDDO

      CLOSE(11) 

      END SUBROUTINE write_conf3D_transport

!=====================================================================================================
!__ SUBROUTINE gauss : 
!              Determine l amplitude de la vitesse d un electron. ____________________________________
!_____________________________________________________________________________________________________
      SUBROUTINE gausss()

      IMPLICIT NONE

      CALL random_number(rdm_vit1)
            vitesse=rdm_vit1*aaa

      END SUBROUTINE gausss

!=====================================================================================================
!__ SUBROUTINE speed gauss : 
!              Determine la direction de delpacement de l electron. __________________________________
!_____________________________________________________________________________________________________
      SUBROUTINE speed_gauss()

      IMPLICIT NONE

!Tirage de cos(theta) biaise coordonnees spheriques
      CALL random_number(rdm_vit1)
      costheta=2.*rdm_vit1-1.
      sintheta=sqrt(1.-costheta**2.)

!Tirage de cos(phi) coordonnees spheriques
      CALL random_number(rdm_vit1)
      phi=2.*pi*rdm_vit1
      cosphi=cos(phi)
      sinphi=sin(phi)

!Calcul des vitesse sur chaque compotab_spin_ante d espace
      vx=vitesse*sintheta*cosphi
      vy=vitesse*sintheta*sinphi
      vz=vitesse*costheta

      END SUBROUTINE speed_gauss

!=====================================================================================================
!__ SUBROUTINE bound_cond_transport_old: 
!              Conditions periodiques et recherche des voisins de l electron non deplace _____________
!_____________________________________________________________________________________________________

      SUBROUTINE bound_cond_transport_old()

      IMPLICIT NONE

!Quy vi tri cua electron theo so maille
      alp=int(tab_pos_e(i_e,1)/aaa)+1     
      bet=alp+natx*int(tab_pos_e(i_e,2)/aaa2)
      gam_old=bet+natx*naty*int(tab_pos_e(i_e,3)/aaa3)

      pix=int(tab_pos_e(i_e,1)/aaa)
      piy=int(tab_pos_e(i_e,2)/aaa2)
      piz=int(tab_pos_e(i_e,3)/aaa3)
      prx=real(pix)*aaa
      pry=real(piy)*aaa
      prz=real(piz)*aaa

      compt=0
      tab_maille_vois_old(:)=0

      DO bcz=-del,del
            DO bcy=-del,del
                  DO bcx=-del,del
                  bcnm=0

!-------------------------------------------------------------------------------
!Condition periodique sur z
!-------------------------------------------------------------------------------
                  IF((prz+real(bcz)*aaa >= 0.).and. &
                     (prz+real(bcz)*aaa < (real(natz)-0.5)*aaa))THEN
                        bcnm=gam_old+natx*naty*bcz
!-------------------------------------------------------------------------------
!Condition periodique sur y
!----------------------------------------------------------------------
                        IF ((pry+real(bcy)*aaa >= 0.) .and. &
                            (pry+real(bcy)*aaa < real(naty)*aaa)) THEN
                            bcnm=bcnm+natx*bcy
                        ELSE 
                              IF (pry+real(bcy)*aaa >= real(naty)*aaa) THEN
                                    bcnm=bcnm-natx*(naty-abs(bcy))
                              ELSE
                                    bcnm=bcnm+natx*(naty-abs(bcy))
                              END IF
                        END IF
                                                           
!-------------------------------------------------------------------------------
!Condition periodique sur x
!-------------------------------------------------------------------------------
                        IF ((prx+real(bcx)*aaa >= 0.) .and. &
                           (prx+real(bcx)*aaa < real(natx)*aaa)) THEN
                              bcnm=bcnm+bcx
                        ELSE 
                              IF (prx+real(bcx)*aaa >= real(natx)*aaa) THEN
                                    bcnm=bcnm-(natx-abs(bcx))          
                              ELSE 
                                    bcnm=bcnm+(natx-abs(bcx))
                              END IF
                        END IF
                        compt=compt+1
                        tab_maille_vois_old(compt)=bcnm
!-------------------------------------------------------------------------------
!Fin des conditions periodiques et de recherche des voisins
                  END IF
                  ENDDO
            ENDDO
      ENDDO

      END SUBROUTINE bound_cond_transport_old


!=====================================================================================================
!__ SUBROUTINE bound cond transport : 
!              Deplacement virtuel & conditions periodiques sur 1 electrons. _________________________
!_____________________________________________________________________________________________________
      SUBROUTINE bound_cond_transport_new()

      IMPLICIT NONE

!----------------------------------------------------------------------------------------------
!Repositionnement de l electron dans la structure
!----------------------------------------------------------------------------------------------

      tab_pos_e_tmp(1)=tab_pos_e(i_e,1)+vx
      tab_pos_e_tmp(2)=tab_pos_e(i_e,2)+vy
      tab_pos_e_tmp(3)=tab_pos_e(i_e,3)+vz
      tab_pos_e_tmp(4)=tab_pos_e(i_e,4)

!-------------------------------------------------------------------------------
!Condition speculaire sur z
!-------------------------------------------------------------------------------
      IF((tab_pos_e_tmp(3) < (real(natz)-0.5)*aaa3) .and. tab_pos_e_tmp(3) >= 0)THEN
                tab_pos_e_tmp(3) =  tab_pos_e_tmp(3)                
      ELSE IF(tab_pos_e_tmp(3) >= (real(natz)-0.5)*aaa3)THEN
                tab_pos_e_tmp(3) = -tab_pos_e_tmp(3)+2.*(real(natz)-0.5)*aaa3
            ELSE IF(tab_pos_e_tmp(3) < 0)THEN
                  tab_pos_e_tmp(3) = -tab_pos_e_tmp(3)
      ENDIF
!-------------------------------------------------------------------------------
!Condition periodique sur y
!-------------------------------------------------------------------------------
      IF((tab_pos_e_tmp(2) < real(naty)*aaa2) .and. tab_pos_e_tmp(2) >= 0)THEN
                tab_pos_e_tmp(2) =  tab_pos_e_tmp(2)
      ELSE IF(tab_pos_e_tmp(2) >= real(naty)*aaa2)THEN
                tab_pos_e_tmp(2) =  tab_pos_e_tmp(2)-real(naty)*aaa2
           ELSE IF(tab_pos_e_tmp(2) < 0)THEN
                tab_pos_e_tmp(2) =  tab_pos_e_tmp(2)+real(naty)*aaa2
      ENDIF
!-------------------------------------------------------------------------------
!Condition periodique sur x
!-------------------------------------------------------------------------------
      IF((tab_pos_e_tmp(1) < real(natx)*aaa) .and. tab_pos_e_tmp(1) >= 0)THEN
                tab_pos_e_tmp(1) =  tab_pos_e_tmp(1)
      ELSE IF(tab_pos_e_tmp(1) >= real(natx)*aaa)THEN
                tab_pos_e_tmp(1) =  tab_pos_e_tmp(1)-real(natx)*aaa
           ELSE IF(tab_pos_e_tmp(1) < 0)THEN
                tab_pos_e_tmp(1) =  tab_pos_e_tmp(1)+real(natx)*aaa
      ENDIF

!Conversion de la position de l electron en numero de maille
      alp=int(tab_pos_e_tmp(1)/aaa)+1     
      bet=alp+natx*int(tab_pos_e_tmp(2)/aaa2)
      gam_new=bet+natx*naty*int(tab_pos_e_tmp(3)/aaa3)

!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------
!Conditions periodiques et recherche des voisins de l electron deplace
!----------------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------------

      pix=int(tab_pos_e_tmp(1)/aaa)
      piy=int(tab_pos_e_tmp(2)/aaa2)
      piz=int(tab_pos_e_tmp(3)/aaa3)
      prx=real(pix)*aaa
      pry=real(piy)*aaa
      prz=real(piz)*aaa

      compt=0

      IF (gam_new==gam_old) THEN
            tab_maille_vois_new(:)=tab_maille_vois_old(:)
      ELSE 
            tab_maille_vois_new(:)=0

            DO bcz=-del,del
                  DO bcy=-del,del
                        DO bcx=-del,del

                         bcnm=0

!-------------------------------------------------------------------------------
!Condition periodique sur z
!-------------------------------------------------------------------------------
                        IF ((prz+real(bcz)*aaa >= 0.) .and. &
                            (prz+real(bcz)*aaa < (real(natz)-0.5)*aaa)) THEN
                              bcnm=gam_new+natx*naty*bcz
!-------------------------------------------------------------------------------
!Condition periodique sur y
!-------------------------------------------------------------------------------
                              IF ((pry+real(bcy)*aaa >= 0.) .and. &
                                  (pry+real(bcy)*aaa < real(naty)*aaa)) THEN
                                  bcnm=bcnm+natx*bcy
                              ELSE 
                                    IF (pry+real(bcy)*aaa >= real(naty)*aaa) THEN
                                          bcnm=bcnm-natx*(naty-abs(bcy))
                                    ELSE 
                                          bcnm=bcnm+natx*(naty-abs(bcy))
                                    END IF
                              END IF
!-------------------------------------------------------------------------------
!Condition periodique sur x
!-------------------------------------------------------------------------------
                              IF ((prx+real(bcx)*aaa >= 0.) .and. &
                                    (prx+real(bcx)*aaa < real(natx)*aaa)) THEN
                                    bcnm=bcnm+bcx
                              ELSE 
                                    IF (prx+real(bcx)*aaa >= real(natx)*aaa) THEN
                                          bcnm=bcnm-(natx-abs(bcx))
                                    ELSE 
                                          bcnm=bcnm+(natx-abs(bcx))
                                    END IF
                              END IF

                              compt=compt+1
                              tab_maille_vois_new(compt)=bcnm
                        END IF
                        ENDDO
                  ENDDO
            ENDDO
      ENDIF

      END SUBROUTINE bound_cond_transport_new

!=====================================================================================================
!__ SUBROUTINE dmod elec elec : 
!              Calcul de la distance electron electron. ______________________________________________
!_____________________________________________________________________________________________________
      SUBROUTINE dmod_elec_elec_old()

      IMPLICIT NONE

      pvx=tab_pos_e(i_e_tmp,1)
      pvy=tab_pos_e(i_e_tmp,2)      

!Repositionnement electron non deplace (sur axe x)
      IF((tab_pos_e(i_e,1)>pvx).and.(abs(tab_pos_e(i_e,1)-pvx)>real(natx)*aaa/2.))THEN
            pvx=tab_pos_e(i_e_tmp,1)+real(natx)*aaa
      ELSE 
      IF((tab_pos_e(i_e,1)<pvx).and.(abs(tab_pos_e(i_e,1)-pvx)>real(natx)*aaa/2.))THEN
            pvx=tab_pos_e(i_e_tmp,1)-real(natx)*aaa
      END IF
      END IF      
!Repositionnement electron non deplace (sur axe y)
      IF((tab_pos_e(i_e,2)>pvy).and.(abs(tab_pos_e(i_e,2)-pvy)>real(naty)*aaa2/2.))THEN
            pvy=tab_pos_e(i_e_tmp,2)+real(naty)*aaa2
      ELSE 
      IF((tab_pos_e(i_e,2)<pvy).and.(abs(tab_pos_e(i_e,2)-pvy)>real(naty)*aaa2/2.))THEN
            pvy=tab_pos_e(i_e_tmp,2)-real(naty)*aaa2
      END IF
      END IF
!Calcul de la distance par rapport a l electron non deplace
      dmod_e_e_old=sqrt((pvx-tab_pos_e(i_e,1))**2.+(pvy-tab_pos_e(i_e,2))**2.+ &
                        (tab_pos_e(i_e,3)-tab_pos_e(i_e_tmp,3))**2.)

      END SUBROUTINE dmod_elec_elec_old  
!_______________________________________________________________________________________________
    
      SUBROUTINE dmod_elec_elec_new()

      IMPLICIT NONE

      pvx=tab_pos_e(i_e_tmp,1)
      pvy=tab_pos_e(i_e_tmp,2)

!Repositionnement electron deplace (sur axe x)
      IF((tab_pos_e_tmp(1)>pvx).and.(abs(tab_pos_e_tmp(1)-pvx)>real(natx)*aaa/2.))THEN
            pvx=tab_pos_e(i_e_tmp,1)+real(natx)*aaa
      ELSE  
      IF((tab_pos_e_tmp(1)<pvx) .and.(abs(tab_pos_e_tmp(1)-pvx)>real(natx)*aaa/2.))THEN
            pvx=tab_pos_e(i_e_tmp,1)-real(natx)*aaa
      END IF
      END IF
!Repositionnement electron deplace (sur axe y)
      IF((tab_pos_e_tmp(2)>pvy).and.(abs(tab_pos_e_tmp(2)-pvy)>real(naty)*aaa2/2.))THEN
            pvy=tab_pos_e(i_e_tmp,2)+real(naty)*aaa2
      ELSE 
      IF((tab_pos_e_tmp(2)<pvy).and.(abs(tab_pos_e_tmp(2)-pvy)>real(naty)*aaa2/2.))THEN
            pvy=tab_pos_e(i_e_tmp,2)-real(naty)*aaa2
      END IF
      END IF      

!Calcul de la distance par rapport a l electron deplace 
      dmod_e_e_new=sqrt((pvx-tab_pos_e_tmp(1))**2.+(pvy-tab_pos_e_tmp(2))**2.+ &
                        (tab_pos_e_tmp(3)-tab_pos_e(i_e_tmp,3))**2.)

      END SUBROUTINE dmod_elec_elec_new

!=====================================================================================================
!__ SUBROUTINE dmod elec atom : 
!              Calcul de la distance electron atome. _________________________________________________
!_____________________________________________________________________________________________________

      SUBROUTINE dmod_elec_atom_old()

      IMPLICIT NONE
      pvx=tab_pos_a(i_a,1)
      pvy=tab_pos_a(i_a,2)       

!Repositionnement electron non deplace (sur axe x)
      IF((tab_pos_e(i_e,1)>pvx).and.(abs(tab_pos_e(i_e,1)-pvx)>real(natx)*aaa/2.))THEN
            pvx=tab_pos_a(i_a,1)+real(natx)*aaa
      ELSE 
      IF((tab_pos_e(i_e,1)<pvx).and.(abs(tab_pos_e(i_e,1)-pvx)>real(natx)*aaa/2.))THEN
            pvx=tab_pos_a(i_a,1)-real(natx)*aaa
      END IF
      END IF
!Repositionnement electron non deplace (sur axe y)
      IF((tab_pos_e(i_e,2)>pvy).and.(abs(tab_pos_e(i_e,2)-pvy)>real(naty)*aaa2/2.))THEN
            pvy=tab_pos_a(i_a,2)+real(naty)*aaa2
      ELSE 
      IF((tab_pos_e(i_e,2)<pvy).and.(abs(tab_pos_e(i_e,2)-pvy)>real(naty)*aaa2/2.))THEN
            pvy=tab_pos_a(i_a,2)-real(naty)*aaa2
      END IF
      END IF

!Calcul de la distance par rapport a l electron non deplace
      dmod_e_a_old=sqrt((pvx-tab_pos_e(i_e,1))**2.+(pvy-tab_pos_e(i_e,2))**2.+ &
                      (tab_pos_e(i_e,3)-tab_pos_a(i_a,3))**2.)

      END SUBROUTINE dmod_elec_atom_old
!_______________________________________________________________________________________________

      SUBROUTINE dmod_elec_atom_new()

      IMPLICIT NONE

      pvx=tab_pos_a(i_a,1)
      pvy=tab_pos_a(i_a,2)

!Repositionnement electron deplace (sur axe x)
      IF((tab_pos_e_tmp(1)>pvx).and.(abs(tab_pos_e_tmp(1)-pvx)>real(natx)*aaa/2.))THEN
            pvx=tab_pos_a(i_a,1)+real(natx)*aaa
      ELSE 
      IF((tab_pos_e_tmp(1)<pvx).and.(abs(tab_pos_e_tmp(1)-pvx)>real(natx)*aaa/2.))THEN
            pvx=tab_pos_a(i_a,1)-real(natx)*aaa
      END IF
      END IF
!Repositionnement electron deplace (sur axe y)
      IF((tab_pos_e_tmp(2)>pvy).and.(abs(tab_pos_e_tmp(2)-pvy)>real(naty)*aaa2/2.))THEN
            pvy=tab_pos_a(i_a,2)+real(naty)*aaa2
      ELSE 
      IF((tab_pos_e_tmp(2)<pvy).and.(abs(tab_pos_e_tmp(2)-pvy)>real(naty)*aaa2/2.))THEN
            pvy=tab_pos_a(i_a,2)-real(naty)*aaa2
      END IF
      END IF

!Calcul de la distance par rapport a l electron deplace 
      dmod_e_a_new=sqrt((pvx-tab_pos_e_tmp(1))**2.+(pvy-tab_pos_e_tmp(2))**2.+ &
                      (tab_pos_e_tmp(3)-tab_pos_a(i_a,3))**2.)

      END SUBROUTINE dmod_elec_atom_new
!=====================================================================================================
!__ SUBROUTINE test pauli elec : 
!              recherche des electron present dans la sphere de Pauli. _______________________________
!_____________________________________________________________________________________________________    

      SUBROUTINE test_pauli_elec_old()

      IMPLICIT NONE
      
      IF(dmod_e_e_old <= 0.05*aaa)THEN
            test_pauli=1
            IF(bugvar/=0)THEN
                        WRITE(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
                        WRITE(*,*) 'BUG : pauli elec old -> test error'
                        WRITE(*,*) 'RUN WAS STOP'
                        WRITE(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
                        STOP
            ENDIF
      ENDIF

      END SUBROUTINE test_pauli_elec_old
!____________________________________________________________________________________________________

      SUBROUTINE test_pauli_elec_new()

      IMPLICIT NONE
      
      IF(dmod_e_e_new <= 0.05*aaa)THEN
              test_pauli=1
      ENDIF

      END SUBROUTINE test_pauli_elec_new
!=====================================================================================================
!__ SUBROUTINE test pauli atom : 
!              recherche des atomes present dans la sphere de Pauli. _________________________________
!_____________________________________________________________________________________________________
      SUBROUTINE test_pauli_atom_old()

      IMPLICIT NONE

      IF(dmod_e_a_old <= 0.05*aaa)THEN
            test_pauli=1
            IF(bugvar/=0)THEN
                        WRITE(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
                        WRITE(*,*) 'BUG : pauli atom old -> test error'
                        WRITE(*,*) 'RUN WAS STOP'
                        WRITE(*,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
                        STOP
            ENDIF
      ENDIF

      END SUBROUTINE test_pauli_atom_old
!____________________________________________________________________________________________________

      SUBROUTINE test_pauli_atom_new()

      IMPLICIT NONE

      IF(dmod_e_a_new <= 0.05*aaa)THEN
            test_pauli=1
      ENDIF

      END SUBROUTINE test_pauli_atom_new

!=====================================================================================================
!__ SUBROUTINE energy elec elec : 
!              Calcul de l energie d interaction electron-electron. __________________________________
!_____________________________________________________________________________________________________
      SUBROUTINE energy_elec_elec_new()

      IMPLICIT NONE

      IF(dmod_e_e_new <= d2*aaa)THEN
      energy_e_e_new=energy_e_e_new  &
                    - K0*exp(-dmod_e_e_new/aaa)*tab_pos_e_tmp(4)*tab_pos_e(i_e_tmp,4)
      density_new=density_new+1.
      ENDIF 

      END SUBROUTINE energy_elec_elec_new

!_______________________________________________________________________________________________

      SUBROUTINE energy_elec_elec_old()

      IMPLICIT NONE

      IF(dmod_e_e_old <= d2*aaa)THEN
      energy_e_e_old=energy_e_e_old-K0*exp(-dmod_e_e_old/aaa)*tab_pos_e(i_e,4)*tab_pos_e(i_e_tmp,4)
      density_old=density_old+1.

      ENDIF 

      END SUBROUTINE energy_elec_elec_old
!=====================================================================================================
!__ SUBROUTINE energy elec atom : 
!              Calcul de l energie d interaction electron-atome. _____________________________________
!_____________________________________________________________________________________________________
      SUBROUTINE energy_elec_atom_new()

      IMPLICIT NONE

      IF(dmod_e_a_new <= d1*aaa)THEN
      energy_e_a_new=energy_e_a_new  &
                    - I0*exp(-dmod_e_a_new/aaa)*tab_pos_e_tmp(4)*tab_pos_a(i_a,4)
      ENDIF

      END SUBROUTINE energy_elec_atom_new
!_______________________________________________________________________________________________

      SUBROUTINE energy_elec_atom_old()

      IMPLICIT NONE

      IF(dmod_e_a_old <= d1*aaa)THEN
      energy_e_a_old=energy_e_a_old &
                    - I0*exp(-dmod_e_a_old/aaa)*tab_pos_e(i_e,4)*tab_pos_a(i_a,4)
      ENDIF

      END SUBROUTINE energy_elec_atom_old
!=====================================================================================================
!__ SUBROUTINE sweep_elec : 
!              Balayage des voisins de l electron selectionne dans les mailles
!              voisines de tab_spin_a nouvelle et ancienne maille. ___________________________________________
!_____________________________________________________________________________________________________
      SUBROUTINE energy_elec_new()

      IMPLICIT NONE

!Initialitab_spin_ation variables
      test_pauli=0
      energy_e_e_new=0.
      energy_e_a_new=0.
      density_new=0.

!----------------------------------------------------------------------------------------------
!Calcul energie sur electron deplace
!----------------------------------------------------------------------------------------------
!Boucle sur les mailles voisines de celle de l electron selectionne
      i_m=1
      DO WHILE ((tab_maille_vois_new(i_m)/=0) .and. (test_pauli==0))

            bcnm=tab_maille_vois_new(i_m)

!Test sur les electrons voisins des mailles voisines de celle de l electron
!selectionne
            itme=1

            DO WHILE ((tab_maille_e(bcnm,itme)/=0) .and. (test_pauli==0))

                  i_e_tmp=tab_maille_e(bcnm,itme)

                  IF((tab_maille_e(bcnm,itme) > 0) .and. (i_e/=i_e_tmp))THEN
                        CALL dmod_elec_elec_new()
                        CALL test_pauli_elec_new()
                        CALL energy_elec_elec_new() 
                  ENDIF

                  itme=itme+1
            ENDDO

            i_m=i_m+1     
      ENDDO

!Boucle sur les mailles voisines de celle de l electron selectionne
      i_m=1
      DO WHILE ((tab_maille_vois_new(i_m)/=0) .and. (test_pauli==0))
            bcnm=tab_maille_vois_new(i_m)
            DO l=1,n_atom_p_m
                  i_a=tab_maille_a(bcnm,l)
                  CALL dmod_elec_atom_new()
                  CALL test_pauli_atom_new()
                  CALL energy_elec_atom_new()
            ENDDO
            i_m=i_m+1  
      ENDDO

      END SUBROUTINE energy_elec_new

!_______________________________________________________________________________________________
      SUBROUTINE energy_elec_old()

      IMPLICIT NONE
      
      test_pauli=0
      energy_e_e_old=0.
      energy_e_a_old=0.
      density_old=0.

!Tinh nang luong e_e cua electron dang xet voi tat ca cac electrons voisins
      i_m=1
      DO WHILE(tab_maille_vois_old(i_m)/=0)
            bcnm=tab_maille_vois_old(i_m)
            
            itme=1
            DO WHILE(tab_maille_e(bcnm,itme)/=0)
                  i_e_tmp=tab_maille_e(bcnm,itme)
                  IF((tab_maille_e(bcnm,itme) > 0) .and. (i_e/=i_e_tmp))THEN
                        CALL dmod_elec_elec_old()
                        CALL test_pauli_elec_old()
                        CALL energy_elec_elec_old() 
                  ENDIF

                  itme=itme+1
            ENDDO

            i_m=i_m+1     
      ENDDO

!Tinh nang luong e_a cua electron dang xet voi tat ca cac atoms voisins
      i_m=1
      DO WHILE(tab_maille_vois_old(i_m)/=0)
            bcnm=tab_maille_vois_old(i_m)
            
            DO l=1,n_atom_p_m
                  i_a=tab_maille_a(bcnm,l)
                  CALL dmod_elec_atom_old()
                  CALL test_pauli_atom_old()
                  CALL energy_elec_atom_old()
            ENDDO
            i_m=i_m+1  
      ENDDO

      END SUBROUTINE energy_elec_old
!=====================================================================================================
!__ SUBROUTINE arrange : 
!              Repositionne les electrons apres deplacement dans les tableaux. _______________________
!_____________________________________________________________________________________________________
      SUBROUTINE arrange()

      IMPLICIT NONE

!Conversion de la position de l electron en numero de maille
      alp=int(tab_pos_e(i_e,1)/aaa)+1     
      bet=alp+natx*int(tab_pos_e(i_e,2)/aaa2)
      gam=bet+natx*naty*int(tab_pos_e(i_e,3)/aaa3)

!Repositionnement de l electron dans tab_maille_e
      itme=1
      DO WHILE(tab_maille_e(gam,itme)/=i_e)
            itme=itme+1
      ENDDO

!Libere l ancienne position dans tab_maille_e
      IF((tab_maille_e(gam,itme+1)/=0) .and. (itme<n_e_itin))THEN
            tab_maille_e(gam,itme)=-1
      ELSE 
      IF((tab_maille_e(gam,itme+1))==0 .and. (itme<n_e_itin))THEN
            tab_maille_e(gam,itme)=0
      ELSE 
      IF (itme==n_e_itin)THEN
           tab_maille_e(gam,itme)=0
      END IF
      END IF
      END IF

!stockage nouvelle position de l electron dans la structure
      tab_pos_e(i_e,1)=tab_pos_e_tmp(1)
      tab_pos_e(i_e,2)=tab_pos_e_tmp(2)
      tab_pos_e(i_e,3)=tab_pos_e_tmp(3)

!Conversion de la position de l electron en numero de maille
      alp=int(tab_pos_e(i_e,1)/aaa)+1     
      bet=alp+natx*int(tab_pos_e(i_e,2)/aaa2)
      gam=bet+natx*naty*int(tab_pos_e(i_e,3)/aaa3)

!Repositionnement de l electron dans tab_maille_e
      itme=1
      DO WHILE(tab_maille_e(gam,itme)>0)
            itme=itme+1
      ENDDO

      tab_maille_e(gam,itme)=i_e 

      END SUBROUTINE arrange
!=====================================================================================================
!__ SUBROUTINE resistance: 
!              Calcul de la resistance du materiau. __________________________________________________
!_____________________________________________________________________________________________________
      SUBROUTINE resistance()

      IMPLICIT NONE

!Comptage parois 1
      nr1=0.

      IF((tab_pos_e(i_e,1) <= (real(natx)*aaa/2.)).and. &
           (tab_pos_e_tmp(1)   > (real(natx)*aaa/2.)).and.  &
           (abs(tab_pos_e(i_e,1)-tab_pos_e_tmp(1))<(real(natx)*aaa/2.)))THEN
                nr1=1.
      ELSE IF((tab_pos_e(i_e,1) > (real(natx)*aaa/2.)).and. &
                (tab_pos_e_tmp(1)    <= (real(natx)*aaa/2.)).and.&
                (abs(tab_pos_e(i_e,1)-tab_pos_e_tmp(1))<(real(natx)*aaa/2.)))THEN
                nr1=-1.
      ENDIF

      nr1_tot=nr1_tot+nr1

!Comptage parois 2
      nr2=0.

      IF((tab_pos_e(i_e,1) <= (real(natx)*aaa/4.)).and. &
           (tab_pos_e_tmp(1)    > (real(natx)*aaa/4.)).and.  &
           (abs(tab_pos_e(i_e,1)-tab_pos_e_tmp(1))<(real(natx)*aaa/2.)))THEN
                nr2=1.
      ELSE IF((tab_pos_e(i_e,1) > (real(natx)*aaa/4.)).and. &
                (tab_pos_e_tmp(1)    <= (real(natx)*aaa/4.)).and.&
                (abs(tab_pos_e(i_e,1)-tab_pos_e_tmp(1))<(real(natx)*aaa/2.)))THEN
                nr2=-1.
      ENDIF

      nr2_tot=nr2_tot+nr2

!Comptage parois 3
      nr3=0.

      IF((tab_pos_e(i_e,1) <= (real(natx)*3.*aaa/4.)).and. &
           (tab_pos_e_tmp(1)    > (real(natx)*3.*aaa/4.)).and.  &
           (abs(tab_pos_e(i_e,1)-tab_pos_e_tmp(1))<(real(natx)*aaa/2.)))THEN
                nr3=1.
      ELSE IF((tab_pos_e(i_e,1) > (real(natx)*3.*aaa/4.)).and. &
                (tab_pos_e_tmp(1)    <= (real(natx)*3.*aaa/4.)).and.&
                (abs(tab_pos_e(i_e,1)-tab_pos_e_tmp(1))<(real(natx)*aaa/2.)))THEN
                nr3=-1.
      ENDIF

      nr3_tot=nr3_tot+nr3

      END SUBROUTINE resistance
!=====================================================================================================
!__ SUBROUTINE vitesse_moy: 
!              Calcul de la vitesse moyenne des electrons. ___________________________________________
!_____________________________________________________________________________________________________
      SUBROUTINE vitesse_moy()

      IMPLICIT NONE

      cpt_vit_moy=cpt_vit_moy+1.
      vit_moy_x=vit_moy_x+vx/aaa
      vit_moy_y=vit_moy_y+vy/aaa
      vit_moy_z=vit_moy_z+vz/aaa
      vit_moy  =vit_moy+vitesse/aaa

      END SUBROUTINE vitesse_moy
!=====================================================================================================
!__ SUBROUTINE magnetization_layer_by_layer: 
!              Calcul du nombre d electron couche par couche selon l axe z. __________________________
!_____________________________________________________________________________________________________
      SUBROUTINE magnetization_layer_by_layer()

      IMPLICIT NONE

!------------------------------------------------------------------------------
!Variables locales
      INTEGER (KIND=4) :: i,j,k
!------------------------------------------------------------------------------

      n_e_layer1=0.
      n_e_layer2=0.

!Comptage premiere demie maille
      midmaille=0.
      DO k=1,natz
            midmaille=midmaille+aaa
            DO j=1,naty
                  DO i=1,natx
                        gam=1+(i-1)+(j-1)*natx+(k-1)*natx*naty
                        itme=1
                        DO WHILE(tab_maille_e(gam,itme)/=0)
                              IF(tab_maille_e(gam,itme)/=-1)THEN
                                    bet=tab_maille_e(gam,itme)
                                    IF((tab_pos_e(bet,3) <= midmaille-aaa/2.) .and. &
                                       (midmaille-aaa/2. < (real(natz)-0.5)*aaa))THEN
                                          n_e_layer1=n_e_layer1+tab_pos_e(bet,4)
                                    ELSE IF((tab_pos_e(bet,3) < midmaille) .and. &
                                          (midmaille < (real(natz)-0.5)*aaa))THEN
                                           n_e_layer2=n_e_layer2+tab_pos_e(bet,4)
                                    ENDIF
                              ENDIF
                              itme=itme+1
                        ENDDO
                  ENDDO
            ENDDO
                
            WRITE(Ligne14,*) n_e_layer1,n_e_layer2
            WRITE(14,'(a)',advance='no') trim(Ligne14)
            n_e_layer1=0.
            n_e_layer2=0.
                
      ENDDO

      WRITE(Ligne14,*) ''
      WRITE(14,'(a)') trim(Ligne14)

      END SUBROUTINE magnetization_layer_by_layer
!=====================================================================================================
!__ SUBROUTINE distance_ee_histogram: 
!              Construction de l histogramme des distances entre electrons. __________________________
!_____________________________________________________________________________________________________
      SUBROUTINE distance_ee_histogram()

      IMPLICIT NONE

      IF(comptmc>1.)THEN
            case_histo=int(dmod_e_e_old/(aaa/100.))+1
            tab_dist_ee_histogram(case_histo)=tab_dist_ee_histogram(case_histo)+1.
      ENDIF

      END SUBROUTINE distance_ee_histogram
!=====================================================================================================
!__ SUBROUTINE cosO_CosP_histogram: 
!              Construction de l histogramme des distribution en fonction de 
!              Cos(theta) et cos(Phi). _______________________________________________________________
!_____________________________________________________________________________________________________
      SUBROUTINE cosO_CosP_histogram()

      IMPLICIT NONE

      case_histo=int(costheta/(1./630.))
      tab_theta_histogram(case_histo)=tab_theta_histogram(case_histo)+1.
      case_histo=int(phi/(2.*pi/630.))
      tab_phi_histogram(case_histo)=tab_phi_histogram(case_histo)+1.

      END SUBROUTINE cosO_CosP_histogram
!=====================================================================================================
!__ SUBROUTINE 3D_cosO_CosP_histogram: 
!              Construction de l histogramme de la distribution 3D en fonction  
!              de Cos(theta) et cos(Phi). ____________________________________________________________
!_____________________________________________________________________________________________________
      SUBROUTINE cosO_CosP_3D_histogram()

      IMPLICIT NONE

      case_histo_O=int(costheta/(1./510.))
      case_histo_P=int(phi/(2.*pi/630.))
      tab_thetaphi_histogram(case_histo_O,case_histo_P)= &
      tab_thetaphi_histogram(case_histo_O,case_histo_P)+1.

      END SUBROUTINE cosO_CosP_3D_histogram
!=====================================================================================================
!__ SUBROUTINE Vxyz_histogram: 
!              Construction de l histogramme de la distribution 3D en fonction  
!              de Cos(theta) et cos(Phi). ____________________________________________________________
!_____________________________________________________________________________________________________
      SUBROUTINE Vxyz_histogram()

      IMPLICIT NONE

      case_histo=int(vx/(2.*aaa/500.))
      tab_vx_histogram(case_histo)=tab_vx_histogram(case_histo)+1.

      case_histo=int(vy/(2.*aaa/500.))
      tab_vy_histogram(case_histo)=tab_vy_histogram(case_histo)+1.
      case_histo=int(vz/(2.*aaa/500.))
      tab_vz_histogram(case_histo)=tab_vz_histogram(case_histo)+1.

      END SUBROUTINE Vxyz_histogram

!=====================================================================================================
!__ SUBROUTINE drift: 
!              Routine de derrive dR/dT. _____________________________________________________________
!_____________________________________________________________________________________________________
      SUBROUTINE drift()

      IMPLICIT NONE

      IF(i_T > 1)THEN
            Rt2=Rho
            R_drift=(Rt2-Rt1)/dT
            Rt1=Rt2
      ELSE IF(i_T==1)THEN
            Rt1=Rho      
            R_drift=0.          
      ENDIF

      END SUBROUTINE drift

!=====================================================================================================
!__ SUBROUTINE equi_elec: 
!              Tao trang thai equilibre cho tat ca electron, 1 lan _________________________________
!_____________________________________________________________________________________________________

      SUBROUTINE equi_elec()
      IMPLICIT NONE
      
      
      wl=(exp(1.))**2.
      DO WHILE(wl/=1.)
            wl=sqrt(wl)
      
            DO i_e=1,n_e_itin
                  denom_gauss=1.
                  CALL bound_cond_transport_old()
                  CALL energy_elec_old()
                  energy_old=energy_e_e_old+energy_e_a_old-fieldB*tab_pos_e(i_e,4)+d*density_old/Vsph
      
                  test_pauli=1
                  DO WHILE(test_pauli/=0)
                        CALL gausss()
                        CALL speed_gauss()
                        CALL bound_cond_transport_new()
                        CALL energy_elec_new()
                  ENDDO
                  energy_new=energy_e_e_new+energy_e_a_new-fieldE*vx/aaa &
                                    -fieldB*tab_pos_e(i_e,4)+d*density_new/Vsph
                              
                  CALL random_number(rdm_mtrp)
                  IF((energy_new <= energy_old).OR.((wl*exp((energy_old-energy_new)/T))>rdm_mtrp))THEN
                        CALL arrange()
                  ENDIF
      
            ENDDO
      ENDDO

      END SUBROUTINE equi_elec
!=====================================================================================================
!__ SUBROUTINE SUBROUTINE equi_elec_1()
! 
!       Tao trang thai equilibre ban dau cua tat ca electron, n_equi_elec_1 lan. 
!_____________________________________________________________________________________________________

      SUBROUTINE equi_elec_1()
      IMPLICIT NONE

      DO i_equi_elec_1=1,n_equi_elec_1
            CALL equi_elec()
      ENDDO
      END SUBROUTINE equi_elec_1

!---------------------------------------------------------------------------------------------------
! SUBROUTINE equi_elec_2(): Tao lai trang thai equilibre cua tat ca electron tab_spin_au moi lan trasport

      SUBROUTINE equi_elec_2()
      IMPLICIT NONE

      DO i_equi_elec_2=1,n_equi_elec_2
            CALL equi_elec()
      ENDDO
      END SUBROUTINE equi_elec_2


!=====================================================================================================
!__ SUBROUTINE valeu_transport(): 
!              Tinh gia tri cua cac dai luong trong qua trinh transport ______________________________
!_____________________________________________________________________________________________________
     
      SUBROUTINE valeu_transport()
      IMPLICIT NONE

      DO i_elec_transport=1,n_transport_elec
            DO i_e=1,n_e_itin
                  min_energy=100000.
      
                  CALL bound_cond_transport_old()
                  CALL energy_elec_old()
                  energy_old=energy_e_e_old+energy_e_a_old &
                                    -fieldB*tab_pos_e(i_e,4)+d*density_old/Vsph

                  CALL gausss()
                  DO i_gauss=1,n_gauss
                        test_pauli=1

                        DO WHILE(test_pauli/=0)
                              CALL speed_gauss()
                              CALL bound_cond_transport_new()
                              CALL energy_elec_new()
                        ENDDO
                        energy_new=energy_e_e_new+energy_e_a_new-fieldE*vx/aaa &
                                         -fieldB*tab_pos_e(i_e,4)+d*density_new/Vsph

                        IF(energy_new < min_energy)THEN
                              vx_tmp=vx
                              vy_tmp=vy
                              vz_tmp=vz
                              min_energy=energy_new
                              costhetatmp=costheta
                              phitmp=phi
                        ENDIF
                  ENDDO

                  ! Chargement de la meilleur vitesse
                  vx=vx_tmp
                  vy=vy_tmp
                  vz=vz_tmp
                  costheta=costhetatmp
                  phi=phitmp

                  ! Tinh so lan eletron so 1 thoat ra khoi he
                  IF ((tab_pos_e(1,3)+vz)>real(natz)*aaa) THEN
                        n_traverse=n_traverse+1
                  END IF  
                  !------------------------------------------
            
                  CALL bound_cond_transport_new()
                  CALL energy_elec_new()
                  energy_new=energy_e_e_new+energy_e_a_new-fieldE*vx/aaa &
                               -fieldB*tab_pos_e(i_e,4)+d*density_new/Vsph             
            
                  CALL random_number(rdm_mtrp)
      
                  IF((energy_new <= energy_old).or. &
                              (exp((energy_old-energy_new)/T)>rdm_mtrp))THEN

                        energy_new_tot=energy_new_tot+energy_new-fieldE*vx/aaa &
                                             -fieldB*tab_pos_e(i_e,4)+d*density_new/Vsph
            
                        energy_old_tot=energy_old_tot+energy_old &
                                               -fieldB*tab_pos_e(i_e,4)+d*density_old/Vsph

                        CALL vitesse_moy()
                        CALL resistance()
                        CALL arrange()

                  ELSE
                        energy_new_tot=energy_new_tot+energy_old &
                                                -fieldB*tab_pos_e(i_e,4)+d*density_old/Vsph
                        energy_old_tot=energy_old_tot+energy_old &
                                                -fieldB*tab_pos_e(i_e,4)+d*density_old/Vsph
                  ENDIF   

            ENDDO
      ENDDO
      END SUBROUTINE valeu_transport


!=====================================================================================================
!__ SUBROUTINE average_transport(): 
!              Tinh gia tri TB cua cac dai luong trong qua trinh transport ___________________________
!_____________________________________________________________________________________________________


      SUBROUTINE average_transport()
      IMPLICIT NONE
      
      n_gauss=10
      vit_moy_x=0.
      vit_moy_y=0.
      vit_moy_z=0.
      Rho1=0.
      Rho2=0.
      Rho3=0.
      Rho=0.
      nr1_tot=0.
      nr2_tot=0.
      nr3_tot=0.
      vit_moy=0.
      cpt_vit_moy=0.
      cptpos=0.
      energy_new_tot=0.
      energy_old_tot=0.
      tab_dist_ee_histogram(:)=0.
      tab_theta_histogram(:)=0.
      tab_phi_histogram(:)=0.
      tab_thetaphi_histogram(:,:)=0.
      tab_vx_histogram(:)=0.
      tab_vy_histogram(:)=0.
      tab_vz_histogram(:)=0.
      n_traverse=0
      !----------------  
      
      CALL equi_elec_1()
      
      DO i_average_transport=1,n_average_transport
            IF (i_average_transport>1) THEN
                  CALL equi_reseau_3()
                  CALL reload_lattice()
                  CALL equi_elec_2()
            END IF

            CALL valeu_transport()
      END DO
      comptmc=real(n_e_itin*n_transport_elec*n_average_transport)
      Rho1=comptmc/nr1_tot
      Rho2=comptmc/nr2_tot
      Rho3=comptmc/nr3_tot

      Rho=(Rho1+Rho2+Rho3)/3.
      vit_moy_x=vit_moy_x/comptmc
      vit_moy_y=vit_moy_y/comptmc
      vit_moy_z=vit_moy_z/comptmc
      vit_moy=vit_moy/comptmc
      energy_new_tot=energy_new_tot/(2.*comptmc)
      energy_old_tot=energy_old_tot/(2.*comptmc)

      CALL drift()        

      WRITE(Ligne13,*) T,Rho,vit_moy_x,vit_moy_y,vit_moy_z,vit_moy, &
                        energy_new_tot,energy_old_tot,R_drift
      WRITE(13,'(a)') trim(Ligne13)
      !WRITE(15,*)T, n_traverse

     ! WRITE(*,*)time_relax,n_equi_reseau_2,n_equi_elec_2,n_transport_elec,n_average_transport

      END SUBROUTINE average_transport
!_____________________________________________________________________________________________________
!_____________________________________________________________________________________________________

     END PROGRAM main
