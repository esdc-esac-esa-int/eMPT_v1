! updated to match eMPT retooled to handle MPT roll assignment
! updated to handle dependency of FORE/OTE transformation on order separation filter
!
!---------------------------------------------------------------------------'
!--- eMPT Modules Include File Development Version 2.0 24 June 2022  -------'
!---------------------------------------------------------------------------'

!-------------------------------------------------------
    module path_to_ref_files

      implicit none
      character*70,save ::  ref_path='./reference_files/'  

    end module path_to_ref_files

!-------------------------------------------------------
    module config_parameters

      implicit none
      integer,save :: ntile
      character*70,save ::  confname
      character*20,save ::  disperser
      integer,save ::  n_disp,n_filt
      real,save ::  sthresh
      real,save ::  raccx,raccy
      character*70,save ::  catfile
      character*100,save ::  imname
      character*100,save ::  segname
      real*8,save ::  cra,cdec
      real,save ::  cpa_ap,cpa_v3
      real,save ::  szone_x,szone_y
      real,save ::  angle_to_target,dva_mag
      integer,save ::  max_c_score1_1,max_c_score23_1,max_c_score45_1
      integer,save ::  max_c_score1_2,max_c_score23_2,max_c_score45_2
      integer,save ::  n_dither
      integer,save ::  mtarget
      logical,save ::  skypad
      
    end module config_parameters
!-------------------------------------------------------
    module derived_parameters
 
      real,save  :: av_xpitch !mas
      real,save  :: av_ypitch !mas
      real,save  :: av_xopen  !mas
      real,save  :: av_yopen  !mas
      real,save  :: phi !deg
      real,save  :: phi_r !rad
      real,save  :: xm_ref !mm
      real,save  :: ym_ref !mm
      real,save  :: ax_refv2,ay_refv3
      real,save  :: avxgap
    
    end module derived_parameters
!-------------------------------------------------------
    module pointings_and_groups
      
      implicit none
      integer,save :: n_p
      integer,allocatable,save :: n_targs_p(:),ids_p(:,:),no_p(:)
      real*8,allocatable,save :: ra_p(:),dec_p(:)
      real,allocatable,save :: pa_ap_p(:),pa_v3_p(:)

      integer,save :: max_targs,min_targs,np_max_targs,np_min_targs,ng_max_targs,ng_min_targs
      integer,save :: n_groups,max_groups
      integer,allocatable,save :: gr_flag(:)
      integer,allocatable,save :: gr_pcount(:),gr_tcount(:),gr_targs(:),gr_indx(:),gr_inv(:),gr_ptcount(:)
      real,allocatable,save :: gr_pcountr(:)
      integer,allocatable,save :: gr_subs(:,:)
      logical,allocatable,save :: p_free(:)
      integer,save :: min_pointings_in_group,min_numb_groups,min_ok_gr_targs,max_ok_gr_targs,ng_max_ok,ng_min_ok,ng_all_ok

    end module pointings_and_groups
!-------------------------------------------------------
    module target_cat

      implicit none
      integer,save :: ntarg
      integer,save :: ntarg_1
      integer,save :: nclass
      integer,save :: n_reach
      integer,allocatable,save :: id_t(:)
      integer,allocatable,save :: pri_t(:)
      real*8,allocatable,save :: ra_t(:),dec_t(:) 
      logical,allocatable,save :: targ_poss(:)
      real*8, save :: ra_c,dec_c !catalog median reference position
                
    end module target_cat
!-------------------------------------------------------
    module k_list
      
      implicit none     
      character*70 k_list_mod_file
      character*70 k_list_file
      integer n_max_k
      integer, allocatable :: nk(:)
      integer, allocatable :: id(:,:),id_cat(:,:)
      integer, allocatable :: pri_k(:,:),kt(:,:),it(:,:),jt(:,:)
      real, allocatable :: rx(:,:),ry(:,:)

    end module k_list
!-------------------------------------------------------
    module m_list
      
      implicit none     
      character*70,save :: m_list_file
      integer,save :: n_max_m

    end module m_list
!-------------------------------------------------------
    module order_and_weights
      
      implicit none     
      integer, allocatable :: order_matrix(:,:)
      integer, allocatable :: place_order(:,:)
      real min_dx,max_dx,min_dy,max_dy
      real, allocatable :: weights(:)

    end module order_and_weights
!-------------------------------------------------------
    module fits_image

      implicit none
      integer,save :: nx,ny
      real,allocatable,save :: image(:,:)
      integer,allocatable,save :: segmap(:,:)
      real,save :: datamin,datamax
      real*8,save :: crpix1,crpix2,crval1,crval2,cd1_1,cd1_2,cd2_1,cd2_2
                
    end module fits_image
!-------------------------------------------------------
    module histograms
      
      implicit none
      integer, allocatable :: hist_k1(:),hist_k2(:),hist_k3(:)
      integer, allocatable :: hist_m1(:),hist_m2(:),hist_m3(:)
      integer, allocatable :: hist_cat(:),hist_com(:)
      real, allocatable :: avexp(:)
      integer n_com,n_spec

      end module histograms
!-------------------------------------------------------
