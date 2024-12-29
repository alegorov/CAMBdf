    module CAMB
    use GaugeInterface
    implicit none
    contains

    !Call this routine with a set of parameters to generate the results you want.
    subroutine CAMB_GetResults(OutData, Params, error)
    use CAMBmain
    type(CAMBdata)  :: OutData
    type(CAMBparams) :: Params
    integer, optional :: error !Zero if OK
    type(CAMBparams) P

    global_error_flag = 0
    call OutData%Free()
    call SetActiveState(OutData)

    if (Params%WantCls .and. Params%WantScalars) then
        P = Params
        P%Max_eta_k=max(min(P%max_l,3000)*2.5_dl,P%Max_eta_k)
        call OutData%SetParams(P)
        if (global_error_flag==0) call cmbmain
        if (global_error_flag/=0) then
            if (present(error)) error =global_error_flag
            return
        end if
    end if

    OutData%CP%WantCls = Params%WantCls
    OutData%CP%WantScalars = Params%WantScalars
    OutData%CP%Accuracy = Params%Accuracy
    OutData%CP%Reion%Reionization = Params%Reion%Reionization

    if (global_error_flag/=0 .and. present(error)) error =global_error_flag

    end subroutine CAMB_GetResults

    end module CAMB


    program CAMBdf
    use CAMB
    use DarkFieldModule
    use InitialPower
    use Recombination
    use Reionization
    implicit none

    ! Variables
    integer :: num_args, ix
    character(256), allocatable :: args(:)
    integer :: error
    type(CAMBdata)  :: OutData
    type(CAMBparams) :: P
    class(TDarkFieldModel), allocatable :: DarkField
    class(TInitialPowerLaw), allocatable :: InitPower
    class(TRecfast), allocatable :: Recomb
    class(TTanhReionization), allocatable :: Reion
    integer, parameter :: lmax = 2500
    real(dl), parameter :: norm_coef = 2725500.0_dl ** 2
    integer :: l
    character(50) :: col1, col2, col3, col4

    num_args = command_argument_count()
    allocate(args(num_args))
    do ix = 1, num_args
        call get_command_argument(ix,args(ix))
    end do

    allocate(TDarkFieldModel::DarkField)
    allocate(TInitialPowerLaw::InitPower)
    allocate(TRecfast::Recomb)
    allocate(TTanhReionization::Reion)

    read(args(1),*) P%ombh2                    !0.0223828_dl
    read(args(2),*) P%omch2                    !0.1201075_dl
    read(args(3),*) P%H0                       !67.32117_dl
    read(args(4),*) Reion%optical_depth        !0.05430842_dl
    read(args(5),*) InitPower%ns               !0.9660499_dl
    read(args(6),*) InitPower%As               !2.100549e-09_dl
    read(args(7),*) DarkField%alpha
    read(args(8),*) DarkField%delta

    P%WantCls = .true.
    P%WantScalars = .true.
    P%min_l = 2
    P%max_l = 2500
    P%max_eta_k = 6250.0_dl
    P%max_eta_k_tensor = 1200.0_dl
    P%omnuh2 = 0.0006451439_dl
    P%TCMB = 2.7255_dl
    P%YHe = 0.2454006_dl
    P%num_nu_massless = 2.046_dl
    P%num_nu_massive = 1
    P%nu_mass_eigenstates = 1
    P%nu_mass_degeneracies(1) = 0.0_dl
    P%nu_mass_fractions(1) = 1.0_dl
    P%nu_mass_numbers(1) = 1
    InitPower%nrun = 0.0_dl
    InitPower%nrunrun = 0.0_dl
    InitPower%pivot_scalar = 0.05_dl
    Recomb%RECFAST_fudge = 1.125_dl
    Recomb%RECFAST_fudge_He = 0.86_dl
    Recomb%RECFAST_Heswitch = 6
    Recomb%RECFAST_Hswitch = .true.
    Recomb%AGauss1 = -0.14_dl
    Recomb%AGauss2 = 0.079_dl
    Recomb%zGauss1 = 7.28_dl
    Recomb%zGauss2 = 6.73_dl
    Recomb%wGauss1 = 0.18_dl
    Recomb%wGauss2 = 0.33_dl
    Reion%Reionization = .true.
    Reion%use_optical_depth = .true.
    Reion%redshift = 10.0_dl
    Reion%fraction = -1.0_dl
    Reion%include_helium_fullreion = .true.
    Reion%helium_redshift = 3.5_dl
    Reion%helium_delta_redshift = 0.5_dl
    Reion%helium_redshiftstart = 6.0_dl
    Reion%tau_solve_accuracy_boost = 1.0_dl
    Reion%timestep_boost = 1.0_dl
    Reion%max_redshift = 50.0_dl
    Reion%delta_redshift = 0.5_dl
    P%Accuracy%AccuracyBoost = 1.0_dl
    P%Accuracy%lSampleBoost = 1.0_dl
    P%Accuracy%lAccuracyBoost = 1.0_dl
    P%Accuracy%TimeStepBoost = 1.0_dl
    P%Accuracy%BackgroundTimeStepBoost = 1.0_dl
    P%Accuracy%IntTolBoost = 1.0_dl
    P%Accuracy%SourcekAccuracyBoost = 1.0_dl
    P%Accuracy%IntkAccuracyBoost = 1.0_dl
    P%Accuracy%BessIntBoost = 1.0_dl
    P%Accuracy%BesselBoost = 1.0_dl
    P%Accuracy%neutrino_q_boost = 1.0_dl
    P%Alens = 1.0_dl
    P%MassiveNuMethod = Nu_trunc
    P%min_l_logl_sampling = 5000

    P%DarkField = DarkField
    P%InitPower = InitPower
    P%Recomb = Recomb
    P%Reion = Reion

    error = 0

    call CAMB_GetResults(OutData, P, error)

    do l = State%CP%Min_l, min(State%CP%Max_l, lmax)
        write(col1,*) l
        write(col2,*) norm_coef * State%CLData%Cl_Scalar(l, C_Temp)
        write(col3,*) norm_coef * State%CLData%Cl_Scalar(l, C_E)
        write(col4,*) norm_coef * State%CLData%Cl_Scalar(l, C_Cross)
        print '(A)', trim(col1)//' '//trim(col2)//' '//trim(col3)//' '//trim(col4)
    end do

    end program CAMBdf

