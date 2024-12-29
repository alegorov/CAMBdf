    module model
    use Precision
    use classes
    use constants, only : COBE_CMBTemp, default_nnu
    use DarkFieldModule
    use MassiveNu
    use config
    use iso_c_binding
    implicit none

    integer, parameter :: outNone=1

    integer, parameter :: neutrino_hierarchy_normal = 1, neutrino_hierarchy_inverted = 2, neutrino_hierarchy_degenerate = 3

    integer, parameter :: Nu_int = 0, Nu_trunc=1, Nu_approx = 2, Nu_best = 3
    !For CAMBparams%MassiveNuMethod
    !Nu_int: always integrate distribution function
    !Nu_trunc: switch to expansion in velocity once non-relativistic
    !Nu_approx: approximate scheme - good for CMB, but not formally correct and no good for matter power
    !Nu_best: automatically use mixture which is fastest and most accurate

    integer, parameter :: max_Nu = 5 !Maximum number of neutrino species

    type AccuracyParams
        !Parameters for checking/changing overall accuracy
        !parameters equal to 1 corresponds to ~0.1% scalar C_l accuracy (at L>600)

        real(dl) :: AccuracyBoost =1._dl
        !Decrease step sizes, etc. by this parameter. Useful for checking accuracy.
        !Can also be used to improve speed significantly if less accuracy is required.
        !or improving accuracy for extreme models.
        !Note this does not increase lSamp%l sampling or massive neutrino q-sampling


        real(dl) :: lSampleBoost=1._dl
        !Increase lSampleBoost to increase sampling in lSamp%l for Cl interpolation

        real(dl) :: lAccuracyBoost=1._dl
        !Boost number of multipoles integrated in Boltzman heirarchy
      
        !The following allow separate tweaking (all also affected by AccuracyBoost above)

        real(dl) :: TimeStepBoost = 1._dl !sampling timesteps

        real(dl) :: BackgroundTimeStepBoost = 1._dl !number of time steps for background thermal history interpolation

        real(dl) :: IntTolBoost = 1._dl !Tolerances for integrating differential equations

        real(dl) :: SourcekAccuracyBoost = 1._dl !Accuracy of k sampling for source time integration

        real(dl) :: IntkAccuracyBoost = 1._dl !Accuracy of k sampling for integration

        real(dl) :: BessIntBoost = 1._dl !Accuracy of bessel integration truncation

        real(dl) :: BesselBoost = 1._dl !accuracy of bessel pre-computation sampling

        real(dl) :: neutrino_q_boost = 1._dl !number of momenta integrated for neutrino perturbations

    end type AccuracyParams

    !Non-linear corrections, either just P(k), or just CMB lensing/sources, or both
    integer, parameter :: NonLinear_none=0, NonLinear_Pk =1, NonLinear_Lens=2
    integer, parameter :: NonLinear_both=3

    ! Main parameters type
    type, extends (TCAMBParameters) :: CAMBparams
        logical   :: WantCls  = .true.

        logical   :: WantScalars = .true.

        integer   :: Min_l = 2 ! 1 or larger, usually 1 or 2
        integer   :: Max_l = 2500
        real(dl)  :: Max_eta_k = 5000
        real(dl)  :: Max_eta_k_tensor = 1200
        ! _tensor settings only used in initialization,
        !Max_l and Max_eta_k are set to the tensor variables if only tensors requested

        real(dl)  :: ombh2 = 0._dl !baryon density Omega_b h^2
        real(dl)  :: omch2 = 0._dl !cold dark matter density Omega_c h^2
        real(dl)  :: omnuh2 = 0._dl !massive neutino Omega_nu h^2
        real(dl)  :: H0 = 67._dl !Hubble parameter in km/s/Mpc
        real(dl)  :: TCMB = COBE_CMBTemp
        real(dl)  :: Yhe = 0.24_dl
        real(dl)  :: Num_Nu_massless = default_nnu
        integer   :: Num_Nu_massive = 0 !sum of Nu_mass_numbers below
        integer   :: Nu_mass_eigenstates = 0  !1 for degenerate masses
        real(dl)  :: Nu_mass_degeneracies(max_nu)
        real(dl)  :: Nu_mass_fractions(max_nu) !The ratios of the total densities
        integer   :: Nu_mass_numbers(max_nu) !physical number per eigenstate

        class(TInitialPower), allocatable :: InitPower
        class(TRecombinationModel), allocatable :: Recomb
        class(TReionizationModel), allocatable :: Reion
        class(TDarkFieldModel), allocatable :: DarkField
        type(AccuracyParams)     :: Accuracy

        real(dl), allocatable :: z_outputs(:) !Redshifts to output background outputs

        real(dl)  :: Alens = 1._dl !Unphysical rescaling parameter of the CMB lensing power

        integer   :: MassiveNuMethod = Nu_best

        integer :: min_l_logl_sampling = 5000 ! increase to use linear sampling for longer
    end type CAMBparams

    end module model
