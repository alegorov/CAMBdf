    !     This this is the main CAMB program module.
    !
    !     Code for Anisotropies in the Microwave Background
    !     by Antony lewis (http://cosmologist.info) and Anthony Challinor
    !     See readme.html for documentation.

    !     Note that though the code is internally parallelised, it is not thread-safe
    !     so you cannot generate more than one model at the same time in different threads.
    !
    !     Based on CMBFAST  by  Uros Seljak and Matias Zaldarriaga, itself based
    !     on Boltzmann code written by Edmund Bertschinger, Chung-Pei Ma and Paul Bode.
    !     Original CMBFAST copyright and disclaimer:
    !
    !     Copyright 1996 by Harvard-Smithsonian Center for Astrophysics and
    !     the Massachusetts Institute of Technology.  All rights reserved.
    !
    !     THIS SOFTWARE IS PROVIDED "AS IS", AND M.I.T. OR C.f.A. MAKE NO
    !     REPRESENTATIONS OR WARRANTIES, EXPRESS OR IMPlIED.
    !     By way of example, but not limitation,
    !     M.I.T. AND C.f.A MAKE NO REPRESENTATIONS OR WARRANTIES OF
    !     MERCHANTABIlITY OR FITNESS FOR ANY PARTICUlAR PURPOSE OR THAT
    !     THE USE OF THE lICENSED SOFTWARE OR DOCUMENTATION WIll NOT INFRINGE
    !     ANY THIRD PARTY PATENTS, COPYRIGHTS, TRADEMARKS OR OTHER RIGHTS.
    !
    !     portions of this software are based on the COSMICS package of
    !     E. Bertschinger.  See the lICENSE file of the COSMICS distribution
    !     for restrictions on the modification and distribution of this software.

    module CAMBmain

    !     This code evolves the linearized perturbation equations of general relativity,
    !     the Boltzmann equations and the fluid equations for perturbations
    !     of a Friedmann-Robertson-Walker universe with a supplied system of gauge-dependent equation
    !     in a modules called GaugeInterface.  The sources for the line of sight integral are
    !     computed at sampled times during the evolution for various of wavenumbers. The sources
    !     are then interpolated to a denser wavenumber sampling for computing the line of
    !     sight integrals of the form Integral d(conformal time) source_k * bessel_k_l.
    !     For _flat models the bessel functions are interpolated from a pre-computed table, for
    !     non-_flat models the hyperspherical Bessel functions are computed by integrating their
    !     differential equation. Both phases ('Evolution' and 'Integration') can do separate
    !     wavenumbers in parallel.

    !     The time variable is conformal  time dtau=dt/a(t) and the spatial dependence is Fourier transformed
    !     with q=sqrt(k**2 + (|m|+1)K), comoving distances are x=r/a(t), with a(t)=1 today.
    !     The units of both length and time are Mpc.

    !    Many elements are part of derived types (to make thread safe or to allow non-sequential code use
    !    CP = CAMB parameters
    !    EV = Time evolution variables
    !    IV = Source integration variables
    !    CT = Cl transfer data
    !    MT = matter transfer data


    use precision
    use results
    use GaugeInterface
    use SpherBessels
    use MassiveNu
    use InitialPower
    use Recombination
    use RangeUtils
    use constants
    use MathUtils
    implicit none
    private

    logical ExactClosedSum  !do all nu values in sum for Cls for Omega_k>0.1

    !Variables for integrating the sources with the bessel functions for each wavenumber
    type IntegrationVars
        integer q_ix
        real(dl) q, dq    !q value we are doing and delta q
        !Contribution to C_l integral from this k
        real(dl), dimension(:,:), pointer :: Source_q, ddSource_q
        !Interpolated sources for this k
        integer SourceSteps !number of steps up to where source is zero
    end type IntegrationVars

    real(dl), dimension(:,:), allocatable :: iCl_scalar, iCl_vector, iCl_tensor
    ! Cls at the l values we actually compute,  iCl_xxx(l_index, Cl_type, initial_power_index)

    real(dl), dimension(:,:,:), allocatable :: iCl_Array
    !Full covariance at each L (alternative more general arrangement to above)

    real(dl),parameter :: qmin0=0.1_dl

    real(dl) :: dtaurec_q

    !     qmax - CP%Max_eta_k/State%tau0, qmin = qmin0/State%tau0 for _flat case
    real(dl) qmin, qmax

    real(dl) max_etak_tensor , max_etak_vector, max_etak_scalar
    !     Will only be calculated if k*tau < max_etak_xx

    integer maximum_l !Max value of l to compute
    real(dl) :: maximum_qeta = 3000._dl

    integer :: l_smooth_sample = 3000 !assume transfer functions effectively small for k*chi0>2*l_smooth_sample

    integer :: max_bessels_l_index  = 1000000
    real(dl) :: max_bessels_etak = 1000000*2

    Type(TTimeSources) , pointer :: ThisSources => null()
    Type(TTimeSources), allocatable, target :: TempSources

    real(dl), dimension(:,:,:), allocatable :: ddScaledSrc !temporary source second derivative array
    real(dl), dimension(:,:,:), pointer ::  ScaledSrc => null() !linear source with optional non-linear scaling

    integer n_source_points ! number of CL source wavenumbers (for use when calculted remaining non-CL transfers)

    procedure(obj_function), private :: dtauda

    public cmbmain, TimeSourcesToCl, ClTransferToCl, InitVars, GetTauStart !InitVars for BAO hack

    contains


    subroutine cmbmain
    integer q_ix
    type(EvolutionVars) EV
    Type(TTimer) :: Timer
    real(dl) starttime
    Type(ClTransferData), pointer :: ThisCT 

    if (CP%WantCls) then
        !Use CAMB_GetResults instead

        maximum_l = CP%Max_l
        maximum_qeta = CP%Max_eta_k
    end if

    if (DebugMsgs .and. Feedbacklevel > 0) call Timer%Start(starttime)

    call InitVars(State) !Most of single thread time spent here (in InitRECFAST)
    if (global_error_flag/=0) return

    if (DebugMsgs .and. Feedbacklevel > 0) then
        call Timer%WriteTime('Timing for InitVars')
    end if

    call CP%InitPower%Init(CP)
    if (global_error_flag/=0) return

    !Calculation of the CMB and other sources.
    if (CP%WantCls) then
        if (CP%WantScalars) then
            !Allow keeping time source for scalars, so e.g. different initial power spectra
            ! and non-linear corrections can be applied later
            if (.not. allocated(State%ScalarTimeSources)) allocate(State%ScalarTimeSources)
            ThisSources => State%ScalarTimeSources
        else
            if (.not. allocated(TempSources)) allocate(TempSources)
            ThisSources => TempSources
        end if
        call SetkValuesForSources
    end if

    !***note that !$ is the prefix for conditional multi-processor compilation***
    !$ if (ThreadNum /=0) call OMP_SET_NUM_THREADS(ThreadNum)

    if (CP%WantCls) then
        if (DebugMsgs .and. Feedbacklevel > 0) call WriteFormat('Set %d source k values', &
            ThisSources%Evolve_q%npoints)

        ThisCT => State%ClData%CTransScal

        call GetSourceMem

        ThisCT%NumSources = ThisSources%SourceNum
        call ThisCT%ls%Init(State,CP%Min_l, maximum_l)

        !$OMP PARALLEL DO DEFAULT(SHARED),SCHEDULE(DYNAMIC), PRIVATE(EV, q_ix)
        do q_ix= ThisSources%Evolve_q%npoints,1,-1
            if (global_error_flag==0) call DoSourcek(EV,q_ix)
        end do
        !$OMP END PARALLEL DO

        if (DebugMsgs .and. Feedbacklevel > 0) call Timer%WriteTime('Timing for source calculation')

    endif !WantCls

    !     if CMB calculations are requested, calculate the Cl by
    !     integrating the sources over time and over k.

    if (CP%WantCls) then
        call TimeSourcesToCl(ThisCT)

        if (CP%WantScalars) then
            deallocate(State%ScalarTimeSources)
        else
            deallocate(TempSources)
        end if
        nullify(ThisSources)
    end if
    if (DebugMsgs .and. Feedbacklevel > 0) then
        call Timer%WriteTime('Timing for whole of cmbmain', starttime)
    end if

    end subroutine cmbmain

    subroutine TimeSourcesToCl(ThisCT)
    Type(ClTransferData) :: ThisCT 
    integer q_ix
    Type(TTimer) :: Timer

    if (CP%WantScalars) ThisSources => State%ScalarTimeSources

    if (DebugMsgs .and. Feedbacklevel > 0) call Timer%Start()

    ScaledSrc => ThisSources%LinearSrc

    if (global_error_flag==0) then
        call InitSourceInterpolation

        ExactClosedSum = State%scale < 0.93_dl

        max_bessels_l_index = ThisCT%ls%nl
        max_bessels_etak  = maximum_qeta

        if (CP%WantScalars) call GetLimberTransfers(ThisCT)
        ThisCT%max_index_nonlimber = max_bessels_l_index

        call InitSpherBessels(ThisCT%ls, CP, max_bessels_l_index,max_bessels_etak )
        !This is only slow if not called before with same (or higher) Max_l, Max_eta_k
        !Preferably stick to Max_l being a multiple of 50

        call SetkValuesForInt(ThisCT)

        if (DebugMsgs .and. Feedbacklevel > 0) call WriteFormat('Set %d integration k values',ThisCT%q%npoints)

        !Begin k-loop and integrate Sources*Bessels over time
        !$OMP PARALLEL DO DEFAULT(SHARED), SCHEDULE(STATIC,4)
        do q_ix=1,ThisCT%q%npoints
            call SourceToTransfers(ThisCT, q_ix)
        end do !q loop
        !$OMP END PARALLEL DO

        if (DebugMsgs .and. Feedbacklevel > 0) call Timer%WriteTime('Timing for Integration')
    end if

    if (allocated(ddScaledSrc)) deallocate(ddScaledSrc)
    if (associated(ScaledSrc) .and. .not. associated(ScaledSrc,ThisSources%LinearSrc)) then
        deallocate(ScaledSrc)
        nullify(ScaledSrc)
    end if

    !Final calculations for CMB output unless want the Cl transfer functions only.
    if (global_error_flag==0) &
        call ClTransferToCl(State)

    if (DebugMsgs .and. Feedbacklevel > 0) call Timer%WriteTime('Timing for final CL output')

    end subroutine TimeSourcesToCl

    subroutine ClTransferToCl(State)
    class(CAMBdata) :: State

    call SetActiveState(State)
    if (State%CP%WantScalars .and. State%CP%WantCls .and. global_error_flag==0) then
        allocate(iCl_Scalar(State%CLdata%CTransScal%ls%nl,C_Temp:State%Scalar_C_last), source=0._dl)

        call CalcScalCls(State%CLdata%CTransScal)
        if (DebugMsgs .and. Feedbacklevel > 0) write (*,*) 'CalcScalCls'
    end if

    if (global_error_flag==0) then
        call State%CLdata%InitCls(State)
        !     Calculating Cls for every l.
        call InterpolateCls()
        if (DebugMsgs .and. Feedbacklevel > 0) write (*,*) 'InterplolateCls'
    end if

    if (CP%WantScalars .and. allocated(iCl_Scalar)) deallocate(iCl_scalar)
    if (CP%WantScalars .and. allocated(iCl_Array)) deallocate(iCl_Array)

    end subroutine ClTransferToCl

    subroutine GetLimberTransfers(ThisCT)
    Type(ClTransferData), target :: ThisCT 
    integer ell, ell_needed
    integer i, s_ix, s_ix_lens
    integer n1,n2,n, ell_limb
    real(dl) int,k, chi
    integer klo, khi
    real(dl) a0,b0,ho2o6,a03,b03,ho,reall
    Type(LimberRec), pointer :: LimbRec

    call Init_Limber(ThisCT,State)

    end subroutine GetLimberTransfers

    subroutine SourceToTransfers(ThisCT, q_ix)
    type(ClTransferData), target :: ThisCT 
    integer q_ix
    type(IntegrationVars) :: IV

    allocate(IV%Source_q(State%TimeSteps%npoints,ThisSources%SourceNum))

    call IntegrationVars_init(IV)

    IV%q_ix = q_ix
    IV%q =ThisCT%q%points(q_ix)
    IV%dq= ThisCT%q%dpoints(q_ix)

    call InterpolateSources(IV)

    call DoSourceIntegration(IV, ThisCT)

    deallocate(IV%Source_q)

    end subroutine SourceToTransfers

    function GetTauStart(q)
    real(dl), intent(IN) :: q
    real(dl) taustart, GetTauStart

    !     Begin when wave is far outside horizon.
    !     Conformal time (in Mpc) in the radiation era, for photons plus 3 species
    !     of relativistic neutrinos.
    taustart=0.001_dl/q

    !     Make sure to start early in the radiation era.
    taustart=min(taustart,0.1_dl)

    !     Start when massive neutrinos are strongly relativistic.
    if (CP%Num_nu_massive>0 .and. any(State%nu_masses(1:CP%Nu_mass_eigenstates)/=0)) then
        taustart=min(taustart,1.d-3/maxval(State%nu_masses(1:CP%Nu_mass_eigenstates))/State%adotrad)
    end if

    GetTauStart=taustart
    end function GetTauStart

    subroutine DoSourcek(EV,q_ix)
    integer q_ix
    real(dl) taustart
    type(EvolutionVars) EV

    EV%q=ThisSources%Evolve_q%points(q_ix)

    EV%q2=EV%q**2

    EV%q_ix = q_ix
    EV%TransferOnly=.false.

    EV%ThermoData => State%ThermoData

    taustart = GetTauStart(EV%q)

    call GetNumEqns(EV)

    if (CP%WantScalars .and. global_error_flag==0) call CalcScalarSources(EV,taustart)

    end subroutine DoSourcek

    subroutine GetSourceMem
    integer :: err

    if (CP%WantScalars) then
        ThisSources%SourceNum=2
        ThisSources%NonCustomSourceNum = ThisSources%SourceNum
        State%Scalar_C_last = C_Cross
    else
        ThisSources%SourceNum=3
        ThisSources%NonCustomSourceNum = ThisSources%SourceNum
    end if

    if (allocated(ThisSources%LinearSrc)) &
        deallocate(ThisSources%LinearSrc)
    allocate(ThisSources%LinearSrc(ThisSources%Evolve_q%npoints,&
        ThisSources%SourceNum,State%TimeSteps%npoints), source=0._dl, stat=err)
    if (err/=0) call GlobalError('Sources requires too much memory to allocate', &
        error_unsupported_params)                                                                               

    end subroutine GetSourceMem



    !  initial variables, number of steps, etc.
    subroutine InitVars(state)
    type(CAMBdata) :: state
    real(dl) taumin, maxq, initAccuracyBoost
    integer itf

    call SetActiveState(state)

    initAccuracyBoost = CP%Accuracy%AccuracyBoost * CP%Accuracy%TimeStepBoost

    ! Maximum and minimum k-values.
    qmax=maximum_qeta/State%tau0
    qmin=qmin0/State%tau0/initAccuracyBoost
    !     Timesteps during recombination (tentative, the actual
    !     timestep is the minimum between this value and taurst/40,
    !     where taurst is the time when recombination starts - see inithermo

    dtaurec_q=4/qmax/initAccuracyBoost
    !AL:Changed Dec 2003, dtaurec feeds back into the non-_flat integration via the step size
    State%dtaurec = dtaurec_q
    !dtau rec may be changed by ThermoData_init

    max_etak_tensor = initAccuracyBoost*maximum_qeta /10
    max_etak_scalar = initAccuracyBoost*max(1700._dl,maximum_qeta) /20
    if (maximum_qeta <3500 .and. CP%Accuracy%AccuracyBoost < 2) max_etak_scalar = max_etak_scalar * 1.5
    !tweak to get large scales right
    max_etak_vector = max_etak_scalar

    if (CP%WantCls) then
        maxq = qmax
    end if

    taumin=GetTauStart(maxq)

    !     Initialize baryon temperature and ionization fractions vs. time.
    !     This subroutine also fixes the timesteps where the sources are
    !     saved in order to do the integration. So TimeSteps is set here.
    !These routines in ThermoData (modules.f90)
    call State%ThermoData%Init(State,taumin)
    if (global_error_flag/=0) return

    if (DebugMsgs .and. Feedbacklevel > 0) write (*,*) 'ThermoData.Init'

    !Do any array initialization for propagation equations
    call GaugeInterface_Init

    if (Feedbacklevel > 0)  &
        write(*,'("tau_recomb/Mpc       = ",f7.2,"  tau_now/Mpc = ",f8.1)') State%tau_maxvis,State%tau0

    end subroutine InitVars

    subroutine SetkValuesForSources
    implicit none
    real(dl) dlnk0, dkn1, dkn2, q_switch, q_cmb, dksmooth
    real(dl) qmax_log
    real(dl) SourceAccuracyBoost
    !     set k values for which the sources for the anisotropy and
    !     polarization will be calculated. For low values of k we
    !     use a logarithmic spacing. _closed case dealt with by SetClosedkValues

    SourceAccuracyBoost = CP%Accuracy%AccuracyBoost * CP%Accuracy%SourcekAccuracyBoost
    if (CP%WantScalars .and. CP%Reion%Reionization) then
        dlnk0=2._dl/10/SourceAccuracyBoost
        !Need this to get accurate low l polarization
    else
        dlnk0=5._dl/10/SourceAccuracyBoost
    end if

    dlnk0 = dlnk0/2

    dkn1=0.6_dl/State%taurst/SourceAccuracyBoost
    dkn2=0.9_dl/State%taurst/SourceAccuracyBoost/1.2

    qmax_log = dkn1/dlnk0
    q_switch = 2*6.3/State%taurst
    !Want linear spacing for wavenumbers which come inside horizon
    !Could use sound horizon, but for tensors that is not relevant

    q_cmb = 2*l_smooth_sample/State%chi0*SourceAccuracyBoost  !assume everything is smooth at l > l_smooth_sample
    if (maximum_l > 5000) q_cmb = q_cmb*1.4
    q_cmb = max(q_switch*2, q_cmb)
    !prevent EE going wild in tail
    dksmooth = q_cmb/2/(SourceAccuracyBoost)**2
    dksmooth = dksmooth/6

    associate(Evolve_q => ThisSources%Evolve_q)
        call Evolve_q%Init()
        call Evolve_q%Add_delta(qmin, qmax_log, dlnk0, IsLog = .true.)
        if (qmax > qmax_log) &
            call Evolve_q%Add_delta(qmax_log, min(qmax,q_switch), dkn1)
        if (qmax > q_switch) then
            call Evolve_q%Add_delta(q_switch, min(q_cmb,qmax), dkn2)
            if (qmax > q_cmb) then
                dksmooth = log(1 + dksmooth/q_cmb)
                call Evolve_q%Add_delta(q_cmb, qmax, dksmooth, IsLog = .true.)
            end if
        end if

        call Evolve_q%GetArray(.false.)
    end associate

    end subroutine SetkValuesForSources


    subroutine CalcScalarSources(EV,taustart)
    use FileUtils
    use MpiUtils
    implicit none
    type(EvolutionVars) EV
    real(dl) tau,tol1,tauend, taustart
    integer j,ind,itf
    real(dl) c(24),w(EV%nvar,9), y(EV%nvar), sources(ThisSources%SourceNum)

    w=0
    y=0
    call initial(EV,y, taustart)
    if (global_error_flag/=0) return

    tau=taustart
    ind=1

    !     Begin timestep loop.
    itf=1
    tol1=tol/exp(CP%Accuracy%AccuracyBoost*CP%Accuracy%IntTolBoost-1)

    do j=2,State%TimeSteps%npoints
        tauend=State%TimeSteps%points(j)

        if (.not. DebugEvolution .and. (EV%q*tauend > max_etak_scalar .and. tauend > State%taurend)) then
            ThisSources%LinearSrc(EV%q_ix,:,j)=0
        else
            !Integrate over time, calulate end point derivs and calc output
            call GaugeInterface_EvolveScal(EV,tau,y,tauend,tol1,ind,c,w)
            if (global_error_flag/=0) return

            call output(EV,y,j, tau,sources)
            ThisSources%LinearSrc(EV%q_ix,:,j)=sources
        end if
    end do !time step loop

    end subroutine CalcScalarSources


    subroutine InitSourceInterpolation
    integer i,j
    !     get the interpolation matrix for the sources to interpolate them
    !     for other k-values

    if (allocated(ddScaledSrc)) deallocate(ddScaledSrc)
    allocate(ddScaledSrc(ThisSources%Evolve_q%npoints,ThisSources%SourceNum,State%TimeSteps%npoints))
    !$OMP PARALLEL DO DEFAULT(SHARED), SCHEDULE(STATIC), PRIVATE(i,j)
    do  i=1,State%TimeSteps%npoints
        do j=1, ThisSources%SourceNum
            call spline_def(ThisSources%Evolve_q%points,ScaledSrc(:,j,i), &
                ThisSources%Evolve_q%npoints, ddScaledSrc(:,j,i))
        end do
    end do
    !$OMP END PARALLEL DO

    end subroutine InitSourceInterpolation


    subroutine SetkValuesForInt(ThisCT)
    Type(ClTransferData) :: ThisCT 
    integer no
    real(dl) dk,dk0,dlnk1, dk2, max_k_dk, k_max_log, k_max_0
    integer lognum
    real(dl)  qmax_int,IntSampleBoost


    qmax_int = min(qmax,max_bessels_etak/State%tau0)

    IntSampleBoost=CP%Accuracy%AccuracyBoost*CP%Accuracy%IntkAccuracyBoost
    if (do_bispectrum) then
        IntSampleBoost = IntSampleBoost * 2
        if (hard_bispectrum) IntSampleBoost = IntSampleBoost * 2
    end if

    !     Fixing the # of k for the integration.

    call ThisCT%q%Init()

    !Split up into logarithmically spaced intervals from qmin up to k=lognum*dk0
    !then no-lognum*dk0 linearly spaced at dk0 up to no*dk0
    !then at dk up to max_k_dk, then dk2 up to qmax_int
    lognum=nint(10*IntSampleBoost)
    dlnk1=1._dl/lognum
    no=nint(600*IntSampleBoost)
    dk0=1.8_dl/State%chi0/IntSampleBoost
    dk=3._dl/State%chi0/IntSampleBoost/1.6

    k_max_log = lognum*dk0
    k_max_0  = no*dk0

    if (do_bispectrum) k_max_0 = max(10.d0,k_max_0)

    dk2 = 0.04/IntSampleBoost  !very small scales

    call ThisCT%q%Add_delta(qmin, k_max_log, dlnk1, IsLog = .true.)
    call ThisCT%q%Add_delta(k_max_log, min(qmax_int,k_max_0), dk0)

    if (qmax_int > k_max_0) then
        max_k_dk = max(3000, 2*maximum_l)/State%tau0

        call ThisCT%q%Add_delta(k_max_0, min(qmax_int, max_k_dk), dk)
        if (qmax_int > max_k_dk) then
            !This allows inclusion of high k modes for computing BB lensed spectrum accurately
            !without taking ages to compute.
            call ThisCT%q%Add_delta(max_k_dk, qmax_int, dk2)
        end if
    end if

    call Init_ClTransfer(ThisCT)

    end subroutine setkValuesForInt

    subroutine InterpolateSources(IV)
    implicit none
    integer i,khi,klo, step
    real(dl) xf,b0,ho,a0,ho2o6,a03,b03
    type(IntegrationVars) IV
    Type(TRanges), pointer :: Evolve_q

    Evolve_q => ThisSources%Evolve_q

    !     finding position of k in table Evolve_q to do the interpolation.

    !Can't use the following in _closed case because regions are not set up (only points)
    !           klo = min(Evolve_q%npoints-1,Evolve_q%IndexOf(IV%q))
    !This is a bit inefficient, but thread safe
    klo=1
    do while ((IV%q > Evolve_q%points(klo+1)).and.(klo < (Evolve_q%npoints-1)))
        klo=klo+1
    end do

    khi=klo+1


    ho=Evolve_q%points(khi)-Evolve_q%points(klo)
    a0=(Evolve_q%points(khi)-IV%q)/ho
    b0=(IV%q-Evolve_q%points(klo))/ho
    ho2o6 = ho**2/6
    a03=(a0**3-a0)
    b03=(b0**3-b0)
    IV%SourceSteps = 0

    !     Interpolating the source as a function of time for the present
    !     wavelength.
    step=2
    do i=2, State%TimeSteps%npoints
        xf=IV%q*(State%tau0-State%TimeSteps%points(i))

        if (CP%WantScalars) then
            if ((DebugEvolution .or. IV%q*State%TimeSteps%points(i) < max_etak_scalar) &
                .and. xf > 1.e-8_dl) then
                step=i
                IV%Source_q(i,:) = a0 * ScaledSrc(klo,:,i) +  b0 * ScaledSrc(khi,:,i) + (a03*ddScaledSrc(klo,:,i) + &
                    b03 * ddScaledSrc(khi,:,i)) * ho2o6
            else
                IV%Source_q(i,:) = 0
            end if
        end if
    end do
    IV%SourceSteps = step

    end subroutine


    subroutine IntegrationVars_Init(IV)
    type(IntegrationVars), intent(INOUT) :: IV

    IV%Source_q(1,:)=0
    IV%Source_q(State%TimeSteps%npoints,:) = 0
    IV%Source_q(State%TimeSteps%npoints-1,:) = 0

    end  subroutine IntegrationVars_Init


    subroutine DoSourceIntegration(IV, ThisCT) !for particular wave number q
    type(IntegrationVars) IV
    Type(ClTransferData) :: ThisCT    
    integer j,ll,llmax
    real(dl) nu
    real(dl) :: sixpibynu

    nu=IV%q
    sixpibynu  = 6._dl*const_pi/nu

    llmax=nint(nu*State%chi0)
    if (llmax<15) then
        llmax=17 !AL Sept2010 changed from 15 to get l=16 smooth
    else
        llmax=nint(nu*State%rofChi(State%tau0 + sixpibynu))
    end if

    call DoFlatIntegration(IV,ThisCT, llmax)

    end subroutine DoSourceIntegration

    !_flat source integration
    subroutine DoFlatIntegration(IV, ThisCT, llmax)
    implicit none
    type(IntegrationVars) IV
    Type(ClTransferData) :: ThisCT 
    integer llmax
    integer j
    logical DoInt
    real(dl) xlim,xlmax1
    real(dl) tmin, tmax
    real(dl) a2, J_l, aa(IV%SourceSteps), fac(IV%SourceSteps)
    real(dl) xf, sums(ThisSources%SourceNum)
    real(dl) qmax_int
    integer bes_ix,n, bes_index(IV%SourceSteps)
    integer custom_source_off, s_ix
    integer nwin
    real(dl) :: BessIntBoost

    BessIntBoost = CP%Accuracy%AccuracyBoost*CP%Accuracy%BessIntBoost
    custom_source_off = 4

    !     Find the position in the xx table for the x correponding to each
    !     timestep

    do j=1,IV%SourceSteps !Precompute arrays for this k
        xf=abs(IV%q*(State%tau0-State%TimeSteps%points(j)))
        bes_index(j)=BessRanges%IndexOf(xf)
        !Precomputed values for the interpolation
        bes_ix= bes_index(j)
        fac(j)=BessRanges%points(bes_ix+1)-BessRanges%points(bes_ix)
        aa(j)=(BessRanges%points(bes_ix+1)-xf)/fac(j)
        fac(j)=fac(j)**2*aa(j)/6
    end do

    do j=1,max_bessels_l_index
        if (ThisCT%ls%l(j) > llmax) return
        xlim=xlimfrac*ThisCT%ls%l(j)
        xlim=max(xlim,xlimmin)
        xlim=ThisCT%ls%l(j)-xlim
        if (full_bessel_integration .or. do_bispectrum) then
            tmin = State%TimeSteps%points(2)
        else
            xlmax1=80*ThisCT%ls%l(j)*BessIntBoost
            tmin=State%tau0-xlmax1/IV%q
            tmin=max(State%TimeSteps%points(2),tmin)
        end if
        tmax=State%tau0-xlim/IV%q
        tmax=min(State%tau0,tmax)
        tmin=max(State%TimeSteps%points(2),tmin)

        if (tmax < State%TimeSteps%points(2)) exit
        sums = 0

        !As long as we sample the source well enough, it is sufficient to
        !interpolate the Bessel functions only

        if (ThisSources%SourceNum==2) then
            !This is the innermost loop, so we separate the no lensing scalar case to optimize it
            do n= State%TimeSteps%IndexOf(tmin),min(IV%SourceSteps,State%TimeSteps%IndexOf(tmax))
                a2=aa(n)
                bes_ix=bes_index(n)

                J_l=a2*ajl(bes_ix,j)+(1-a2)*(ajl(bes_ix+1,j) - ((a2+1) &
                    *ajlpr(bes_ix,j)+(2-a2)*ajlpr(bes_ix+1,j))* fac(n)) !cubic spline

                J_l = J_l*State%TimeSteps%dpoints(n)
                sums(1) = sums(1) + IV%Source_q(n,1)*J_l
                sums(2) = sums(2) + IV%Source_q(n,2)*J_l
            end do
        else
            qmax_int= max(850,ThisCT%ls%l(j))*3*BessIntBoost/State%tau0*1.2
            DoInt = .not. CP%WantScalars .or. IV%q < qmax_int
            !Do integral if any useful contribution to the CMB, or large scale effects

            if (DoInt) then
                do n= State%TimeSteps%IndexOf(tmin),min(IV%SourceSteps,State%TimeSteps%IndexOf(tmax))
                    !Full Bessel integration
                    a2=aa(n)
                    bes_ix=bes_index(n)

                    J_l=a2*ajl(bes_ix,j)+(1-a2)*(ajl(bes_ix+1,j) - ((a2+1) &
                        *ajlpr(bes_ix,j)+(2-a2)*ajlpr(bes_ix+1,j))* fac(n)) !cubic spline
                    J_l = J_l*State%TimeSteps%dpoints(n)

                    !The unwrapped form is faster
                    sums(1) = sums(1) + IV%Source_q(n,1)*J_l
                    sums(2) = sums(2) + IV%Source_q(n,2)*J_l
                    sums(3) = sums(3) + IV%Source_q(n,3)*J_l
                end do
            end if
            if (.not. DoInt .and. ThisSources%NonCustomSourceNum>3) then
                if (any(ThisCT%limber_l_min(4:ThisSources%NonCustomSourceNum)==0 .or. &
                    ThisCT%limber_l_min(4:ThisSources%NonCustomSourceNum) > j)) then
                    !When CMB does not need integral but other sources do
                    do n= State%TimeSteps%IndexOf(State%ThermoData%tau_start_redshiftwindows), &
                        min(IV%SourceSteps, State%TimeSteps%IndexOf(tmax))
                        !Full Bessel integration
                        a2 = aa(n)
                        bes_ix = bes_index(n)

                        J_l = a2 * ajl(bes_ix, j) + (1 - a2) * (ajl(bes_ix + 1, j) -&
                            ((a2 + 1) * ajlpr(bes_ix, j) + (2 - a2) * &
                            ajlpr(bes_ix + 1, j)) * fac(n)) !cubic spline
                        J_l = J_l * State%TimeSteps%dpoints(n)

                        sums(4) = sums(4) + IV%Source_q(n, 4) * J_l
                        do s_ix = 5, ThisSources%NonCustomSourceNum
                            sums(s_ix) = sums(s_ix) + IV%Source_q(n, s_ix) * J_l
                        end do
                    end do
                end if
            end if
        end if

        ThisCT%Delta_p_l_k(:,j,IV%q_ix) = ThisCT%Delta_p_l_k(:,j,IV%q_ix) + sums
    end do

    end subroutine DoFlatIntegration

    subroutine CalcScalCls(CTrans)
    implicit none
    Type(ClTransferData) :: CTrans
    integer j, q_ix, w_ix2
    real(dl) apowers
    real(dl) dlnk, ell, ctnorm, dbletmp, Delta1, Delta2
    real(dl), allocatable :: ks(:), dlnks(:), pows(:)
    real(dl) fac(3)
    integer nscal, i

    allocate(ks(CTrans%q%npoints),dlnks(CTrans%q%npoints), pows(CTrans%q%npoints))
    do q_ix = 1, CTrans%q%npoints
        ks(q_ix) = CTrans%q%points(q_ix)
        dlnks(q_ix) = CTrans%q%dpoints(q_ix)/CTrans%q%points(q_ix)
        pows(q_ix) = CP%InitPower%ScalarPower(ks(q_ix))
        if (global_error_flag/=0) return
    end do

    !Can't use OpenMP here on ifort. Does not terminate.
    !$OMP PARALLEL DO DEFAULT(SHARED), SCHEDULE(STATIC,4), &
    !$OMP PRIVATE(ell,q_ix,dlnk,apowers,ctnorm,dbletmp,Delta1,Delta2,w_ix,w_ix2,fac, nscal, i)
    do j=1,CTrans%ls%nl
        !Integrate dk/k Delta_l_q**2 * Power(k)
        ell = real(CTrans%ls%l(j),dl)
        if (j<= CTrans%max_index_nonlimber) then
            do q_ix = 1, CTrans%q%npoints
                !cut off at nu = l + 1
                dlnk = dlnks(q_ix)
                apowers = pows(q_ix)

                iCl_scalar(j,C_Temp:C_E) = iCl_scalar(j,C_Temp:C_E) +  &
                    apowers*CTrans%Delta_p_l_k(1:2,j,q_ix)**2*dlnk
                iCl_scalar(j,C_Cross) = iCl_scalar(j,C_Cross) + &
                    apowers*CTrans%Delta_p_l_k(1,j,q_ix)*CTrans%Delta_p_l_k(2,j,q_ix)*dlnk


                if (CTrans%NumSources>2 ) then
                    if (CTrans%limber_l_min(3)== 0 .or. j<CTrans%limber_l_min(3)) then
                        iCl_scalar(j,C_Phi) = iCl_scalar(j,C_Phi) +  &
                            apowers*CTrans%Delta_p_l_k(3,j,q_ix)**2*dlnk
                        iCl_scalar(j,C_PhiTemp) = iCl_scalar(j,C_PhiTemp) +  &
                            apowers*CTrans%Delta_p_l_k(3,j,q_ix)*CTrans%Delta_p_l_k(1,j,q_ix)*dlnk
                        iCl_scalar(j,C_PhiE) = iCl_scalar(j,C_PhiE) +  &
                            apowers*CTrans%Delta_p_l_k(3,j,q_ix)*CTrans%Delta_p_l_k(2,j,q_ix)*dlnk
                    end if
                end if
            end do

        end if !limber (j<= max_bessels_l_index)

        !Output l(l+1)C_l/OutputDenominator
        ctnorm=(ell*ell-1)*(ell+2)*ell
        dbletmp=(ell*(ell+1))/OutputDenominator*const_fourpi

        iCl_scalar(j,C_Temp)  =  iCl_scalar(j,C_Temp)*dbletmp
        iCl_scalar(j,C_E) =  iCl_scalar(j,C_E)*dbletmp*ctnorm
        iCl_scalar(j,C_Cross) =  iCl_scalar(j,C_Cross)*dbletmp*sqrt(ctnorm)
        if (CTrans%NumSources>2) then
            iCl_scalar(j,C_Phi) = CP%ALens*iCl_scalar(j,C_Phi)*const_fourpi*ell**4
            !The lensing power spectrum computed is l^4 C_l^{\phi\phi}
            !We put pix extra factors of l here to improve interpolation in l
            iCl_scalar(j,C_PhiTemp) = sqrt(CP%ALens)*  iCl_scalar(j,C_PhiTemp)*const_fourpi*ell**3
            !Cross-correlation is CTrans%ls%l^3 C_l^{\phi T}
            iCl_scalar(j,C_PhiE) = sqrt(CP%ALens)*  iCl_scalar(j,C_PhiE)*const_fourpi*ell**3*sqrt(ctnorm)
            !Cross-correlation is CTrans%ls%l^3 C_l^{\phi E}
        end if
    end do
    !$OMP END PARALLEL DO

    end subroutine CalcScalCls


    subroutine InterpolateCls()
    implicit none
    integer i,j
    integer, dimension(2,2), parameter :: ind = reshape( (/ 1,3,3,2 /), shape(ind))
    !use verbose form above for gfortran consistency  [[1,3],[3,2]]

    !Note using log interpolation is worse

    if (CP%WantScalars) then
        associate(lSet=>State%CLdata%CTransScal%ls)
            do i = C_Temp, State%Scalar_C_last
                call lSet%InterpolateClArrTemplated(iCl_scalar(1,i),State%CLData%Cl_scalar(lSet%lmin, i), &
                    lSet%nl,i)
            end do
        end associate
    end if

    end subroutine InterpolateCls


    end module CAMBmain
