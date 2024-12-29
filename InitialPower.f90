    !This module provides the initial power spectra

    !TInitialPowerLaw is parameterized as an expansion in ln k
    !
    ! ln P_s = ln A_s + (n_s -1)*ln(k/k_0_scalar) + n_{run}/2 * ln(k/k_0_scalar)^2 + n_{runrun}/6 * ln(k/k_0_scalar)^3
    !
    ! so if n_{run} = 0, n_{runrun}=0
    !
    ! P_s = A_s (k/k_0_scalar)^(n_s-1)
    !
    !for the scalar spectrum, when n_s=an(in) is the in'th spectral index. k_0_scalar
    !is a pivot scale, fixed here to 0.05/Mpc (change it below as desired or via .ini file).
    !
    !The tensor spectrum has three different supported parameterizations giving
    !
    ! ln P_t = ln A_t + n_t*ln(k/k_0_tensor) + n_{t,run}/2 * ln(k/k_0_tensor)^2
    !
    ! _tensor_parameterization==_tensor_param_indeptilt (=1) (default, same as CAMB pre-April 2014)
    !
    ! A_t = r A_s
    !
    ! _tensor_parameterization==_tensor_param_rpivot (=2)
    !
    ! A_t = r P_s(k_0_tensor)
    !
    ! _tensor_parameterization==_tensor_param_AT (=3)
    !
    ! A_t =  tensor_amp
    !
    !The absolute normalization of the Cls is unimportant here, but the relative ratio
    !of the tensor and scalar Cls generated with this module will be correct for general models

    module InitialPower
    use Precision
    use MpiUtils, only : MpiStop
    use classes
    implicit none

    private

    Type, extends(TInitialPower) :: TInitialPowerLaw
        !For the default implementation return power spectra based on spectral indices
        real(dl) :: ns = 1._dl !scalar spectral indices
        real(dl) :: nrun = 0._dl !running of spectral index
        real(dl) :: nrunrun  = 0._dl !running of spectral index
        real(dl) :: pivot_scalar = 0.05_dl !pivot scales in Mpc^{-1}
        real(dl) :: As = 1._dl
    contains
    procedure :: Init => TInitialPowerLaw_Init
    procedure :: ScalarPower => TInitialPowerLaw_ScalarPower
    end Type TInitialPowerLaw

    public TInitialPowerLaw
    contains

    subroutine TInitialPowerLaw_Init(this, Params)
    class(TInitialPowerLaw) :: this
    class(TCAMBParameters), intent(in) :: Params


    end subroutine TInitialPowerLaw_Init

    function TInitialPowerLaw_ScalarPower(this, k)
    class(TInitialPowerLaw) :: this
    real(dl), intent(in) :: k
    real(dl) TInitialPowerLaw_ScalarPower
    real(dl) lnrat
    !ScalarPower = const for scale invariant spectrum
    !The normalization is defined so that for adiabatic perturbations the gradient of the 3-Ricci
    !scalar on co-moving hypersurfaces receives power
    ! < |D_a R^{(3)}|^2 > = int dk/k 16 k^6/S^6 (1-3K/k^2)^2 ScalarPower(k)
    !In other words ScalarPower is the power spectrum of the conserved curvature perturbation given by
    !-chi = \Phi + 2/3*\Omega^{-1} \frac{H^{-1}\Phi' - \Psi}{1+w}
    !(w=p/\rho), so < |\chi(x)|^2 > = \int dk/k ScalarPower(k).
    !Near the end of inflation chi is equal to 3/2 Psi.
    !Here nu^2 = (k^2 + _curv)/|_curv|

    !This power spectrum is also used for isocurvature modes where
    !< |\Delta(x)|^2 > = \int dk/k ScalarPower(k)
    !For the isocurvture velocity mode ScalarPower is the power in the neutrino heat flux.


    lnrat = log(k/this%pivot_scalar)
    TInitialPowerLaw_ScalarPower = this%As * exp(lnrat * (this%ns - 1 + &
        &             lnrat * (this%nrun / 2 + this%nrunrun / 6 * lnrat)))

    end function TInitialPowerLaw_ScalarPower


    end module InitialPower
