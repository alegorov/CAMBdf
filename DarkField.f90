    module DarkFieldModule
    use precision
    use classes
    use Interpolation
    implicit none

    private

    type :: TDarkFieldSplines
        type(TCubicSpline) :: dadtau, W, N
    contains
    procedure :: Clear
    end type TDarkFieldSplines

    type, extends(TCambComponent) :: TDarkFieldModel
        type(TDarkFieldSplines), allocatable :: splines
        real(dl) :: delta
        real(dl) :: W0
        real(dl) :: p(66)

    contains
    procedure :: Init
    procedure :: dfDeriv
    procedure :: a2X
    procedure :: ia2M
    end type TDarkFieldModel

    public TDarkFieldSplines
    public TDarkFieldModel
    contains

    subroutine Clear(this)
    class(TDarkFieldSplines), intent(inout) :: this

    call this%dadtau%Clear()
    call this%W%Clear()
    call this%N%Clear()

    end subroutine Clear

    subroutine Init(this)
    class(TDarkFieldModel), intent(inout) :: this
    real(dl) :: alpha,delta

    delta = this%delta
    alpha = -7 - 3*delta

    this%p = 0._dl

    this%p(1) = -2
    this%p(2) = delta + 2
    this%p(3) = -8
    this%p(4) = 8
    this%p(5) = -(8*delta + 19)/3
    this%p(6) = (delta + 2)/3
    this%p(7) = 4 + delta
    this%p(8) = -20
    this%p(9) = 13._dl/4
    this%p(10) = -2*(7*delta**2 + 33*delta + 39)/alpha
    this%p(11) = 4*(8 + 3*delta)*(delta + 2)/alpha
    this%p(12) = -2
    this%p(13) = 4*(delta + 2)/alpha
    this%p(14) = -8/(3*alpha)*(delta + 2)
    this%p(15) = (8 + 3*delta)/alpha
    this%p(16) = (3*delta + 10)/alpha
    this%p(17) = (5 + delta)/(3*alpha)
    this%p(18) = 4*(delta + 1)/alpha
    this%p(19) = 2*(21*delta + 68)/alpha
    this%p(20) = -2*(3*delta + 10)/alpha
    this%p(21) = -6*(8 + 3*delta)*(7*delta + 16)/alpha**2
    this%p(22) = 6*(8 + 3*delta)*(9*delta + 22)/alpha**2
    this%p(23) = -2*(7*delta**2 + 27*delta + 24)/alpha**2
    this%p(24) = 6*this%p(21)
    this%p(25) = 6*this%p(22)
    this%p(26) = 4*(8 + 3*delta)*(5 + delta)/alpha**2
    this%p(27) = -16*(delta + 2)/(3*(alpha + delta + 1))
    this%p(28) = 8*(delta + 2)/(3*(alpha + delta + 1))
    this%p(29) = -2
    this%p(30) = -16*delta/(3*(alpha + delta + 1))
    this%p(31) = -4._dl/3
    this%p(32) = 2*(delta + 8)/(alpha + delta + 1)
    this%p(33) = 2*(delta - 2)/(alpha + delta + 1)
    this%p(34) = -1
    this%p(35) = -4*(7*delta**2 + 32*delta + 37)/(alpha*(alpha + delta + 1))
    this%p(36) = 8*(delta + 2)*(11 + 3*delta)/(alpha*(alpha + delta + 1))
    this%p(37) = 16*(7*delta**2 + 33*delta + 40)/(alpha*(alpha + delta + 1))
    this%p(38) = -16*(delta + 2)**2/(alpha*(alpha + delta + 1))
    this%p(39) = -32*(8 + 3*delta)*(7*delta + 13)/(alpha*(alpha + delta + 1))
    this%p(40) = 32*(8 + 3*delta)*(delta + 2)/(alpha*(alpha + delta + 1))
    this%p(41) = 8
    this%p(42) = -8*(3 + delta)
    this%p(43) = -320
    this%p(44) = 64 + (64*delta)/3
    this%p(45) = 64
    this%p(46) = -(8*delta)/3
    this%p(47) = -8*(delta + 8)
    this%p(48) = 8
    this%p(49) = -80
    this%p(50) = 640
    this%p(51) = 8/(3*alpha)*(7*alpha**2 + 7*alpha*delta + alpha - 2*delta)
    this%p(52) = -16
    this%p(53) = 32*delta*(1 - 1/alpha)
    this%p(54) = 64._dl/3
    this%p(55) = 8*(delta + 8._dl/3)
    this%p(56) = -128._dl/3
    this%p(57) = 40._dl/3
    this%p(58) = 8*(delta - 4._dl/3)
    this%p(59) = 16._dl/3
    this%p(60) = 8*(3 + delta)
    this%p(61) = 16._dl/3*(1/alpha - 8)
    this%p(62) = 32*(1/alpha - 3)
    this%p(63) = 16/(3*alpha)*(-7*alpha**2 + 106*alpha - 19)
    this%p(64) = 8*(alpha**2 + alpha*delta - 9*alpha + 2)/alpha
    this%p(65) = 32*(1/alpha - 1)*(7*alpha - 19)
    this%p(66) = 32*(1/alpha - 1)*(3 - alpha)

    end subroutine Init

    ! Dark Field equation of motion
    subroutine dfDeriv(this,ix_df,ix_etak,ay,ayprime,k,k2,W,N,B,dB,B2,BN,kN,kB,k2N)
    class(TDarkFieldModel), intent(inout) :: this
    integer ix_df,ix_etak
    real(dl), target :: ay(:),ayprime(:)
    real(dl) :: k,k2,W,N,B,dB,B2,BN,kN,kB,k2N
    real(dl) :: x,y,m,dx,dy,dm
    real(dl) :: h,eta,dh,deta
    real(dl) :: p1,p2,p3,p4,p5,p6,p7,p8,p9,p10
    real(dl) :: p11,p12,p13,p14,p15,p16,p17,p18,p19,p20
    real(dl) :: p21,p22,p23,p24,p25,p26,p27,p28,p29,p30
    real(dl) :: p31,p32,p33,p34,p35,p36,p37,p38,p39,p40

    x = ay(ix_df + 1)
    y = ay(ix_df + 2)
    m = ay(ix_df + 3)
    dx = ay(ix_df + 4)
    dy = ay(ix_df + 5)
    dm = ay(ix_df + 6)

    h = ay(ix_df)
    eta = ay(ix_etak)/k

    dh = ayprime(ix_df)
    deta = ayprime(ix_etak)/k

    p1 = this%p(1)
    p2 = this%p(2)
    p3 = this%p(3)
    p4 = this%p(4)
    p5 = this%p(5)
    p6 = this%p(6)
    p7 = this%p(7)
    p8 = this%p(8)
    p9 = this%p(9)
    p10 = this%p(10)
    p11 = this%p(11)
    p12 = this%p(12)
    p13 = this%p(13)
    p14 = this%p(14)
    p15 = this%p(15)
    p16 = this%p(16)
    p17 = this%p(17)
    p18 = this%p(18)
    p19 = this%p(19)
    p20 = this%p(20)
    p21 = this%p(21)
    p22 = this%p(22)
    p23 = this%p(23)
    p24 = this%p(24)
    p25 = this%p(25)
    p26 = this%p(26)
    p27 = this%p(27)
    p28 = this%p(28)
    p29 = this%p(29)
    p30 = this%p(30)
    p31 = this%p(31)
    p32 = this%p(32)
    p33 = this%p(33)
    p34 = this%p(34)
    p35 = this%p(35)
    p36 = this%p(36)
    p37 = this%p(37)
    p38 = this%p(38)
    p39 = this%p(39)
    p40 = this%p(40)

    ayprime(ix_df + 1) = dx
    ayprime(ix_df + 2) = dy
    ayprime(ix_df + 3) = dm
    ayprime(ix_df + 4) = p1*B*dx + p2*k*dm + (p3*dB + p4*B2 + p5*k2)*x + p6*k2*y + p7*kB*m + (p8*BN + p9*W)*dh + k2N*(p10*h + p11*eta)
    ayprime(ix_df + 5) = p12*B*dy + p13*k*dm + p14*k2*x + (p15*dB + p16*B2 + p17*k2)*y + p18*kB*m + (p19*BN + p20*W)*(dh + 6*deta) + N*((p21*dB + p22*B2 + p23*k2)*h + (p24*dB + p25*B2 + p26*k2)*eta)
    ayprime(ix_df + 6) = k*(p27*dx + p28*dy) + B*(p29*dm + k*(p30*x + p31*y)) + (p32*dB + p33*B2 + p34*k2)*m + kN*(p35*dh + p36*deta) + k*((p37*BN + p38*W)*h + (p39*BN + p40*W)*eta)


    end subroutine dfDeriv

    function a2X(this,ix_df,ix_etak,ay,k,k2,W,N,B2,W2,N2,kN,BW,BN,k2N,k2N2)
    class(TDarkFieldModel), intent(inout) :: this
    integer ix_df,ix_etak
    real(dl), target :: ay(:)
    real(dl) :: k,k2,W,N,B2,W2,N2,kN,BW,BN,k2N,k2N2
    real(dl) :: a2X
    real(dl) :: x,y,m,dx,dm
    real(dl) :: h,eta
    real(dl) :: p41,p42,p43,p44,p45,p46,p47,p48,p50,p51,p52,p53

    x = ay(ix_df + 1)
    y = ay(ix_df + 2)
    m = ay(ix_df + 3)
    dx = ay(ix_df + 4)
    !dy = ay(ix_df + 5)
    dm = ay(ix_df + 6)

    h = ay(ix_df)
    eta = ay(ix_etak)/k

    p41 = this%p(41)
    p42 = this%p(42)
    p43 = this%p(43)
    p44 = this%p(44)
    p45 = this%p(45)
    p46 = this%p(46)
    p47 = this%p(47)
    p48 = this%p(48)
    !p49 = this%p(49)
    p50 = this%p(50)
    p51 = this%p(51)
    p52 = this%p(52)
    p53 = this%p(53)

    a2X = p41*W*dx + p42*kN*dm + ((p43*B2 + p44*k2)*N + p45*BW)*x + p46*k2N*y + (p47*BN + p48*W)*k*m + ((p50*B2 + p51*k2)*N2 + p52*W2)*h + p53*k2N2*eta

    end function a2X

    function ia2M(this,ix_df,ix_etak,ay,ayprime,k,W,N,N2,BN,kN)
    class(TDarkFieldModel), intent(inout) :: this
    integer ix_df,ix_etak
    real(dl), target :: ay(:),ayprime(:)
    real(dl) :: k,W,N,N2,BN,kN
    real(dl) :: ia2M
    real(dl) :: x,y,m,dx,dy
    real(dl) :: h,eta,dh
    real(dl) :: p54,p55,p56,p57,p58,p59,p60,p61,p63,p64,p65,p66

    x = ay(ix_df + 1)
    y = ay(ix_df + 2)
    m = ay(ix_df + 3)
    dx = ay(ix_df + 4)
    dy = ay(ix_df + 5)
    !dm = ay(ix_df + 6)

    h = ay(ix_df)
    eta = ay(ix_etak)/k

    dh = ayprime(ix_df)

    p54 = this%p(54)
    p55 = this%p(55)
    p56 = this%p(56)
    p57 = this%p(57)
    p58 = this%p(58)
    p59 = this%p(59)
    p60 = this%p(60)
    p61 = this%p(61)
    !p62 = this%p(62)
    p63 = this%p(63)
    p64 = this%p(64)
    p65 = this%p(65)
    p66 = this%p(66)

    ia2M = N*(p54*dx + p55*dy) + (p56*BN + p57*W)*x + (p58*BN + p59*W)*y + p60*kN*m + N2*(p61*dh) + N*((p63*BN + p64*W)*h + (p65*BN + p66*W)*eta)

    end function ia2M

    end module DarkFieldModule
