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
        real(dl) :: delta,epsilon,Lambda
        real(dl) :: W0
        real(dl) :: p(99)

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
    real(dl) :: alpha,delta,epsilon,Lambda

    delta = this%delta
    epsilon = this%epsilon
    Lambda = this%Lambda

    alpha = (12*delta**2 + (48 - 12*Lambda + 12*epsilon)*delta + 48 + 24*epsilon - 28*Lambda + 3*epsilon**2)/(4*Lambda)

    this%p = 0._dl

    this%p(1) = (8*Lambda*alpha - 12*Lambda*epsilon + 24*Lambda - 18*delta - 9*epsilon - 36)/(12 + 6*delta - 8*Lambda + 3*epsilon)
    this%p(2) = 3*(4*Lambda - 3)*(4 + 2*delta + epsilon)/(12 + 6*delta - 8*Lambda + 3*epsilon)
    this%p(3) = (epsilon - alpha + delta + 1)/2
    this%p(4) = 6*(epsilon + 2*delta)*(epsilon + 2*delta - 2*Lambda + 4)/(12 + 6*delta - 8*Lambda + 3*epsilon)
    this%p(5) = 12*(2*Lambda - epsilon - 2*delta - 4)*(alpha + delta - epsilon + 1)/(12 + 6*delta - 8*Lambda + 3*epsilon)
    this%p(6) = 2*(2*Lambda - epsilon - 2*delta - 4)*(epsilon - 2*alpha)/(12 + 6*delta - 8*Lambda + 3*epsilon)
    this%p(7) = -epsilon/6
    this%p(8) = 12*Lambda*(epsilon + 2*delta)/(12 + 6*delta - 8*Lambda + 3*epsilon)
    this%p(9) = 24*Lambda*(epsilon - alpha - delta - 1)/(12 + 6*delta - 8*Lambda + 3*epsilon)
    this%p(10) = (8*Lambda*alpha - 4*Lambda*epsilon + 8*Lambda - 6*delta - 3*epsilon - 12)/(12 + 6*delta - 8*Lambda + 3*epsilon)
    this%p(11) = (3*epsilon - 3*alpha - delta - 1)/2
    this%p(12) = -20
    this%p(13) = (5 - 4*delta - 2*epsilon)/4
    this%p(14) = 20*(4 + 2*delta + epsilon)
    this%p(15) = (-14*alpha**2 + 7*alpha*epsilon - 2*alpha - 2*epsilon)/(6*alpha)
    this%p(16) = alpha - 5*delta - 4*epsilon - 9
    this%p(17) = 2*epsilon*(1 - 1/alpha)
    this%p(18) = -2
    this%p(19) = (1 - alpha + delta)/alpha
    this%p(20) = 2*(8*Lambda*delta - 8*delta**2 - 4*delta*epsilon + 16*Lambda - 30*delta - 7*epsilon - 28)/(alpha*(12 + 6*delta - 8*Lambda + 3*epsilon))
    this%p(21) = 1/alpha - 1
    this%p(22) = 3/alpha - 1
    this%p(23) = -(2 + alpha + 2*delta)/(3*alpha)
    this%p(24) = 2*(-8*Lambda*delta - 16*Lambda + 6*delta + 3*epsilon + 12)/(alpha*(12 + 6*delta - 8*Lambda + 3*epsilon))
    this%p(25) = 4*(delta + 1)/alpha
    this%p(26) = 2*(19/alpha - 7)
    this%p(27) = 2*(1 - 3/alpha)
    this%p(28) = 2*(1/alpha - 1)*(1/alpha + 7)
    this%p(29) = 6*(1 - 1/alpha)*(3 - 1/alpha)
    this%p(30) = 2*(7*alpha*delta + 10*alpha - 2*delta - 2)/(3*alpha**2)
    this%p(31) = 6*this%p(28)
    this%p(32) = 6*this%p(29)
    this%p(33) = 4*(1 - 1/alpha)*(2 + alpha + 2*delta)/alpha
    this%p(34) = (12*delta + 12*epsilon + 21 - 12*Lambda - 4*alpha)*(4 + 2*delta + epsilon)/(12 + 6*delta - 8*Lambda + 3*epsilon)
    this%p(35) = (-8*Lambda*alpha + 12*Lambda*epsilon + 8*Lambda - 6*delta - 3*epsilon - 12)/(12 + 6*delta - 8*Lambda + 3*epsilon)
    this%p(36) = (1 - alpha + delta + 3*epsilon)/6
    this%p(37) = 2*(epsilon + 2*delta - 2*Lambda + 4)*(3*epsilon - 2*alpha + 2)/(12 + 6*delta - 8*Lambda + 3*epsilon)
    this%p(38) = 4*(epsilon + 2*delta - 2*Lambda + 4)*(3 + 3*epsilon - alpha + 3*delta)/(12 + 6*delta - 8*Lambda + 3*epsilon)
    this%p(39) = (8*Lambda*delta + 12*Lambda*epsilon - 8*delta**2 - 16*delta*epsilon - 6*epsilon**2 + 16*Lambda - 30*delta - 31*epsilon - 28)/(3*(12 + 6*delta - 8*Lambda + 3*epsilon))
    this%p(40) = (2*alpha - 2*delta - 2 - 3*epsilon)/18
    this%p(41) = 4*Lambda*(3*epsilon - 2*alpha + 2)/(12 + 6*delta - 8*Lambda + 3*epsilon)
    this%p(42) = 8*Lambda*(3 + 3*epsilon - alpha + 3*delta)/(12 + 6*delta - 8*Lambda + 3*epsilon)
    this%p(43) = 2*(-4*Lambda*delta - 6*Lambda*epsilon + 4*Lambda - 6*delta - 3*epsilon - 12)/(3*(12 + 6*delta - 8*Lambda + 3*epsilon))
    this%p(44) = (7 - 3*alpha + 7*delta + 9*epsilon)/6
    this%p(45) = 20._dl/3
    this%p(46) = (4*alpha - 6*epsilon - 9)/12
    this%p(47) = 20*(3*epsilon - 2*alpha - 2)/3
    this%p(48) = (14*alpha*delta + 21*alpha*epsilon + 20*alpha - 4*delta - 6*epsilon - 4)/(18*alpha)
    this%p(49) = (1 - 3*delta - 12*epsilon + 7*alpha)/3
    this%p(50) = (2._dl/3)*(1 - 1/alpha)*(2 - 2*alpha + 2*delta + 3*epsilon)
    this%p(51) = 4*(-4*Lambda*alpha - 4*Lambda*delta + 4*alpha*delta + 2*alpha*epsilon + 4*delta**2 + 2*delta*epsilon - 12*Lambda + 8*alpha + 18*delta + 5*epsilon + 20)/((alpha + delta + 1)*(12 + 6*delta - 8*Lambda + 3*epsilon))
    this%p(52) = 2*(1 - alpha + delta)/(3*(alpha + delta + 1))
    this%p(53) = 4*(4*Lambda*alpha + 4*Lambda*delta + 12*Lambda - 6*delta - 3*epsilon - 12)/((alpha + delta + 1)*(12 + 6*delta - 8*Lambda + 3*epsilon))
    this%p(54) = -2
    this%p(55) = 4*(4*Lambda*alpha + 20*Lambda*delta - 4*alpha*delta - 2*alpha*epsilon - 20*delta**2 - 10*delta*epsilon + 28*Lambda - 8*alpha - 66*delta - 13*epsilon - 52)/((alpha + delta + 1)*(12 + 6*delta - 8*Lambda + 3*epsilon))
    this%p(56) = -4._dl/3
    this%p(57) = 4*(-4*Lambda*alpha - 20*Lambda*delta - 28*Lambda + 6*delta + 3*epsilon + 12)/((alpha + delta + 1)*(12 + 6*delta - 8*Lambda + 3*epsilon))
    this%p(58) = 2*(1 - alpha - 2*delta)/(alpha + delta + 1)
    this%p(59) = 2*(alpha + 4*delta + 5)/(alpha + delta + 1)
    this%p(60) = -1
    this%p(61) = 2*(-7*alpha**2 - 7*alpha*delta - 17*alpha + 2*delta + 2)/(3*alpha*(alpha + delta + 1))
    this%p(62) = 8*(-alpha*delta - 3*alpha + delta + 1)/(alpha*(alpha + delta + 1))
    this%p(63) = 4*(21*alpha**2 + 35*alpha*delta + 77*alpha - 10*delta - 10)/(3*alpha*(alpha + delta + 1))
    this%p(64) = 4*(-2*alpha**2 - 2*alpha*delta - 7*alpha + delta + 1)/(3*alpha*(alpha + delta + 1))
    this%p(65) = 16*(-1 + alpha)*(5 - 3*alpha + 5*delta)/(alpha*(alpha + delta + 1))
    this%p(66) = 8*(-1 + alpha)*(alpha - delta - 1)/(alpha*(alpha + delta + 1))
    this%p(67) = (24*epsilon + 48*delta - 48*Lambda + 96)/(12 + 6*delta - 8*Lambda + 3*epsilon)
    this%p(68) = 48*Lambda/(12 + 6*delta - 8*Lambda + 3*epsilon)
    this%p(69) = 4*(alpha + delta + 1)
    this%p(70) = 960*(2*Lambda - epsilon - 2*delta - 4)/(12 + 6*delta - 8*Lambda + 3*epsilon)
    this%p(71) = 8*(8*Lambda*alpha + 8*Lambda*delta - 8*alpha*delta - 4*alpha*epsilon - 8*delta**2 - 4*delta*epsilon + 8*Lambda - 16*alpha - 22*delta - 3*epsilon - 12)/(12 + 6*delta - 8*Lambda + 3*epsilon)
    this%p(72) = 192*(epsilon + 2*delta - 2*Lambda + 4)/(12 + 6*delta - 8*Lambda + 3*epsilon)
    this%p(73) = -(8*delta)/3
    this%p(74) = (-1920*Lambda)/(12 + 6*delta - 8*Lambda + 3*epsilon)
    this%p(75) = 8*(-8*Lambda*alpha - 8*Lambda*delta - 8*Lambda + 6*delta + 3*epsilon + 12)/(12 + 6*delta - 8*Lambda + 3*epsilon)
    this%p(76) = 384*Lambda/(12 + 6*delta - 8*Lambda + 3*epsilon)
    this%p(77) = 8*(alpha + 2*delta - 1)
    this%p(78) = 8
    this%p(79) = -80
    this%p(80) = 640
    this%p(81) = (8/(3*alpha))*(7*alpha**2 + 7*alpha*delta + alpha - 2*delta)
    this%p(82) = -16
    this%p(83) = 32*delta*(1 - 1/alpha)
    this%p(84) = (56*epsilon + 112*delta - 128*Lambda + 224)/(12 + 6*delta - 8*Lambda + 3*epsilon)
    this%p(85) = (8._dl/3)*(1 - alpha)
    this%p(86) = (128*Lambda - 24*epsilon - 48*delta - 96)/(12 + 6*delta - 8*Lambda + 3*epsilon)
    this%p(87) = 32*(2*Lambda - epsilon - 2*delta - 4)*(alpha + 3*delta + 11)/(12 + 6*delta - 8*Lambda + 3*epsilon)
    this%p(88) = 40*(epsilon + 2*delta - 2*Lambda + 4)/(12 + 6*delta - 8*Lambda + 3*epsilon)
    this%p(89) = -(8._dl/3)*(alpha + 11)
    this%p(90) = 16._dl/3
    this%p(91) = -64*Lambda*(alpha + 3*delta + 11)/(12 + 6*delta - 8*Lambda + 3*epsilon)
    this%p(92) = 80*Lambda/(12 + 6*delta - 8*Lambda + 3*epsilon)
    this%p(93) = -4*(alpha + delta + 1)
    this%p(94) = (16._dl/3)*(1/alpha - 8)
    this%p(95) = 32*(1/alpha - 3)
    this%p(96) = (16/(3*alpha))*(-7*alpha**2 + 106*alpha - 19)
    this%p(97) = 8*(alpha**2 + alpha*delta - 9*alpha + 2)/alpha
    this%p(98) = 32*(1/alpha - 1)*(7*alpha - 19)
    this%p(99) = 32*(1/alpha - 1)*(3 - alpha)

    end subroutine Init

    ! Dark Field equation of motion
    subroutine dfDeriv(this,ix_df,ix_etak,ay,ayprime,k,k2,W,N,B,dB,B2,BW,BN,kN,kB,k2N)
    class(TDarkFieldModel), intent(inout) :: this
    integer ix_df,ix_etak
    real(dl), target :: ay(:),ayprime(:)
    real(dl) :: k,k2,W,N,B,dB,B2,BW,BN,kN,kB,k2N
    real(dl) :: x,y,z,m,dx,dy,dz,dm
    real(dl) :: h,eta,dh,deta
    real(dl) :: p1,p2,p3,p4,p5,p6,p7,p8,p9,p10
    real(dl) :: p11,p12,p13,p14,p15,p16,p17,p18,p19,p20
    real(dl) :: p21,p22,p23,p24,p25,p26,p27,p28,p29,p30
    real(dl) :: p31,p32,p33,p34,p35,p36,p37,p38,p39,p40
    real(dl) :: p41,p42,p43,p44,p45,p46,p47,p48,p49,p50
    real(dl) :: p51,p52,p53,p54,p55,p56,p57,p58,p59,p60
    real(dl) :: p61,p62,p63,p64,p65,p66

    x = ay(ix_df + 1)
    y = ay(ix_df + 2)
    z = ay(ix_df + 3)
    m = ay(ix_df + 4)
    dx = ay(ix_df + 5)
    dy = ay(ix_df + 6)
    dz = ay(ix_df + 7)
    dm = ay(ix_df + 8)

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
    p41 = this%p(41)
    p42 = this%p(42)
    p43 = this%p(43)
    p44 = this%p(44)
    p45 = this%p(45)
    p46 = this%p(46)
    p47 = this%p(47)
    p48 = this%p(48)
    p49 = this%p(49)
    p50 = this%p(50)
    p51 = this%p(51)
    p52 = this%p(52)
    p53 = this%p(53)
    p54 = this%p(54)
    p55 = this%p(55)
    p56 = this%p(56)
    p57 = this%p(57)
    p58 = this%p(58)
    p59 = this%p(59)
    p60 = this%p(60)
    p61 = this%p(61)
    p62 = this%p(62)
    p63 = this%p(63)
    p64 = this%p(64)
    p65 = this%p(65)
    p66 = this%p(66)

    ayprime(ix_df + 1) = dx
    ayprime(ix_df + 2) = dy
    ayprime(ix_df + 3) = dz
    ayprime(ix_df + 4) = dm
    ayprime(ix_df + 5) = B*(p1*dx + p2*dz) + p3*k*dm + (p4*dB + p5*B2 + p6*k2)*x + p7*k2*y + (p8*dB + p9*B2 + p10*k2)*z + p11*kB*m + (p12*BN + p13*W)*dh + ((p14*B2 + p15*k2)*N + p16*BW)*h + p17*k2N*eta
    ayprime(ix_df + 6) = p18*B*dy + p19*k*dm + p20*k2*x + (p21*dB + p22*B2 + p23*k2)*y + p24*k2*z + p25*kB*m + (p26*BN + p27*W)*(dh + 6*deta) + N*((p28*dB + p29*B2 + p30*k2)*h + (p31*dB + p32*B2 + p33*k2)*eta)
    ayprime(ix_df + 7) = B*(p34*dx + p35*dz) + p36*k*dm + (p37*dB + p38*B2 + p39*k2)*x + p40*k2*y + (p41*dB + p42*B2 + p43*k2)*z + p44*kB*m + (p45*BN + p46*W)*dh + ((p47*B2 + p48*k2)*N + p49*BW)*h + p50*k2N*eta
    ayprime(ix_df + 8) = k*(p51*dx + p52*dy + p53*dz) + B*(p54*dm + k*(p55*x + p56*y + p57*z)) + (p58*dB + p59*B2 + p60*k2)*m + kN*(p61*dh + p62*deta) + k*((p63*BN + p64*W)*h + (p65*BN + p66*W)*eta)


    end subroutine dfDeriv

    function a2X(this,ix_df,ix_etak,ay,k,k2,W,N,B2,W2,N2,kN,BW,BN,k2N,k2N2)
    class(TDarkFieldModel), intent(inout) :: this
    integer ix_df,ix_etak
    real(dl), target :: ay(:)
    real(dl) :: k,k2,W,N,B2,W2,N2,kN,BW,BN,k2N,k2N2
    real(dl) :: a2X
    real(dl) :: x,y,z,m,dx,dz,dm
    real(dl) :: h,eta
    real(dl) :: p67,p68,p69,p70,p71,p72,p73,p74,p75,p76,p77,p78,p80,p81,p82,p83

    x = ay(ix_df + 1)
    y = ay(ix_df + 2)
    z = ay(ix_df + 3)
    m = ay(ix_df + 4)
    dx = ay(ix_df + 5)
    !dy = ay(ix_df + 6)
    dz = ay(ix_df + 7)
    dm = ay(ix_df + 8)

    h = ay(ix_df)
    eta = ay(ix_etak)/k

    p67 = this%p(67)
    p68 = this%p(68)
    p69 = this%p(69)
    p70 = this%p(70)
    p71 = this%p(71)
    p72 = this%p(72)
    p73 = this%p(73)
    p74 = this%p(74)
    p75 = this%p(75)
    p76 = this%p(76)
    p77 = this%p(77)
    p78 = this%p(78)
    !p79 = this%p(79)
    p80 = this%p(80)
    p81 = this%p(81)
    p82 = this%p(82)
    p83 = this%p(83)

    a2X = W*(p67*dx + p68*dz) + p69*kN*dm + ((p70*B2 + p71*k2)*N + p72*BW)*x + p73*k2N*y + ((p74*B2 + p75*k2)*N + p76*BW)*z + (p77*BN + p78*W)*k*m + ((p80*B2 + p81*k2)*N2 + p82*W2)*h + p83*k2N2*eta

    end function a2X

    function ia2M(this,ix_df,ix_etak,ay,ayprime,k,W,N,N2,BN,kN)
    class(TDarkFieldModel), intent(inout) :: this
    integer ix_df,ix_etak
    real(dl), target :: ay(:),ayprime(:)
    real(dl) :: k,W,N,N2,BN,kN
    real(dl) :: ia2M
    real(dl) :: x,y,z,m,dx,dy,dz
    real(dl) :: h,eta,dh
    real(dl) :: p84,p85,p86,p87,p88,p89,p90,p91,p92,p93,p94,p96,p97,p98,p99

    x = ay(ix_df + 1)
    y = ay(ix_df + 2)
    z = ay(ix_df + 3)
    m = ay(ix_df + 4)
    dx = ay(ix_df + 5)
    dy = ay(ix_df + 6)
    dz = ay(ix_df + 7)
    !dm = ay(ix_df + 8)

    h = ay(ix_df)
    eta = ay(ix_etak)/k

    dh = ayprime(ix_df)

    p84 = this%p(84)
    p85 = this%p(85)
    p86 = this%p(86)
    p87 = this%p(87)
    p88 = this%p(88)
    p89 = this%p(89)
    p90 = this%p(90)
    p91 = this%p(91)
    p92 = this%p(92)
    p93 = this%p(93)
    p94 = this%p(94)
    !p95 = this%p(95)
    p96 = this%p(96)
    p97 = this%p(97)
    p98 = this%p(98)
    p99 = this%p(99)

    ia2M = N*(p84*dx + p85*dy + p86*dz) + (p87*BN + p88*W)*x + (p89*BN + p90*W)*y + (p91*BN + p92*W)*z + p93*kN*m + N2*(p94*dh) + N*((p96*BN + p97*W)*h + (p98*BN + p99*W)*eta)

    end function ia2M

    end module DarkFieldModule
