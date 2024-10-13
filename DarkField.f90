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
        real(dl) :: alpha,delta
        real(dl) :: W0
        real(dl) :: p(88)

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

    alpha = this%alpha
    delta = this%delta

    this%p = 0._dl

    this%p(1) = -3
    this%p(2) = -(3 + delta + alpha)/2
    this%p(3) = -8
    this%p(4) = -4*(5 + 3*delta + alpha)
    this%p(5) = 4._dl/3*(alpha + delta + 2)
    this%p(6) = (2 + delta)/3
    this%p(7) = -1
    this%p(8) = -(7*delta + 3*alpha + 13)/2
    this%p(9) = -20
    this%p(10) = 13._dl/4
    this%p(11) = (-7*alpha**2 - 7*alpha*delta - 15*alpha + 2*delta + 4)/(3*alpha)
    this%p(12) = 7 + alpha + 3*delta
    this%p(13) = 4*(2 + delta)*(1/alpha - 1)
    this%p(14) = -2
    this%p(15) = (1 - alpha + delta)/alpha
    this%p(16) = -2*(4*delta + 7)/(3*alpha)
    this%p(17) = 1/alpha - 1
    this%p(18) = 3/alpha - 1
    this%p(19) = -(2 + alpha + 2*delta)/(3*alpha)
    this%p(20) = 2/alpha
    this%p(21) = 4*(1 + delta)/alpha
    this%p(22) = 2*(19/alpha - 7)
    this%p(23) = 2*(1 - 3/alpha)
    this%p(24) = 2*(1/alpha - 1)*(1/alpha + 7)
    this%p(25) = 6*(-1/alpha + 1)*(3 - 1/alpha)
    this%p(26) = 2*(7*alpha*delta + 10*alpha - 2*delta - 2)/(3*alpha**2)
    this%p(27) = 6*this%p(24)
    this%p(28) = 6*this%p(25)
    this%p(29) = 4*(-1/alpha + 1)*(2 + alpha + 2*delta)/alpha
    this%p(30) = -(27 + 12*delta + 4*alpha)/3
    this%p(31) = -1
    this%p(32) = -(11 + 5*delta + alpha)/6
    this%p(33) = -4._dl/3*(5 + 3*delta + alpha)
    this%p(34) = -4._dl/3*(alpha + 3*delta + 9)
    this%p(35) = (8*delta + 17)/9
    this%p(36) = (alpha + 2*delta + 5)/9
    this%p(37) = -2._dl/3
    this%p(38) = -(3*alpha + 29 + 11*delta)/6
    this%p(39) = 20._dl/3
    this%p(40) = (4*alpha + 15 + 12*delta)/12
    this%p(41) = -40._dl/3*(7 + alpha + 3*delta)
    this%p(42) = 2/(9*alpha)*(-7*alpha*delta - 16*alpha + 2*delta + 5)
    this%p(43) = 7._dl/3*(7 + alpha + 3*delta)
    this%p(44) = 4._dl/3*(1/alpha - 1)*(alpha + 2*delta + 5)
    this%p(45) = 4*(5 + 2*alpha + 2*delta)/(3*(alpha + delta + 1))
    this%p(46) = 2*(1 - alpha + delta)/(3*(alpha + delta + 1))
    this%p(47) = -4/(alpha + delta + 1)
    this%p(48) = -2
    this%p(49) = -4*(2*alpha + 10*delta + 13)/(3*(alpha + delta + 1))
    this%p(50) = -4._dl/3
    this%p(51) = 4/(alpha + delta + 1)
    this%p(52) = 2*(1 - alpha - 2*delta)/(alpha + delta + 1)
    this%p(53) = 2*(alpha + 4*delta + 5)/(alpha + delta + 1)
    this%p(54) = -1
    this%p(55) = 2*(-7*alpha**2 - 7*alpha*delta - 17*alpha + 2*delta + 2)/(3*alpha*(alpha + delta + 1))
    this%p(56) = 8*(-alpha*delta - 3*alpha + delta + 1)/(alpha*(alpha + delta + 1))
    this%p(57) = 4*(21*alpha**2 + 35*alpha*delta + 77*alpha - 10*delta - 10)/(3*alpha*(alpha + delta + 1))
    this%p(58) = 4*(-2*alpha**2 - 2*alpha*delta - 7*alpha + delta + 1)/(3*alpha*(alpha + delta + 1))
    this%p(59) = 16*(-1 + alpha)*(5 - 3*alpha + 5*delta)/(alpha*(alpha + delta + 1))
    this%p(60) = 8*(-1 + alpha)*(-1 + alpha - delta)/(alpha*(alpha + delta + 1))
    this%p(61) = 8
    this%p(62) = 4*(alpha + delta + 1)
    this%p(63) = -320
    this%p(64) = -8._dl/3*(4*alpha + 4*delta + 3)
    this%p(65) = 64
    this%p(66) = -(8*delta)/3
    this%p(67) = 8
    this%p(68) = 8*(-1 + alpha + 2*delta)
    this%p(69) = 8
    this%p(70) = -80
    this%p(71) = 640
    this%p(72) = 8/(3*alpha)*(7*alpha**2 + 7*alpha*delta + alpha - 2*delta)
    this%p(73) = -16
    this%p(74) = 32*delta*(-1/alpha + 1)
    this%p(75) = 56._dl/3
    this%p(76) = 8._dl/3*(1 - alpha)
    this%p(77) = -8
    this%p(78) = -32._dl/3*(11 + alpha + 3*delta)
    this%p(79) = 40._dl/3
    this%p(80) = -8._dl/3*(alpha + 11)
    this%p(81) = 16._dl/3
    this%p(82) = -4*(alpha + delta + 1)
    this%p(83) = 16._dl/3*(1/alpha - 8)
    this%p(84) = 32*(1/alpha - 3)
    this%p(85) = 16/(3*alpha)*(-7*alpha**2 + 106*alpha - 19)
    this%p(86) = 8*(alpha**2 + alpha*delta - 9*alpha + 2)/alpha
    this%p(87) = 32*(1/alpha - 1)*(7*alpha - 19)
    this%p(88) = 32*(1/alpha - 1)*(3 - alpha)

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

    ayprime(ix_df + 1) = dx
    ayprime(ix_df + 2) = dy
    ayprime(ix_df + 3) = dz
    ayprime(ix_df + 4) = dm
    ayprime(ix_df + 5) = p1*B*(dx + dz) + p2*k*dm + (p3*dB + p4*B2 + p5*k2)*x + k2*(p6*y + p7*z) + p8*kB*m + (p9*BN + p10*W)*dh + (p11*k2N + p12*BW)*h + p13*k2N*eta
    ayprime(ix_df + 6) = p14*B*dy + p15*k*dm + p16*k2*x + (p17*dB + p18*B2 + p19*k2)*y + p20*k2*z + p21*kB*m + (p22*BN + p23*W)*(dh + 6*deta) + N*((p24*dB + p25*B2 + p26*k2)*h + (p27*dB + p28*B2 + p29*k2)*eta)
    ayprime(ix_df + 7) = B*(p30*dx + p31*dz) + p32*k*dm + (p33*dB + p34*B2 + p35*k2)*x + k2*(p36*y + p37*z) + p38*kB*m + (p39*BN + p40*W)*dh + ((p41*B2 + p42*k2)*N + p43*BW)*h + p44*k2N*eta
    ayprime(ix_df + 8) = k*(p45*dx + p46*dy + p47*dz) + B*(p48*dm + k*(p49*x + p50*y + p51*z)) + (p52*dB + p53*B2 + p54*k2)*m + kN*(p55*dh + p56*deta) + k*((p57*BN + p58*W)*h + (p59*BN + p60*W)*eta)


    end subroutine dfDeriv

    function a2X(this,ix_df,ix_etak,ay,k,k2,W,N,B2,W2,N2,kN,BW,BN,k2N,k2N2)
    class(TDarkFieldModel), intent(inout) :: this
    integer ix_df,ix_etak
    real(dl), target :: ay(:)
    real(dl) :: k,k2,W,N,B2,W2,N2,kN,BW,BN,k2N,k2N2
    real(dl) :: a2X
    real(dl) :: x,y,z,m,dx,dm
    real(dl) :: h,eta
    real(dl) :: p61,p62,p63,p64,p65,p66,p67,p68,p69,p71,p72,p73,p74

    x = ay(ix_df + 1)
    y = ay(ix_df + 2)
    z = ay(ix_df + 3)
    m = ay(ix_df + 4)
    dx = ay(ix_df + 5)
    !dy = ay(ix_df + 6)
    !dz = ay(ix_df + 7)
    dm = ay(ix_df + 8)

    h = ay(ix_df)
    eta = ay(ix_etak)/k

    p61 = this%p(61)
    p62 = this%p(62)
    p63 = this%p(63)
    p64 = this%p(64)
    p65 = this%p(65)
    p66 = this%p(66)
    p67 = this%p(67)
    p68 = this%p(68)
    p69 = this%p(69)
    !p70 = this%p(70)
    p71 = this%p(71)
    p72 = this%p(72)
    p73 = this%p(73)
    p74 = this%p(74)

    a2X = p61*W*dx + p62*kN*dm + ((p63*B2 + p64*k2)*N + p65*BW)*x + k2N*(p66*y + p67*z) + (p68*BN + p69*W)*k*m + ((p71*B2 + p72*k2)*N2 + p73*W2)*h + p74*k2N2*eta

    end function a2X

    function ia2M(this,ix_df,ix_etak,ay,ayprime,k,W,N,N2,BN,kN)
    class(TDarkFieldModel), intent(inout) :: this
    integer ix_df,ix_etak
    real(dl), target :: ay(:),ayprime(:)
    real(dl) :: k,W,N,N2,BN,kN
    real(dl) :: ia2M
    real(dl) :: x,y,m,dx,dy,dz
    real(dl) :: h,eta,dh
    real(dl) :: p75,p76,p77,p78,p79,p80,p81,p82,p83,p85,p86,p87,p88

    x = ay(ix_df + 1)
    y = ay(ix_df + 2)
    !z = ay(ix_df + 3)
    m = ay(ix_df + 4)
    dx = ay(ix_df + 5)
    dy = ay(ix_df + 6)
    dz = ay(ix_df + 7)
    !dm = ay(ix_df + 8)

    h = ay(ix_df)
    eta = ay(ix_etak)/k

    dh = ayprime(ix_df)

    p75 = this%p(75)
    p76 = this%p(76)
    p77 = this%p(77)
    p78 = this%p(78)
    p79 = this%p(79)
    p80 = this%p(80)
    p81 = this%p(81)
    p82 = this%p(82)
    p83 = this%p(83)
    !p84 = this%p(84)
    p85 = this%p(85)
    p86 = this%p(86)
    p87 = this%p(87)
    p88 = this%p(88)

    ia2M = N*(p75*dx + p76*dy + p77*dz) + (p78*BN + p79*W)*x + (p80*BN + p81*W)*y + p82*kN*m + N2*(p83*dh) + N*((p85*BN + p86*W)*h + (p87*BN + p88*W)*eta)

    end function ia2M

    end module DarkFieldModule
