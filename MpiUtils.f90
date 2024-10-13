    module MpiUtils
    implicit none
    ! use mpi leads to .mod compiler incompatibility errors unless you are very careful
    ! so stick to old method and add manual interface for gcc10+ compatibility
    include "mpif.h"
    interface
    subroutine MPI_BCAST(BUFFER, COUNT, DATATYPE, ROOT, COMM, IERROR)
     !GCC$ ATTRIBUTES NO_ARG_CHECK :: BUFFER
    Type(*) BUFFER
    INTEGER COUNT, DATATYPE, ROOT, COMM, IERROR
    end subroutine
    end interface

    integer, parameter :: TTimer_dp = Kind(1.d0)

    Type TTimer
        real(TTimer_dp) start_time
    contains
    procedure :: Start => TTimer_Start
    procedure :: Time => TTimer_Time
    procedure :: WriteTime => TTimer_WriteTime
    end type TTimer

    contains

    function GetMpiRank()
    integer GetMpiRank
    integer ierror
    call mpi_comm_rank(mpi_comm_world,GetMPIrank,ierror)
    if (ierror/=MPI_SUCCESS) call MpiStop('MPI fail')

    end function GetMpiRank

    function IsMainMPI()
    logical IsMainMPI

    IsMainMPI =  GetMpiRank() == 0

    end function IsMainMPI

    subroutine MpiStop(Msg)
    character(LEN=*), intent(in), optional :: Msg
    integer i
    integer ierror, MpiRank

    if (present(Msg)) write(*,*) trim(Msg)

    call mpi_comm_rank(mpi_comm_world,MPIrank,ierror)
    write (*,*) 'MpiStop: ', MpiRank
    call MPI_ABORT(MPI_COMM_WORLD,i, ierror)
    i=1     !put breakpoint on this line to debug
    stop

    end subroutine MpiStop

    subroutine MpiStat(MpiID, MpiSize)
    implicit none
    integer MpiID,MpiSize
    integer ierror
    call mpi_comm_rank(mpi_comm_world,MpiID,ierror)
    if (ierror/=MPI_SUCCESS) call MpiStop('MpiStat: MPI rank')
    call mpi_comm_size(mpi_comm_world,MpiSize,ierror)
    end subroutine MpiStat

    subroutine MpiQuietWait
    !Set MPI thread to sleep, e.g. so can run openmp on cpu instead
    integer ierr, STATUS(MPI_STATUS_SIZE)
    logical flag
    integer i, MpiId, MpiSize

    call MpiStat(MpiID, MpiSize)
    if (MpiID/=0) then
        do
            call MPI_IPROBE(0,0,MPI_COMM_WORLD,flag, MPI_STATUS_IGNORE,ierr)
            if (flag) then
            call MPI_RECV(i,1,MPI_INTEGER, 0,0,MPI_COMM_WORLD,status,ierr)
            exit
        end if
        call sleep(1)
    end do
    end if
    end subroutine

    subroutine MpiWakeQuietWait
    integer j, MpiId, MpiSize, ierr,r

    call MpiStat(MpiID, MpiSize)
    if (MpiID==0) then
        do j=1, MpiSize-1
            call MPI_ISSEND(MpiId,1,MPI_INTEGER, j,0,MPI_COMM_WORLD,r,ierr)
        end do
    end if
    end subroutine MpiWakeQuietWait

    subroutine MpiShareString(S, from)
    character(LEN=:), allocatable :: S
    integer from
    integer inlen, rank, ierror

    rank = GetMpiRank()

    if (rank==from) inlen=len(S)

    CALL MPI_Bcast(inlen, 1, MPI_INTEGER, from, MPI_COMM_WORLD, ierror)
    if (ierror/=MPI_SUCCESS) call MpiStop('MpiShareString: fail')

    if (rank /= from ) allocate(character(inlen)::S)
    CALL MPI_Bcast(S, LEN(S), MPI_CHARACTER, from, MPI_COMM_WORLD, ierror)
    end subroutine MpiShareString


    function TimerTime()
    real(TTimer_dp) time
    real(TTimer_dp) :: TimerTime
    !$ real(TTimer_dp), external :: omp_get_wtime
    TimerTime = MPI_WTime()
    end function TimerTime

    subroutine TTimer_Start(this, time)
    class(TTimer) :: this
    real(TTimer_dp), intent(out), optional :: time
    this%start_time = TimerTime()
    if (present(time)) time = this%start_time
    end subroutine TTimer_Start

    real(TTimer_dp) function TTimer_Time(this)
    class(TTimer) :: this
    TTimer_Time =  TimerTime() - this%start_time
    end function TTimer_Time

    subroutine TTimer_WriteTime(this,Msg, start)
    class(TTimer) :: this
    character(LEN=*), intent(in), optional :: Msg
    real(TTimer_dp), optional :: start
    real(TTimer_dp) T, DeltaT
    character(LEN=:), allocatable :: tmp

    if (present(start)) then
        T=start
    else
        T=this%start_time
    end if

    DeltaT = TimerTime() - T
    if (present(Msg)) then
        tmp = trim(Msg)//': '
        if (DeltaT > 0.00002 .and. DeltaT < 1000 .and. len_trim(tmp)<24) then
            write (*,'(a25,f10.5)') tmp, DeltaT
        else
            write (*,*) trim(Msg)//': ', DeltaT
        end if
    end if
    if (.not. present(start)) this%start_time = TimerTime()

    end subroutine TTimer_WriteTime


    end module MpiUtils
