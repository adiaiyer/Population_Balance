      PROGRAM READ_VEL_OUT

      IMPLICIT NONE

      INTEGER,PARAMETER :: nx = 168
      INTEGER,PARAMETER :: ny = 84
      INTEGER,PARAMETER :: nz = 121
      INTEGER,PARAMETER :: lh=nx/2+1,ld=2*lh

      CHARACTER(len=30) :: fname
      REAL(KIND=8),DIMENSION(ld,ny,nz) :: u
      REAL(KIND=8),DIMENSION(ld,ny,nz) :: v
      REAL(KIND=8),DIMENSION(ld,ny,nz) :: w
      REAL(KIND=8),DIMENSION(ld,ny,nz) :: theta
      REAL(KIND=8),DIMENSION(ld,ny,nz) :: RHSx
      REAL(KIND=8),DIMENSION(ld,ny,nz) :: RHSy
      REAL(KIND=8),DIMENSION(ld,ny,nz) :: RHSz
      REAL(KIND=8),DIMENSION(ld,ny,nz) :: RHS_T
      REAL(KIND=8),DIMENSION(ld,ny,1)  :: sgs_t3
      REAL(KIND=8),DIMENSION(nx,ny)    :: psi_m
      REAL(KIND=8),DIMENSION(ld,ny,nz) :: Cs_opt2
      REAL(KIND=8),DIMENSION(ld,ny,nz) :: F_LM
      REAL(KIND=8),DIMENSION(ld,ny,nz) :: F_MM
      REAL(KIND=8),DIMENSION(ld,ny,nz) :: F_QN
      REAL(KIND=8),DIMENSION(ld,ny,nz) :: F_NN

      fname='vel.out'
      open(10,file=fname,form='unformatted')
      READ (10) u(:,:,1:nz),v(:,:,1:nz),w(:,:,1:nz),   &
                RHSx(:,:,1:nz),RHSy(:,:,1:nz),RHSz(:,:,1:nz),   &
                Cs_opt2,F_LM,F_MM, &
                F_QN,F_NN

      PRINT*,v(:,1,120)

      END PROGRAM READ_VEL_OUT
