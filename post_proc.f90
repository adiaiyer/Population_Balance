program post_proc

  implicit none
  
  integer, PARAMETER :: nlines
  real(kind=rprec),dimension(nlines) :: u,v,w,dissip
  open(unit=12,file="fort.12001000")
  do i=1:nlines
    READ(12,*) x(i) , y(i) , z ,u , v , w 
  enddo   
  dissip = u**3+v**3
