program main
  use mod_cmLLGS
  implicit none
  type(cmLLGS) cal
  
  call cal%init()
  call cal%readparam()

  call cal%calcmh()
  

  end program
