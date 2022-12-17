program simple
    implicit real*8(a-h,o-z)
    dimension a(2,2)
    a= reshape((/1.d0,2.d0,3.d0,4.d0/),shape(a))
    write(6,*) a
    write(6,*) a(1,1), a(1,2)
end program simple