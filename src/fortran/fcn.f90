subroutine fcn(t,qq,dqq)
  implicit double precision(a-h,o-z)
  integer n
  parameter (n=3)
  dimension dq(n), qq(1 + n, n), dqq(1 + n, n)
  b = 4.0
  sg = 16.0
  r = 45.92

  x = qq(1, 1)
  y = qq(1, 2)
  z = qq(1, 3)

  dx = sg*(y - x)
  dy = -x*z + r*x - y
  dz = x*y - b*z

  dqq(1, 1) = dx
  dqq(1, 2) = dy
  dqq(1, 3) = dz

  do j=1,n
     xp = qq(1 + j, 1)
     yp = qq(1 + j, 2)
     zp = qq(1 + j, 3)

     dxp = sg*(yp - xp)
     dyp = (r - z)*xp - yp - x*zp
     dzp = y*xp + x*yp - b*zp

     dqq(1 + j, 1) = dxp
     dqq(1 + j, 2) = dyp
     dqq(1 + j, 3) = dzp
  end do

end subroutine fcn
