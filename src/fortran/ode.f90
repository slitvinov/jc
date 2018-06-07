program ode
  implicit none
  integer n, j, k, l, m, nstep, io, irate
  parameter (n=3)
  real(kind=8) A, B, C, D, cum, gsc, dt, t, tme
  real(kind=8) y, znorm, norm, q, qq
  external fcn
  dimension y(1 + n, n), A(1 + n, n), B(1 + n, n), C(1 + n, n), D(1 + n, n)

  dimension q(n), qq(1 + n, n)
  dimension cum(n),znorm(n),gsc(n)
  write(*,111)
111 format(1x,'nstep, dt, irate, io :')
  read(*,*)  nstep,dt,irate,io

  ! initial conditions for nonlinear ODEs
  tme=0.0
  qq(1, :) =  1.0

  qq(2 : , :) = 0
  do j = 1, n; qq(1 + j, j) = 1
  end do

  cum = 0
  do 100 m=1,nstep
     do j=1,irate
        y=qq
        t=tme
        call fcn(t, y, A)
        y=qq+(dt*A)/2.0
        t=tme+dt/2.0
        call fcn(t, y, B)
        y=qq+(dt*B)/2.0
        t=tme+dt/2.0
        call fcn(t, y, C)
        y=qq+(dt*C)
        t=tme+dt
        call fcn(t, y, D)
        qq=qq+dt*(A+D+2.0*(B+C))/6.0
        tme=tme+dt
     end do

     do j=1,n
        q = qq(1 + j, :)
        do k=1,j - 1; gsc(k)= sum(q*qq(1 +  k, :))
        end do
        do l=1,j - 1; q=q - gsc(l)*qq(1 + l, :)
        end do
        norm = sqrt(sum(q*q))
        q = q/norm
        znorm(j) = norm
        qq(1 + j, :) = q
     end do

     !     update running vector magnitudes
     cum = cum + log(znorm)/log(2.0)
     if(mod(m,io).ne.0) goto 100
     write(*,334) tme,(cum(k)/tme,k=1,n)
334  format(1x,f12.6,2x,e12.6,2x,e12.6,2x,e12.6)
100  continue
   end program
