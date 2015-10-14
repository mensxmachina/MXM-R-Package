subroutine cor(n, x, y, output)

!dimension x(*), y(*)
!real, dimension(n):: x,y

integer n
double precision x(n), y(n)


double precision:: output
double precision:: sxx, sxy, syy, xsum, ysum, xysum, x2sum, y2sum, rn
      !eps = epsilon(4.0)
!
      !ier = 0
!      do 10 i = 1 , 7
!         output(i) = 0.0
! 10   continue
!
        output = 0.0;
        
      !if (n .le. 2) then
         !ier = 1
         !return
      !endif
!
!      Initialize accumulators
!
!      xsum --  sum of the x observations
!      ysum --  sum of the y observations
!      xysum--  sum of the products of the observations of x and y
!      sxx  --  sum of the squared deviations of x from the mean
!      syy  --  sum of the squared deviations of y from the mean
!      sxy  --  sum of the product of the deviations in x and in y
!
      sxx   = 0.0
      sxy   = 0.0
      syy   = 0.0
      xsum  = 0.0
      ysum  = 0.0
      xysum = 0.0
      x2sum = 0.0
      y2sum = 0.0
      rn    = float(n)
!
!      Accumulate observations, squares, and calculate the means
!
      do 20 i = 1,n
           xsum = xsum + x(i)
           ysum = ysum + y(i)
           xysum = xysum + (x(i) * y(i))
 20   continue

      do 30 i = 1,n
           x2sum = x2sum + x(i) * x(i)
           y2sum = y2sum + y(i) * y(i)
 30   continue
      xsum2 = xsum * xsum
      ysum2 = ysum * ysum
!
!     Set output   to the Pearson product-moment correlation coefficien
!
      cornum = xysum - ((xsum * ysum) / rn)
      corden = sqrt ((x2sum - (xsum2 / rn)) * (y2sum - (ysum2 / rn)))
!      if (corden .le. eps) then
!         ier = 3
!         return
!      endif
      output = cornum / corden
      return
end subroutine cor
