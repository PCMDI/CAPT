      subroutine filterG(f,w,n,m,kindex)
      real f(n,m), w(n,m)
      nn = n - 1
      mm = m - 1
      
      do 10 k = 1,kindex
         do 20 j = 1,m
            do 30 i = 2,nn
               w(i,j) = 0.25*(f(i-1,j) + f(i+1,j)) + 0.5*f(i,j)
   30       continue
            w(1,j) = f(1,j)
            w(n,j) = f(n,j)
   20    continue
         do 40 i = 1,n
            do 50 j = 2,mm
               f(i,j) = 0.25*(w(i,j+1) + w(i,j-1)) + 0.5*w(i,j)
   50       continue
            f(i,1) = w(i,1)
            f(i,m) = w(i,m)
   40    continue
   10 continue

      xnu = -0.5*float(kindex)
      xend = 0.5*xnu
      xmid = 1.0 - 2.0*xend
      do 60 j=1,m
         do 70 i = 2,nn
            w(i,j) = xend*(f(i+1,j) + f(i-1,j)) + xmid*f(i,j)
   70    continue
         w(1,j) = f(1,j)
         w(n,j) = f(n,j)
   60 continue
      do 80 i = 1,n
         do 90 j=2,mm
            f(i,j) = xend*(w(i, j+1) + w(i,j-1)) + xmid*w(i,j)
   90    continue
         f(i,1) = w(i,1)
         f(i,m) = w(i,m)
   80 continue
      return
      end
