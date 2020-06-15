      subroutine filter25(zi,z,li,lj,nfilt)
c zi is the input/output array, z is a work array
c 25 point filter designed by Rainer Bleck (his Masters work)
      real z(li, lj), zi(li, lj)
      li1 = li - 1
      lj1 = lj - 1
      li2 = li - 2
      lj2 = lj - 2

      do 1 n = 1,nfilt
c smooth left and right boundaries with 3pt 1 2 1 smoother
         do 10 j = 2,lj1
            z(1,j)  = 0.25*zi(1, j+1) + 0.25*zi(1, j-1) + 0.5*zi(1, j)
            z(li,j) = 0.25*zi(li,j+1) + 0.25*zi(li,j-1) + 0.5*zi(li,j)
   10    continue
c smooth top and bottom with 3pt
         do 20 i=2,li1
            z(i,1)  = 0.25*zi(i+1, 1) + 0.25*zi(i-1, 1) + 0.5*zi(i, 1)
            z(i,lj) = 0.25*zi(i+1,lj) + 0.25*zi(i-1,lj) + 0.5*zi(i,lj)
   20    continue
c left/right first interior with 9pt 1 2 1
         do 30 j =2,lj1
            z(2,j) = ((zi(3,j+1) + zi(1,j-1) + zi(3,j-1) + zi(1,j+1))
     $              + 2.0*(zi(3,j) + zi(2,j+1) + zi(2,j-1) + zi(1,j))
     $              + 4.0*zi(2,j))/16.0
            z(li1,j) =((zi(li,j+1)+zi(li2,j-1)+zi(li,j-1)+zi(li2,j+1))
     $              + 2.0*(zi(li,j)+zi(li1,j+1)+zi(li1,j-1)+zi(li2,j))
     $              + 4.0*zi(li1,j))/16.0
   30    continue
c top/bottom first interior pts with p pt 1 2 1
         do 40 i = 2,li1
            z(i,2)=((zi(i+1,3)+zi(i-1,1)+zi(i+1,1)+zi(i-1,3))
     $             +2.0*(zi(i+1,2)+zi(i,3)+zi(i,1)+zi(i-1,2))
     $             +4.0*zi(i,2))/16.0
            z(i,lj1)=((zi(i+1,lj)+zi(i-1,lj2)+zi(i+1,lj2)+zi(i-1,lj))
     $            +2.0*(zi(i+1,lj1)+zi(i,lj)+zi(i,lj2)+zi(i-1,lj1))
     $             +4.0*zi(i,lj1))/16.0
   40    continue
c corner points
         z(1,1) = (z(2,1) + z(1,2) + z(2,2))/3.0
         z(1,lj) = (z(1,lj1) + z(2,lj) + z(2,lj1))/3.0
         z(li,lj) = (z(li1,lj) + z(li, lj1) + z(li1,lj1))/3.0
         z(li,1) = (z(li,2) + z(li1, 1) + z(li1,2))/3.0
C now for the rest
         do 50 i=3,li2
            do 50 j = 3,lj2
               z(i,j)=0.279372*zi(i,j)
     $               +0.171943*(zi(i-1,j)+zi(i,j-1)+zi(i+1,j)+zi(i,j+1))
     $               -0.006918*(zi(i-2,j)+zi(i,j-2)+zi(i+2,j)+zi(i,j+2))
     $       +0.077458*(zi(i-1,j-1)+zi(i+1,j+1)+zi(i+1,j-1)+zi(i-1,j+1))
     $       -0.024693*(zi(i-1,j-2)+zi(i+1,j-2)+zi(i-2,j-1)+zi(i+2,j-1)+
     $                  zi(i-2,j+1)+zi(i+2,j+1)+zi(i-1,j+2)+zi(i+1,j+2))
     $       -0.012940*(zi(i-2,j-2)+zi(i+2,j-2)+zi(i-2,j+2)+zi(i+2,j+2))
   50    continue
         do 60 i = 1,li
            do 60 j = 1,lj
               zi(i,j) = z(i,j)
   60   continue
    1 continue
      return
      end
