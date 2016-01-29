      real*8 function interp_project(PTX,PTY,PTZ,pos,XIMAX,ZIMAX,YIMAX, &
     &array,INTA,zsign,offx,offy,offz)
      implicit  none
      integer :: XIMAX,YIMAX,ZIMAX,I,J,K,N,M
      integer,dimension(0:4) :: PTX,PTY,PTZ
      integer,dimension(0:2) :: INTA
      real*8 :: offx,offy,offz,ret
      real*8,dimension(1:XIMAX+2,1:ZIMAX+2,1:YIMAX) :: array
      real*8,dimension(0:2) :: pos
      real*8,dimension(0:4) :: zsign
      real*8,dimension(0:4) :: xvalues,yvalues,array_interp_2
      real*8,dimension(0:24) :: array_interp_1
      real*8 :: interp
      N=0      
        do K=0,4,1    
          do J=0,4,1      
          M=0      
            do I=0,4,1        
            xvalues(M)=INTA(0)-2+I+offx
            yvalues(M)=zsign(K)*(array(PTX(I)+2,PTZ(K)+2,PTY(J)+1))
            M=M+1
            enddo
          array_interp_1(N)=interp(yvalues,xvalues,pos(0),5)
          N=N+1
          enddo     
        enddo
      N=0
      M=0
           
        do K=0,4,1
          do J=0,4,1   
          xvalues(J)=INTA(1)-2+J+offy
          yvalues(J)=array_interp_1(N)
          N=N+1
          enddo        
        array_interp_2(M)=interp(yvalues,xvalues,pos(1),5)
        M=M+1
        enddo
      I=0
      
        do K=0,4,1        
        xvalues(I)=INTA(2)-2+K+offz
        yvalues(I)=array_interp_2(I)
        I=I+1
        enddo
      interp_project=interp(yvalues,xvalues,pos(2),5)
      return
      endfunction interp_project

