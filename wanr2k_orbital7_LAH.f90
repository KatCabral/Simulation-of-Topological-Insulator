      Program Wannier_band_structure      
      Implicit None
!--------to be midified by the usere
      character(len=80):: prefix="BiTeI"
      integer,parameter::nkpath=3,np=1000
      !real*8,parameter::ef= 4.18903772
      real*8 ef, eg
      real*8 max_band12
      real*8 min_band13
!------------------------------------------------------
      integer*4 ik,ikmax
      real*8 kz
      real*8 alpha
      character(len=30)::klabel(nkpath)
      character(len=80) hamil_file_triv, hamil_file_top,nnkp,line
      integer*4,parameter::nk=(nkpath-1)*np+1
      integer*4 i,j,k,nr,nr_top,i1,i2,i1_top,i2_top,nb,nb_top,lwork,info
      real*8,parameter::third=1d0/3d0!,kz=0d0
      real*8 phase,phase_top,pi2,jk,a_triv,b_triv, a_top, b_top
      real*8 klist(3,1:nk),xk(nk),kpath(3,np),bvec(3,3),ktemp1(3),ktemp2(3),xkl(nkpath)
      real*8,allocatable:: HK2(:,:),rvec(:,:),rvec_top(:,:),ene(:,:),rwork(:),&
              prob_te(:,:), prob_bi(:,:), prob_i(:,:)
      integer*4,allocatable:: ndeg(:), ndeg_top(:)
      integer*4 bi(3),te(3),id(3),orb(3)
      complex*16,allocatable:: Hk_triv(:,:), Hk_top(:,:), Hamr_triv(:,:,:),Hamr_top(:,:,:), &
              work(:), HK_alpha(:,:)
      complex*16 temp1,temp2
!-----------------------------------------------------
      interface write_plt
          subroutine write_plt_interface(nkp, xkl, kl, ef, prob_te, prob_bi, prob_i)
              implicit none
              integer nkp, i 
              real*8 xkl(nkp),ef, prob_te(:,:), prob_bi(:,:),prob_i(:,:)
              character(len=30)kl(nkp)
          end subroutine write_plt_interface      
      end interface 
!------------------------------------------------------
      data bi / 255, 000,  000 /
      data te / 000, 000,  255 /
      data id / 000, 255,  000 /
!------------------------------------------------------
      write(hamil_file_triv,'(a,a)')trim(adjustl(prefix)),"_hr_triv.dat"
      write(hamil_file_top,'(a,a)')trim(adjustl(prefix)),"_hr_top.dat"
      write(nnkp,'(a,a)')      trim(adjustl(prefix)),".nnkp"

      pi2=4.0d0*atan(1.0d0)*2.0d0
!---------------  reciprocal vectors
      open(98,file=trim(adjustl(nnkp)),err=333)
111   read(98,'(a)')line
      if(trim(adjustl(line)).ne."begin recip_lattice") goto 111
      
      read(98,*)bvec
!---------------kpath
      data kpath(:,1) /     0.5d0,      0.0d0,    0.5d0/  !L
      data kpath(:,2) /     0.0d0,      0.0d0,    0.5d0/  !A
      data kpath(:,3) /     third,      third,    0.5d0/  !H
!      open(77,file='tmp')
 !     read(77,*)ik,ikmax
 !     kz=float(ik)*0.5d0/float(ikmax)
 !     kpath(3,:)=kz

      data klabel     /'L','A','H'/

      ktemp1(:)=(kpath(1,1)-kpath(1,2))*bvec(:,1)+(kpath(2,1)-kpath(2,2))*bvec(:,2)+(kpath(3,1)-kpath(3,2))*bvec(:,3)

!      xk(1)= 0d0 !-sqrt(dot_product(ktemp1,ktemp1))
      xk(1)= -sqrt(dot_product(ktemp1,ktemp1))
      xkl(1)=xk(1)
      

      k=0
      ktemp1=0d0
      do i=1,nkpath-1
       do j=1,np
        k=k+1
        jk=dfloat(j-1)/dfloat(np)
        klist(:,k)=kpath(:,i)+jk*(kpath(:,i+1)-kpath(:,i))
        ktemp2=klist(1,k)*bvec(:,1)+klist(2,k)*bvec(:,2)+klist(3,k)*bvec(:,3)
        if(k.gt.1) xk(k)=xk(k-1)+sqrt(dot_product(ktemp2-ktemp1,ktemp2-ktemp1))
        if(j.eq.1) xkl(i)=xk(k)
        ktemp1=ktemp2
       enddo
      enddo
      klist(:,nk)=kpath(:,nkpath)
      ktemp2=klist(1,nk)*bvec(:,1)+klist(2,nk)*bvec(:,2)+klist(3,nk)*bvec(:,3)
      xk(nk)=xk(nk-1)+sqrt(dot_product(ktemp2-ktemp1,ktemp2-ktemp1))
      xkl(nkpath)=xk(nk)
!      write(*,*)klist
      klist=klist*pi2

!------read H(R) trivial
      open(99,file=trim(adjustl(hamil_file_triv)),err=444)
      read(99,*)
      read(99,*)nb,nr
      allocate(rvec(3,nr),Hk_triv(nb,nb),Hamr_triv(nb,nb,nr),ndeg(nr))
      read(99,*)ndeg
      do k=1,nr
         do i=1,nb
            do j=1,nb
               read(99,*)rvec(1,k),rvec(2,k),rvec(3,k),i1,i2,a_triv,b_triv
               hamr_triv(i1,i2,k)=dcmplx(a_triv,b_triv)
            enddo
         enddo
      enddo

      

      lwork=max(1,2*nb-1)
      allocate(work(max(1,lwork)),rwork(max(1,3*nb-2)))

      

 
!------read H(R) topological
      open(103,file=trim(adjustl(hamil_file_top)),err=666)
      read(103,*)
      read(103,*)nb_top,nr_top
      allocate(rvec_top(3,nr),Hk_top(nb,nb),Hamr_top(nb,nb,nr),ndeg_top(nr_top))
      read(103,*)ndeg_top
      do k=1,nr_top
         do i=1,nb_top
            do j=1,nb_top
               read(103,*)rvec_top(1,k),rvec_top(2,k),rvec_top(3,k),i1_top,i2_top,a_top,b_top
               hamr_top(i1_top,i2_top,k)=dcmplx(a_top,b_top)
            enddo
         enddo
      enddo
 
      

 
!---- Fourrier transform H(R) H(k)
      allocate(HK2(nb,nb),HK_alpha(nb,nb),ene(nb,nk),prob_te(nb,nk),prob_bi(nb,nk),prob_i(nb,nk))
      open(100,file='bandtop.dat')
      open(111,file='bandrgbtop.dat')
      ene=0d0
      prob_te = 0d0
      prob_bi = 0d0
      prob_i = 0d0
      alpha =0.7754d0
      do k=1,nk
         HK_triv=(0d0,0d0)
         HK_top=(0d0,0d0)
         HK_alpha=(0d0,0d0)
         
         do j=1,nr
            phase=0.0d0
            phase_top=0.0d0
            do i=1,3
               phase=phase+klist(i,k)*rvec(i,j)
               phase_top=phase_top+klist(i,k)*rvec_top(i,j) 
            enddo
      
            do i1=1,nb
               do i2=1,nb
                  Hk_triv(i1,i2)=Hk_triv(i1,i2)+Hamr_triv(i1,i2,j)* &
                  dcmplx(cos(phase),sin(phase))/float(ndeg(j))
                 
                  Hk_top(i1,i2)=Hk_top(i1,i2)+Hamr_top(i1,i2,j)* &
                  dcmplx(cos(phase_top),sin(phase_top))/float(ndeg_top(j))
                  
               enddo
            enddo
       
         enddo
       
         
         HK_alpha= alpha*Hk_top + (1-alpha)*Hk_triv
         
         call zheev('V','U',nb,HK_alpha,nb,ene(:,k),work,lwork,rwork,info)
                  
        
           
         HK2=dconjg(HK_alpha)*HK_alpha
                  
         do i1=1,nb
            prob_te(i1,k)=sum(HK2(1:3,i1))+sum(HK2(10:12,i1))
            prob_bi(i1,k)=sum(HK2(4:6,i1))+sum(HK2(13:15,i1))
            prob_i(i1,k) =sum(HK2(7:9,i1))+sum(HK2(16:18,i1))
           
              
         enddo
      enddo
      

      max_band12 = maxval(ene(12,:),dim=1)
      print *, "max_band12", max_band12
      min_band13 = minval(ene(13,:),dim=1)
      print *, "min_band13", min_band13
      ef = (max_band12 + min_band13) / 2
      print *, "ef", ef
      eg = min_band13 - max_band12
      print *, "eg", eg
      do i=1,nb
         do k=1,nk

           orb = prob_bi(i,k)*bi+prob_te(i,k)*te+prob_i(i,k)*id
           write(100,'(5(x,f12.6))') xk(k),ene(i,k), prob_te(i,k), prob_bi(i,k), prob_i(i,k)
           write(111,'(2(x,f12.6),3(x,i3))') xk(k),ene(i,k),orb
         enddo
           write(100,*)
           write(100,*)
           write(111,*)
           write(111,*)
      enddo
      call write_plt_interface(nkpath,xkl,klabel,ef,prob_te,prob_bi,prob_i)
      stop
333   write(*,'(3a)')'ERROR: input file "',trim(adjustl(nnkp)),' not found'
      stop
444   write(*,'(3a)')'ERROR: input file "',trim(adjustl(hamil_file_triv)),' not found'
      stop
666   write(*,'(3a)')'ERROR: input file "',trim(adjustl(hamil_file_top)),' not found'
      stop
      end

      subroutine write_plt_interface(nkp,xkl,kl,ef,prob_te,prob_bi,prob_i)
      implicit none
      integer nkp, i
      real*8 xkl(nkp),ef, prob_te(:,:), prob_bi(:,:),prob_i(:,:)
      character(len=30)kl(nkp)
      close(99)
      open(99,file='bandrgb.plt')
      write(99,'(a,f12.8)')'ef=',ef
     
      write(99,'(a)') 'set xtics ( \'
      do i=1,nkp
         if(trim(adjustl(kl(i))).eq.'g'.or.trim(adjustl(kl(i))).eq.'G')kl(i)="{/Symbol \107}"
         if(i.ne.nkp) write(99,'(3a,f12.6,a)')'"',trim(adjustl(kl(i))),'"',xkl(i),", \"
         if(i.eq.nkp) write(99,'(3a,f12.6,a)')'"',trim(adjustl(kl(i))),'"',xkl(i)," )"
      enddo
      write(99,'(a)') 'set ytics'

      
      write(99,'(a,f12.6,a,f12.6,a)') 'set xrange [',xkl(1)*0.2,':',xkl(nkp)*0.2,']'
      write(99,'(a)') &
           'set terminal pdfcairo enhanced font "DejaVu"  transparent fontscale 1 size 5.00in, 7.50in'
      write(99,'(a,f4.2,a)')'set output "bandrgb.pdf"'
      write(99,'(11(a,/),a)') &
            'set style line 10 lt 1 lc rgb "black" lw 2',&
            'set encoding iso_8859_1',&
            'set size ratio 0 1.0,1.0',&
            'set ylabel "E-Ef (eV)"',&
            'set yrange [ -0.3 : 0.3 ]',&
            'unset key',&
            'set ytics 1.0 scale 1 nomirror out',&
            'set mytics 2',&
            'set parametric',&
            'set trange [-10:10]',&
            'rgb(r,g,b) = int(r)*65536 + int(g)*256 + int(b)'
      write(99,'(a)') &
            'plot "bandrgb.dat" using 1:($2-ef):(rgb($3,$4,$5)) with l lw 10 lc rgb variable,\'
      do i=2,nkp-1
            write(99,'(f12.6,a)')xkl(i),',t with l ls 10,\'
      enddo
      write(99,'(a)') 't,0 with l ls 10'
 
      end subroutine write_plt_interface

