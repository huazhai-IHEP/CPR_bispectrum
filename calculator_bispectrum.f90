    ! command:ifort
    ! f2py -c -m calculator_bispectrum calculator_bispectrum.f90 calculator_bispectrum_temp.f90 &
    ! wignerSymbols-fortran.f --opt='-O3' --fcompiler=intelem --f90flags="-qopenmp -D__OPENMP" --f77flags="-qopenmp -D__OPENMP" -liomp5
   
    !add cpt transformation inherited from camb/lensing.f90 by me
  
    !gfortran   
    !f2py -c -m calculator_bispectrum calculator_bispectrum.f90 --opt='-O3 -ffast-math'  ---ffree-line-length-none --fcompiler=gnu95 --f90flags="-fopenmp -D__OPENMP" --f77flags="-fopenmp -D__OPENMP" -lgomp
 
    ! if contains two files
	!f2py -c -m calculator_bispectrum calculator_bispectrum.f90 calculator_bispectrum_temp.f90 &
	!wignerSymbols-fortran.f --opt='-O3 -ffast-math' --fcompiler=gnu95 --f90flags="-fopenmp -D__OPENMP" --f77flags="-fopenmp -D__OPENMP" -lgomp
    
    !---------------------------------------------------------------------------------------------------------------------------------------
	SUBROUTINE calculator_bispectrum(t_or_p,  eod, CPT_ang_deg, Caa_in, CEa, CTa, CEE, CTE, bispectrum,  &
		lmax_alpha, lmax_ea, lstart, lmax_bis, lmax_teb, lmax_computed_cl )
    
    implicit none
       
    integer, parameter :: mcp=KIND(1.d0)
    real(mcp), parameter :: pi = 3.14159265358979323846264338328d0    
	INTEGER, intent(in) ::  lmax_alpha, lmax_teb, lstart, lmax_bis,lmax_ea,lmax_computed_cl
	! lmax_computed_cl is for temperature integration
	! a for alpha related power spectra, teb for normal power spectra, bis for output bispectrum_polar in one dimenssion. ususlly cut lmax_ea=500, \
	! cut lmax_teb=2500, lmax_teb >= lmax_ea+ lmax_bis
    real(mcp), intent(in) :: CEE(lmax_teb), CTE(lmax_teb), CEa(lmax_ea), Caa_in(lmax_ea), CTa(lmax_ea)
	! input alpha related power spectra are all scale invariant form.
	! note there is no lensing contribution on CBB,
        character(4), intent(in):: eod! even or odd
	real(mcp), intent(in) :: CPT_ang_deg
        character(2), intent(in)::t_or_p
	real(mcp), intent(out) :: bispectrum(9, lmax_bis)
        real(mcp) Caa(lmax_ea)
	integer :: i, j, j1, i1,i2,i3, j2, j3, IER
	real(mcp) CPT_ang, corr0
	! degree unit	 
	real(mcp) cs2, sn2, cs2sq, sn2sq, cs4, sn4
        real(mcp) :: bispectrum_temp(5, lmax_bis),  bispectrum_polar2(4, lmax_bis), bispectrum_polar3(4, lmax_bis)
	real(mcp) ::  rec_polar2(4), order1_polar(4)
	! two parts, one in two loops, another in three loops
        real(mcp) bispectrum_polar_order1(4, lmax_bis)
        ! store the first order of polarization
	real(mcp) Caa_corr0(lmax_ea)
	real(mcp) geom(lmax_bis)
	!geometrical factors, h(l1, l2,l3) and wigner 3j symbols I^{0-22}_{l1, l2, l3}, I^{0-22}_{l2, l3, l1}, I^{0-22}_{l3, l2, l1}
    ! bispectrum_polar order,4:EEE, 3:EEB, 2:EBB, 1:BBB	
	integer :: tp1, tp2
	! summation types indices
    integer :: l, lb, lb1, lb2, lb3, l1, l2, l3
	character(3) fields_list(4),  perm_list(6)
	character(3) :: perm_inds 
        character(4) :: num_file 
        character(55) :: file_name
	integer :: lblist(3)
	integer :: perm_lb1,perm_lb2, perm_lb3
	! lb1, lb2, lb3 for bispectrum_polar indices, l1, l2, l3 for summing indices
	! first index for fields, second for permutations
	real(mcp) fac
	integer:: botl, topl, nl
	real(mcp)  threej_m22(3, lmax_ea, 0:lmax_bis+lmax_ea+5), threej_20(3,lmax_ea, 0:lmax_bis+lmax_ea+5)
	real(mcp) threej_00(3,lmax_ea, 0:lmax_bis+lmax_ea+5)
	! precalculated 3j arrays to store, 
	!real(mcp) wig3j_000(0:2*lmax_bis+1)
	! for even geometrical factor
	real(mcp) wig3j_112(0:2*lmax_bis+5),wig3j_121(0:2*lmax_bis+5),wig3j_211(0:2*lmax_bis+5)
	! for odd geometrical factor
	! geometrical factor h(l1, l2, l3) and coeff I^{0, -2, 2}_{l1, l2, l3}, I^{0, -2, 2}_{l2, l1, l3}, I^{0, -2, 2}_{l3, l2, l2}
	real(mcp)  ypol(4, 6), sym_ypol(4, 6), perm_ypol(4, 6), tp_sum(4), tp1_sum(4)
	! store triple loop results, 4 means BBB, EBB, EEB, EEE. 6 for six perms, 2 for two perms
	character(10):: current_time
	

	!if (lmax_bis.gt.50)then
	!	lstart = lmax_bis-50
	!else
	!	lstart=3
	!end if 
	!lstart=3
        ! default from 3, can change free values

        print *, 'tp type; eod', t_or_p, eod
	print *, 'angle, lmax_ea, lstart, lmax_bis, lmax_teb, lmax_computed_cl', CPT_ang_deg,  lmax_ea, lstart, lmax_bis, lmax_teb, lmax_computed_cl

        if(lmax_teb.lt.lmax_ea)then
          print*,'CEE multipole index must be the lagest one!'
          stop
	end if 
        if(lmax_alpha .gt. lmax_ea)then
                print *, 'alphaalpha cut should be smaller than Ealpha cut'
                stop
        end if 
	
        ! make cut for C_ell^alphaalpha, set zero larger than lmax_alpha
        Caa = Caa_in
        Caa(lmax_alpha+1:) = 0

        CPT_ang = CPT_ang_deg*pi/180._mcp
        ! make small approximation when small than 0.5 deg to avoid trunck
        if(CPT_ang .lt. 0.0087)then
                sn2 = 2*CPT_ang-4*CPT_ang**3/3._mcp+4*CPT_ang**5/15._mcp
                cs2 = 1-2*CPT_ang**2+2*CPT_ang**4/3._mcp
        else
                sn2 = sin(2*CPT_ang)
                cs2 = cos(2*CPT_ang)
        end if 
        cs2sq = cs2**2
        sn2sq = sn2**2
        cs4 = cs2sq-sn2sq
        sn4 = 2*sn2*cs2
        ! variabce of caa
	! lcut =lmax_alpha
	Caa_corr0=0
        do l=2, lmax_alpha
          Caa_corr0(l) = 0.25_mcp*( 2*l + 1.)/pi*Caa(l)
        end do    
        corr0 =  sum(Caa_corr0(2:lmax_alpha))
        print *, 'corr0', corr0
        !fields_list = (/'BBB', 'EEB', 'EBB', 'EEE'/)

	call date_and_time(time=current_time)
	print *, 'begin calculation! at ', current_time
        
        bispectrum_temp=0
        if(t_or_p=='T'.or.t_or_p=='TP')then 
	! first calculate the temperature bispectrum
	        call  temp_bispectrum(eod, CPT_ang, CEE, CTE, Caa, CEa, CTa, bispectrum_temp,lmax_alpha,  lmax_ea, lstart, lmax_bis, lmax_teb, lmax_computed_cl)
        end if
	
        bispectrum=0
        do i=1, 5
             bispectrum(i, :) = bispectrum_temp(i, :)
        end do 
        geom = 0
	
        if(t_or_p=='P'.or.t_or_p=='TP')then 
                print *, 'calculate polarization bispectrum'
                do lb2 = lstart, lmax_bis-1
                ! specify lb1+1 = lb2= lb3-1, so lb2 is even number to satisfy lb1+lb2+lb3=even selection rule	
                ! usually select lmax_bis = 2000,

                ! to avoid the wigner 3j error, lb2 should start fronm 3, I^{0-22}_{lb1lb2lb3} loops.
                    !lb2 = l
                        lb1= 4
                        if(eod=='even')then
                                lb3 = lb2+4
                        else
                                lb3 = lb2+3
                        end if 
                        !precalaulated 3j arrays
                        threej_m22=0
                        threej_20=0
                        threej_00 = 0
                        !call cpu_time(time1)
                        
                        ! precalculate the wigner 3j symbols
                        call threej_store(lb1, lb2, lb3, threej_m22, threej_20,threej_00,  lmax_ea, lmax_bis)
                         
                        wig3j_112=0
                        wig3j_121=0
                        wig3j_211=0

                        topl =  lb1+lb3
                        !geometrical factor
                        if (eod=='even')then
                                botl = max(0, abs(lb1-lb3))
                                nl = topl- botl+1
                                ! for 3j(lb3,lb1,lb2;0,0,0) 
                                geom(lb2) = threej_00(3, lb1, lb2)
                        else
                                ! below 3 for odd geometrical factor
                                botl = max(1, abs(lb1-lb3))
                              nl = topl- botl+1	 
                                !3j(lb1,lb2,lb3; 1,1,-2)
                                call DRC3JJ(dble(lb3), dble(lb1), dble(-2), dble(1), dble(botl), dble(topl), wig3j_112(botl:), &
                                                nl, IER)

                                botl = max(2, abs(lb1-lb3))
                                nl = topl- botl+1	 		
                                !3j(lb1,lb2,lb3; 1,-2, 1)
                                call DRC3JJ(dble(lb3), dble(lb1), dble(1), dble(1), dble(botl), dble(topl), wig3j_121(botl:), &
                                                nl, IER)

                                botl = max(1, abs(lb1-lb3))
                                nl = topl- botl+1	 		
                                !3j(lb1,lb2,lb3;-2, 1, 1)
                                call DRC3JJ(dble(lb3), dble(lb1), dble(1), dble(-2), dble(botl), dble(topl), wig3j_211(botl:), &
                                                nl, IER)		
                                
                                !print*, 'threej_123(1:)',threej_123(1:lb2+lb3)	
                                geom(lb2) = (wig3j_112(lb2)+wig3j_121(lb2)+wig3j_211(lb2))/3._mcp
                        end if 


                        !print *, 'geom', geom(lb2)
                                              
                        rec_polar2=0
                        order1_polar= 0
                        ! the double loop results
                        call double_loop_polar( CPT_ang, lb1, lb2, lb3,  threej_m22, threej_20, threej_00, CEa, Caa, &
                         CEE, rec_polar2, order1_polar, lmax_alpha,  lmax_ea, lmax_bis, lmax_teb)

                        ! add geometrical factors
                        do i=1, 4
                                fac = 4*exp(-6*corr0)/geom(lb2)
                                bispectrum_polar_order1(i, lb2) = order1_polar(i)
                                !only first order
                                !bispectrum_polar2(i, lb2) = (order1_polar(i))*fac
                                bispectrum_polar2(i, lb2) = (order1_polar(i) + rec_polar2(i))*fac
                        end do 

                        print *, ' double ebb,', lb2, bispectrum_polar2(2, lb2)

                        !--------------------------------------------------------------------------------------------------------------
                        ! begin three loop part
                        perm_list = (/'123','132', '231', '213', '321', '312'/)
                        lblist=(/lb1, lb2, lb3/)
                       
                         
                        ypol = 0
                        sym_ypol = 0	
                        ! first dimension 6 perms, second 4 types 	
                        do j=1, 6
                                perm_inds = trim(perm_list(j))
                                read(perm_inds(1:1),'(I2)') i1
                                read(perm_inds(2:2),'(I2)') i2
                                read(perm_inds(3:3),'(I2)') i3
                        
                                perm_lb1 = lblist(i1)
                                perm_lb2 = lblist(i2)
                                perm_lb3 = lblist(i3)

                                                       
                                !if(perm_lb1==0)then
                                !print*,'perm_inds, lb1,lb2,lb3 ', perm_inds, perm_lb1,perm_lb2,perm_lb3
                                !stop
                                !end if
                
                                call triple_loop_polar(CPT_ang, tp2, perm_lb1, perm_lb2, perm_lb3, perm_inds, threej_m22, threej_20, CEa, Caa, CEE, &
                                                                ypol(:, j), lmax_alpha, lmax_ea, lmax_bis, lmax_teb)
                                                                                                
                                !for first type term, EEB has symmetry (l1,l2), corresponding to j=(1,4), others has symmetry (l2,l3) corresponding to j=(1,2)
                                if((j==1).or.(j==2).or.(j==4))then
                                        call triple_loop_polar(CPT_ang, tp1, perm_lb1, perm_lb2, perm_lb3, perm_inds, threej_m22, threej_20, CEa, Caa, CEE, &
                                                                        sym_ypol(:, j),lmax_alpha,  lmax_ea, lmax_bis, lmax_teb)					
                                end if
                                
                        end do
                        !print *, 'ypol(ebb-eee)', ypol(2, :)-ypol(4, :)

                        
                        fac = 8*exp(-6*corr0)/geom(lb2) 
                        if (eod=='even')then
                                
                                do i =1, 4
                                        ! we only need two perms for type 1
                                        if(i/=3)then
                                                tp_sum(i) = sum(ypol(i, :)) + sum(sym_ypol(i, 1:2))
                                        else
                                                ! for eeb
                                                tp_sum(3) = sum(ypol(3, :)) + sum(sym_ypol(3, [1,4]))
                                        end if 

                                        bispectrum_polar3(i, lb2) =  tp_sum(i)*fac*(1- 2*mod((lb1+lb2+lb3)/2, 2))
                                end do 
                        else
                                !odd parity change sign of geometrical factor when make single exchange perm, but invariant when double exchange.
                                ! for five perms, single exchange change sign, 123, 231, 312 positive, negative:132, 213, 321
                                ! or here vonventions, positive:1, 3, 6, negative:2, 4, 5

                                ! type 1, two perms, always change sign
                                ! bbb, eee, ebb, taken 23 perms
                                tp1_sum([1,2, 4]) =  sym_ypol([1,2, 4], 1) - sym_ypol([1,2, 4], 2)
                                !eeb take 12 perm 
                                tp1_sum(3) =  sym_ypol(3, 1) - sym_ypol(3, 4)


                                perm_ypol = ypol
                                perm_ypol(:, [2,4,5]) = -ypol(:, [2,4,5])
                                bispectrum_polar3(1, lb2) = ( sum(perm_ypol(1, :))+ tp1_sum(1) )*fac*(1- 2*mod((lb1+lb2+lb3-1)/2, 2))
                                bispectrum_polar3(4, lb2) = (sum(perm_ypol(4,:)) + tp1_sum(4))*fac*(1- 2*mod((lb1+lb2+lb3-1)/2, 2))
                                
                                !elemental one is 1<->2=213, 1<->3=321, the 2<->3 part change signs, 1,4, 5 invariant,  2,3,6 change sign
                                perm_ypol = ypol
                                perm_ypol(2, [2,3,6]) = -ypol(2, [2,3,6])
                                bispectrum_polar3(2, lb2) =  (sum(perm_ypol(2, :)) +tp1_sum(2) )*fac*(1- 2*mod((lb1+lb2+lb3-1)/2, 2))
                                
                                !elemental one is 2<->3=132, 1<->3=321, the 1<->2 part change signs . 1, 2, 5 invariant, 3, 4,6 change sign
                                perm_ypol = ypol
                                perm_ypol(3, [3,4,6]) = -ypol(3, [3,4,6])
                                bispectrum_polar3(3, lb2) =  ( sum(perm_ypol(3, :)) + tp1_sum(3) )*fac*(1- 2*mod((lb1+lb2+lb3-1)/2, 2))

                        end if 

                        !print *, 'triple  ebb-eee, ', lb2,   ypol(1, :), sym_ypol(1, :)

                        print *, 'triple ebb', bispectrum_polar3(2, lb2)

                        
                        do i= 6, 9
                                bispectrum(i, lb2) = bispectrum_polar2(i-5, lb2) + bispectrum_polar3(i-5, lb2)
                        end do 

                        !print *, ' double ebb,', lb2, bispectrum(7, lb2)
                        
                end do 
        end if 

        write(num_file, '(i4)')lstart
	write(file_name, '(i4)')lmax_bis
        file_name = './cpt_bispectrum_'//trim(t_or_p)//'_'//trim(eod)//'_'//trim(adjustl(num_file))//'-'//trim(adjustl(file_name))//'.dat'
        print *, 'file: ', file_name
        open(11, file=file_name, status='replace')
	write(11,"(10A5)") "# L", "TTB", "TTE", "TBB", "TEB","TEE", "BBB", "EBB", "EEB", "EEE"
	do l=1, lmax_bis
		write(11,"(I4, 9E20.10)") l,  bispectrum(1:9, l)                                            
	end do
	close(11)

        ! store first order of sepsilon result
        if(t_or_p=='P'.or.t_or_p=='TP')then
                write(num_file, '(i4)')lstart
                write(file_name, '(i4)')lmax_bis
                file_name = './cpt_bispectrum_order1'//trim(t_or_p)//'_'//trim(eod)//'_'//trim(adjustl(num_file))//'-'//trim(adjustl(file_name))//'.dat'
                print *, 'file: ', file_name
                open(13, file=file_name, status='replace')
                write(13,"(5A5)") "# L","BBB", "EBB", "EEB", "EEE"
                do l=1, lmax_bis
                        write(13,"(I4, 9E20.10)") l,  bispectrum_polar_order1(1:4, lb2)                                          
                end do
                close(13)
        end if 

	end subroutine  calculator_bispectrum
!---------------------------------------------------------------------------------------------------------------------------------

    subroutine threej_store(lb1, lb2, lb3, threej_m22, threej_20, threej_00, lmax_ea, lmax_bis)
	implicit none
		
	integer, parameter :: mcp=KIND(1.d0)
	real(mcp), parameter :: pi = 3.14159265358979323846264338328d0    
	INTEGER, intent(in) ::  lb1,lb2,lb3, lmax_ea, lmax_bis
	real(mcp), intent(out) ::  threej_m22(3, lmax_ea, 0:lmax_bis+lmax_ea+5), threej_20(3,lmax_ea, 0:lmax_bis+lmax_ea+5)
	real(mcp) , intent(out) :: threej_00(3,lmax_ea,0:lmax_bis+lmax_ea+5)
	! a for alpha related power spectra, teb for normal power spectra, bis for output bispectrum_polar in one dimenssion. ususlly cut lmax_ea=500, \
	! cut lmax_teb=2500, lmax_teb >= lmax_ea+ lmax_bis
	integer :: l, lb, l1, i
	! lb1, lb2, lb3 for bispectrum_polar indices, l1, l2, l3 for summing indices
	integer :: thread_ix, ix, IER
	integer :: botl, topl, nl
	real(mcp) lblist(3)
	real(mcp), allocatable ::   Wig3j_m22(:, :),Wig3j_20(:, :), wig3j_00(:,:)
	character(10):: current_time 
	
	! precalculated 3j arrays to store, 

	
	!$ integer  OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS
	!$ external OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS 

	thread_ix = 1
	!$ thread_ix = OMP_GET_MAX_THREADS()
	
	allocate(Wig3j_m22(0:lmax_bis+lmax_ea+5, thread_ix) )
	allocate(Wig3j_20(0:lmax_bis+lmax_ea+5, thread_ix) )
	allocate(Wig3j_00(0:lmax_bis+lmax_ea+5, thread_ix) )
		
	Wig3j_m22= 0
	Wig3j_20 = 0
	Wig3j_00 = 0
	
	lblist = (/lb1, lb2, lb3/)
	

	!$OMP PARALLEL DO DEFAULT(PRIVATE),  &
	!$OMP SHARED( lblist,  Wig3j_m22, Wig3j_20, Wig3j_00, threej_m22, threej_20, threej_00, lmax_ea), schedule(dynamic)      
			
	do l1 = 2, lmax_ea
		
		!$   thread_ix = OMP_GET_THREAD_NUM()+1
		
		do i=1, 3
			lb = lblist(i)
			botl = max(0, abs(lb-l1))
			topl =  lb+l1
			nl = topl- botl+1
			
			call DRC3JJ(dble(l1), dble(lb), dble(-2), dble(2), dble(botl), dble(topl), Wig3j_m22(botl:, thread_ix), &
				nl, IER)
			threej_m22(i, l1, botl:) = Wig3j_m22(botl:, thread_ix)
			! store 3j(lfree, l1, lb; 0, -2, 2)
	
			botl = max(2, abs(lb-l1))		
			topl =  lb+l1
			nl = topl- botl+1
			call DRC3JJ( dble(lb), dble(l1),  dble(2), dble(0), dble(botl), dble(topl), Wig3j_20(botl:, thread_ix),&
				nl, IER)
			threej_20(i, l1, botl:) = Wig3j_20(botl:, thread_ix)
			! store 3j(l1,  lfree, lb; 0, -2, 2)
							
	
			botl = max(0, abs(lb-l1))		
			topl =  lb+l1
			nl = topl- botl+1
			call DRC3JJ( dble(lb), dble(l1),  dble(0), dble(0), dble(botl), dble(topl), Wig3j_00(botl:, thread_ix),&
				nl, IER)
			threej_00(i, l1, botl:) = Wig3j_00(botl:, thread_ix)
			! store 3j(l1,  lfree, lb; 0, -2, 2)
                        
			end do 
	end do			 
	!$OMP END PARALLEL DO
        

	deallocate( Wig3j_m22)
	deallocate( Wig3j_20)
	deallocate( Wig3j_00)

	!call date_and_time(time = current_time)
	!print *, 'precalculate wigner 3j symbol completed at',  current_time
	
	end subroutine threej_store


	! containing double loops on six parts , and containing all six permutations, no need to  repeat calling this subroutine
	subroutine double_loop_polar( CPT_ang, lb1, lb2, lb3,  threej_m22, threej_20, threej_00, CEa, Caa, &
		 CEE, out_polar2, order1_polar, lmax_alpha, lmax_ea, lmax_bis, lmax_teb)

	implicit none
       
	integer, parameter :: mcp=KIND(1.d0)
	real(mcp), parameter :: pi = 3.14159265358979323846264338328d0    
	INTEGER, intent(in) :: lmax_alpha,  lmax_teb,lmax_bis,lmax_ea
	! a for alpha related power spectra, teb for normal power spectra, bis for output bispectrum_polar in one dimenssion. ususlly cut lmax_ea=500, \
	! cut lmax_teb=2500, lmax_teb >= lmax_ea+ lmax_bis
	real(mcp), intent(in) :: CEE(lmax_teb),  CEa(lmax_ea), Caa(lmax_ea)
	! input alpha related power spectra are all scale invariant form.
	! note there is no lensing contribution on CBB,
	real(mcp), intent(in):: CPT_ang
	real(mcp) cs2, sn2, cs2sq, sn2sq, cs4, sn4
	! triangle functions
	real(mcp), intent(in) ::  threej_m22(3, lmax_ea, 0:lmax_bis+lmax_ea+5), threej_20(3,lmax_ea, 0:lmax_bis+lmax_ea+5)
	real(mcp), intent(in) ::  threej_00(3, lmax_ea, 0:lmax_bis+lmax_ea+5)
	! degree unit	 
	real(mcp), intent(out) :: out_polar2(4), order1_polar(4)
	! two parts, one in two loops, another in three loops
	!geometrical factors, h(l1, l2,l3) and wigner 3j symbols I^{0-22}_{l1, l2, l3}, I^{0-22}_{l2, l3, l1}, I^{0-22}_{l3, l2, l1}
	! bispectrum_polar order,4:EEE, 3:EEB, 2:EBB, 1:BBB	
	integer :: i, j, j1, i1,i2,i3, j2, j3, IER, idd1, idd2, idd3
	! summation types indices
	integer :: l, lb, lb1, lb2, lb3, l1, l2, l3, thread_ix, l2start, l2end
	character(3) fields_list(4),  perm_list(6)
	character(3) :: perm_inds 
	integer :: lblist(3)
	integer :: perm_lb1,perm_lb2, perm_lb3
	! lb1, lb2, lb3 for bispectrum_polar indices, l1, l2, l3 for summing indices
	real(mcp) px(5, lmax_teb), xpol(6, 6), perm_xpol(6, 6)
	! xpol, first 6 is six parts, second 6 is for six perms
	real(mcp) xbbb(6), xebb(6), xeeb(6), xeee(6)
	! six for six parts
	real(mcp) , allocatable :: ppx(:, :)
	! first index for fields, second for permutations
	real(mcp) fac, tp_sum
	integer:: botl, topl, nl
	! precalculated 3j arrays to store, 
	!real(mcp) wig3j_000(0:2*lmax_teb)
	real(mcp) wig3j022_121(0:lmax_bis+lmax_ea+5), wig3j022_122(0:lmax_bis+lmax_ea+5), wig3j022_123(0:lmax_bis+lmax_ea+5), wig3j022_211(0:lmax_bis+lmax_ea+5)
	real(mcp) Wig3j000_123(0:lmax_bis+lmax_ea+5)
 	! for odd geometrical factor
	real(mcp) coeff_312
	character(10):: current_time
	! geometrical factor h(l1, l2, l3) and coeff I^{0, -2, 2}_{l1, l2, l3}, I^{0, -2, 2}_{l2, l1, l3}, I^{0, -2, 2}_{l3, l2, l2}
	
	!$ integer  OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS
	!$ external OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS 

	!cs2 = cos(2*CPT_ang)
	!sn2 = sin(2*CPT_ang)
	!cs2sq = cos(2*CPT_ang)**2
    	!sn2sq = sin(2*CPT_ang)**2
    	!cs4 = cos(4*CPT_ang)
	!sn4 = sin(4*CPT_ang)
        
        ! make small approximation when small than 0.5 deg to avoid trunck
        if(CPT_ang .lt. 0.0087)then
                ! five order approximations
                sn2 = 2*CPT_ang-4*CPT_ang**3/3._mcp+4*CPT_ang**5/15._mcp
                cs2 = 1-2*CPT_ang**2+2*CPT_ang**4/3._mcp
        else
                sn2 = sin(2*CPT_ang)
                cs2 = cos(2*CPT_ang)
        end if 
        cs2sq = cs2**2
        sn2sq = sn2**2
        cs4 = cs2sq-sn2sq
        sn4 = 2*sn2*cs2
 
	! loop of six permutations

	out_polar2 = 0
	perm_list = (/'123', '132', '231', '213', '321', '312'/)
	lblist = (/lb1, lb2, lb3/)
	
	thread_ix = 1
	!$ thread_ix = OMP_GET_MAX_THREADS()

	!print *, 'max thread_ix', thread_ix
	!print *, 'thread_ix', thread_ix
	allocate(ppx(5, thread_ix))
	xpol=0
	
	
	
	!$OMP PARALLEL DO DEFAULT(PRIVATE),  &
	!$OMP SHARED( perm_list, lblist, lmax_alpha, lmax_teb, lmax_ea), &
	!$OMP SHARED(threej_m22, threej_20, threej_00,  CEa, Caa, CEE, xpol), schedule(dynamic) 

	do j=1, 6
		perm_inds = trim(perm_list(j))
		read(perm_inds(1:1),'(I2)') idd1
		read(perm_inds(2:2),'(I2)') idd2
		read(perm_inds(3:3),'(I2)') idd3
	
		perm_lb1 = lblist(idd1)
		perm_lb2 = lblist(idd2)
		perm_lb3 = lblist(idd3)

		botl = max(2, abs(perm_lb1-perm_lb2))
		topl =  perm_lb1+perm_lb2
		nl = topl- botl+1
		! for 3j(l3,l1,l2;0,-2,2)
		coeff_312 = threej_m22(idd2, perm_lb1, perm_lb3)
		! 3j(lb3, lb1, lb2; 0, -2, 2)

		ppx = 0
		Wig3j022_121 = 0
		Wig3j022_122 = 0
		Wig3j022_123=0
		Wig3j000_123=0
		Wig3j022_211=0

		do l1 = 2, lmax_ea


			botl = max(2, abs(perm_lb1-l1))
			topl =  perm_lb1+l1
			nl = topl- botl+1
			
			Wig3j022_121(botl:) = threej_20(idd1, l1, botl:)
			! for 3j(L1,L2,lb1;0,-2,2)

			botl = max(2, abs(perm_lb2-l1))
			topl =  perm_lb2+l1
			nl = topl- botl+1
			Wig3j022_122(botl:) = threej_20(idd2, l1, botl:)
			! for 3j(L1,L2,lb2;0,-2,2)

			botl = max(2, abs(perm_lb3-l1))
			topl =  perm_lb3+l1
			nl = topl- botl+1
			Wig3j022_123(botl:) = threej_20(idd3, l1, botl:)
			! for 3j(L1,L2,lb3;0,-2,2)

			botl = max(0, abs(perm_lb3-l1))
			topl =  perm_lb3+l1
			nl = topl- botl+1
			Wig3j000_123(botl:) = threej_00(idd3, l1, botl:)
			! for 3j(L1,L2,lb3;0,0,0)

			!if(any(Wig3j000_123>HUge(1.d200)))then
			!	print *, 'perm_lb3, l1,, Wig3j000_123', perm_lb3, l1,  Wig3j000_123(:)	
			!	stop		
			!end if 	
			botl = max(0, abs(perm_lb1-l1))
			topl =  perm_lb1+l1
			nl = topl- botl+1
			Wig3j022_211(botl:) = threej_m22(idd1, l1, botl:)
			! for 3j(L2,L1,lb1;0,-2,2)


			l2start = max(2, min(abs(perm_lb1-l1), abs(perm_lb2-l1), abs(perm_lb3-l1)))
			l2end = min(lmax_teb, max(perm_lb1+l1, perm_lb2+l1, perm_lb3+l1))

			px = 0
			
			do l2 = l2start, l2end
				!first part 
				px(1, l2) = Wig3j022_122(l2)**2*CEE(l2)*(1-2*mod(l1+l2+perm_lb2, 2))*(2*l2+1)/4./pi	
				!the (2l2+1)/4pi factor from the I^{0-22} devided by 3j symbols			
				!second part
				px(2, l2) = Wig3j022_122(l2)**2*CEE(l2)*(2*l2+1)/4./pi	 
				!third part 
				if(l2 .le. lmax_ea)then
					px(3, l2) = Wig3j000_123(l2)*Wig3j022_123(l2)*CEa(l2)*(2*l2+1)/4./pi
					! four th part
					px(4, l2) = Wig3j022_211(l2)*Wig3j022_121(l2)*CEa(l2)*(2*l2+1)/4./pi
					! fifth part
					px(5, l2) = Wig3j022_211(l2)*Wig3j022_121(l2)*CEa(l2)*(1-2*mod(l1+l2+perm_lb1, 2))*(2*l2+1)/4./pi
				else
					px(3:5, l2)=0
				end if 	

				!if(px(3, l2)>Huge(1.d200))then
				!	print *, 'l1, l2, px3', l1, l2, px(3, l2), Wig3j000_123(l2), Wig3j022_123(l2), CEa(l2), pi
				!	stop
				!end if 
			end do 

			!$   thread_ix = OMP_GET_THREAD_NUM()+1
			
			!print *, 'thread_ix', thread_ix
			
			ppx(1,  thread_ix) =ppx(1,  thread_ix)+  sum(px(1, :))*Caa(l1)*(2*l1+1)
			ppx(2,  thread_ix) =ppx(2,  thread_ix) + sum(px(2, :))*Caa(l1)*(2*l1+1)
			ppx(3,  thread_ix) =ppx(3,  thread_ix) + sum(px(3, :))*Caa(l1)*(2*l1+1)
			ppx(4,  thread_ix) =ppx(4,  thread_ix) + sum(px(4, :))*CEa(l1)*(2*l1+1)
			ppx(5,  thread_ix) =ppx(5,  thread_ix) + sum(px(5, :))*CEa(l1)*(2*l1+1)

			!print *, 'ppx(3,  thread_ix), sum(px(3, :)), Caa(l1)', ppx(3,  thread_ix), sum(px(3, :)), Caa(l1)				


		end do 
					
		
		!print *, 'ppx(3, :)', ppx(3, thread_ix)

		! 1, 2,3,4,5 parts are related with loops, fist dimention is part, second is permutations
		if (perm_lb3 .le. lmax_ea)then
		    xpol(1, j) = sum(ppx(1, :))*CEa(perm_lb3)*coeff_312
		    xpol(2, j) = sum(ppx(2, :))*CEa(perm_lb3)*coeff_312
		    xpol(4, j) = sum(ppx(4, :))*CEa(perm_lb3)*coeff_312
		    xpol(5, j) = sum(ppx(5, :))*CEa(perm_lb3)*coeff_312	
		    !the last part are not involved in loop
		    xpol(6, j) = 0.25_mcp*CEE(perm_lb2)*CEa(perm_lb3)	
		end if 
		xpol(3, j) = sum(ppx(3, :))*CEE(perm_lb2)*coeff_312
		!print *, 'coeff_312', coeff_312

	end do 
	!$OMP END PARALLEL DO

	!fac = 4*exp(-6*corr0)/geom(lb2)*coeff_312

	!print *, ' xpol',  xpol(3, :)

	!parity even case
	if(mod(lb1+lb2+lb3, 2)==0)then
		! bbb
		! for different parts
		out_polar2(1) =  -2*cs2*sn2sq*(sum(xpol(1, : )) + sum(xpol(3, :))&
		 + sum(xpol(4, :)) )

                order1_polar(1) = 2*cs2*sn2sq*sum(xpol(6, :))
		
		!ebb
		xebb(6) = 2*(-sn2sq*sum( xpol(6, 1:2)) + cs2sq*sum(xpol(6, 3:6)) )
		xebb(1) = - 2*cs2sq*sum(xpol(1, [3,5])) - cs4*sum(xpol(1, [1,2,4,6])) 
		xebb(2) =  - sum(xpol(2,1:2))+ sum(xpol(2, [4,6]))
		xebb(3) = 2*( sn2sq*sum(xpol(3, [3,5]))  - cs2sq*sum(xpol(3, [1,2,4,6])))
		xebb(4) = 2*cs2sq*sum(xpol(4, [3,5]))-cs4*(sum(xpol(4, [1,2,4,6])))
		xebb(5) = sum(xpol(5,1:2)) - sum(xpol(5, [4,6])) 
		out_polar2(2) = sn2*sum(xebb(1:5))
                
                order1_polar(2) = sn2*xebb(6) 


		!eeb
		xeeb(6) = 2*( cs2sq*sum(xpol(6, [5,6]))- sn2sq*sum(xpol(6, 1:4)) )
		xeeb(1) = 2*sn2sq*sum(xpol(1, [1,4]))-cs4*( sum(xpol(2, [5,6])) - sum(xpol(1, [2,3])) )
		xeeb(2) = sum(xpol(2, [5,6]))- sum(xpol(2, [2,3]))
		xeeb(3) = 2*(  -cs2sq*sum(xpol(3, [1,4]))+sn2sq*sum(xpol(3, [2,3,5,6]))  )
		xeeb(4) = 2*sn2sq*sum(xpol(4, [1,4]))-cs4*sum(xpol(4, [2,3,5,6]))
		xeeb(5) =  sum(xpol(5, [2,3]))- sum(xpol(5, [5,6]))
		out_polar2(3) = cs2*sum(xeeb(1:5))
		
                order1_polar(3) = cs2*xeeb(6)
		
                
                !eee
		out_polar2(4)= 2*sn2*cs2sq*( sum(xpol(1, :)) +sum(xpol(3, :))+sum(xpol(4, :)) )
                
                order1_polar(4) = -2*sn2*cs2sq*sum(xpol(6, :))
	else
		
		! bbb 
                ! orders: 123, 132, 231, 213, 321, 312
		! single exchange change sign, 123, 231, 312 positive, negative:132, 213, 321
		! or here vonventions, positive:1, 3, 6, negative:2, 4, 5
		! ebb pair(2<->3) are: (1,2), (3,5), (4,6), eeb pair(1<->2): (1,4), (2,3), (5,6) 
		perm_xpol = xpol
		perm_xpol(:, [2,4,5]) = -xpol(:, [2,4,5])
		!wrong code , previous
                !out_polar2(1) = sn2*( 2*cs2sq*sum(perm_xpol(:, 3)) &
		!+ cs4*(sum(perm_xpol(:, 1)) - sum(perm_xpol(:, 4))) + sum(perm_xpol(:, 2))  - sum(perm_xpol(:, 5)) )
                
                !order1_polar(1) = sn2*2*sn2sq*sum(perm_xpol(:,6))
                ! corrected in 2020/4/27
                 out_polar2(1) = sn2*( 2*cs2sq*sum(perm_xpol(3, :)) &
                            + cs4*(sum(perm_xpol(1, :)) - sum(perm_xpol(4, :))) + sum(perm_xpol(2, :))  - sum(perm_xpol(5, :)) )
                 order1_polar(1) = sn2*2*sn2sq*sum(perm_xpol(6, :)) 

		!ebb
		!elemental one is 1<->2=213, 1<->3=321, the 2<->3 part change signs, 1,4, 5 invariant,  2,3,6 change sign
		perm_xpol = xpol
		perm_xpol(:, [2,3,6]) = -xpol(:, [2,3,6])
		xebb(6) = 2*sn2sq*(sum(perm_xpol(6, 1:2)) -sum(perm_xpol(6, 3:6))) 
		xebb(1) = -(cs4*sum(perm_xpol(1, [3,5]))+2*sn2sq*(sum(perm_xpol(1, 1:2)) - sum(perm_xpol(1, [4, 6])) ) )
		xebb(2) = -sum(perm_xpol(2, [3,5]))
		xebb(3) = -2*( cs2sq*sum(perm_xpol(3, [4,6])) + sn2sq*(sum(perm_xpol(3, 1:2)) - sum(perm_xpol(3, [3, 5])))  )
		xebb(4) = -cs4*sum(perm_xpol(4, [3,5]))+2*sn2sq*( sum(perm_xpol(4, 1:2)) - sum(perm_xpol(4, [4, 6]) ) )
		xebb(5) = -sum(perm_xpol(5, [3,5]))
		out_polar2(2) = cs2*sum(xebb(1:5))
                
                order1_polar(2) = cs2*xebb(6)

		!eeb
		!elemental one is 2<->3=132, 1<->3=321, the 1<->2 part change signs . 1, 2, 5 invariant, 3, 4,6 change sign
		perm_xpol = xpol
		perm_xpol(:, [3,4,6]) = -xpol(:, [3,4,6])
		xeeb(6) = 2*cs2sq*( sum(perm_xpol(6, [1,4])) -sum(perm_xpol(6, [2,3,5,6])) )
		xeeb(1) = -cs4*sum(perm_xpol(1, [1,4]))+2*cs2sq*sum(perm_xpol(1, [2,3,5,6]))
		xeeb(2) = sum(perm_xpol(2, [1,4]))
		xeeb(3) = 2*(  -sn2sq*sum(perm_xpol(3, [2,3])) +cs2sq*(sum(perm_xpol(3, [5,6])) -sum(perm_xpol(3, [1,4])) ) )
		xeeb(4) = cs4*sum(perm_xpol(4, [1,4]))-2*cs2sq*sum(perm_xpol(4, [2,3,5,6]))
		xeeb(5) = - sum(perm_xpol(5, [1,4]))
		out_polar2(3) = sn2*sum(xeeb(1:5))
                
                order1_polar(3) = sn2*xeeb(6)

		!eee
		! or here vonventions, positive:1, 3, 6, negative:2, 4, 5
		perm_xpol = xpol
		perm_xpol(:, [2,4,5]) = -xpol(:, [2,4,5])
		out_polar2(4) =cs2*(   2*sn2sq*sum(perm_xpol(3, :) )&
		+cs4*(-sum(perm_xpol(1, :) ) +sum(perm_xpol(4, :) )) +sum(perm_xpol(2, :) )-sum(perm_xpol(5, :) )	)
                
                order1_polar(4) = cs2*2*cs2sq*sum(perm_xpol(6, :) )
	end if 

	!do i=1, 4
	!	fac = 4*exp(-6*corr0)*geom(lb2)
		!or I can drop the factor here
	!	out_polar2[i, lb2] = out_polar2[i, lb2]*fac
	!end do 

	!call date_and_time(time = current_time)
	!print *, 'double loop calculation of polar bispectrum is completed at time',  current_time
	

	end subroutine double_loop_polar

	
	! containing triple loops on two types , and containing only one permutations because they are complex, &
	!  need to  repeat calling this subroutine
	subroutine triple_loop_polar(CPT_ang, tp, perm_lb1, perm_lb2, perm_lb3, perm_inds,  threej_m22, threej_20, CEa, Caa, &
		CEE, out_polar,lmax_alpha,  lmax_ea, lmax_bis, lmax_teb)

	implicit none

	integer, parameter :: mcp=KIND(1.d0)
	real(mcp), parameter :: pi = 3.14159265358979323846264338328d0
	character(len=3), intent(in)::   perm_inds
	integer , intent(in):: tp, perm_lb1,perm_lb2, perm_lb3
	real(mcp), intent(in) :: CPT_ang
	! perm_lb1, perm_lb2, perm_lb3 for bispectrum indices, l1, l2, l3 for summing indices
	INTEGER, intent(in) :: lmax_alpha, lmax_ea, lmax_bis, lmax_teb
	real(mcp), intent(in) :: threej_m22(3, lmax_ea, 0:lmax_bis+lmax_ea+5), threej_20(3, lmax_ea, 0:lmax_bis+lmax_ea+5)
	real(mcp), intent(in) :: CEE(lmax_teb), CEa(lmax_ea), Caa(lmax_ea)
	! a for alpha related power spectra, teb for normal power spectra, bis for output bispectrum in one dimenssion. ususlly cut lmax_ea=500, \
	! cut lmax_teb=2500
	real(mcp), intent(out) :: out_polar(4)
	real(mcp) ly(4, lmax_teb), llly(4, lmax_ea)
	! ly has  long dimentions due to EE
	real(mcp), allocatable :: lly(:, :)
	! summing arrays, make mpi on second loop
	integer :: i, j, j1, j2, j3
	integer :: l, l1, l2, l3, l2start, l2end, l3start, l3end
	integer :: id1, id2, id3
	! eod, even or odd judgement number, even:0, odd:1
	integer :: thread_ix, ix, IER
	integer :: botl, topl, nl, botl_6j, topl_6j, nl_6j
	real(mcp) fac, u1, fac_cs1, fac_sn1, fac_cs2, fac_sn2, fac_cs3, fac_sn3
        real(mcp) cs2, sn2, cs2sq, sn2sq, cs4, sn4
	real(mcp) ::   Wigloop_131(0:lmax_bis+lmax_ea+5 ),wigloop_322(0:lmax_bis+lmax_ea+5 ),Wigloop_213(0:lmax_bis+lmax_ea+5 )&
	,Wigloop_232(0:lmax_bis+lmax_ea+5 )
	real(mcp), allocatable ::   Wig6j(:,:)
	real(mcp) coe_eee, coe_bbb
	real(mcp) coe_ebb(3), coe_eeb(3), beta(3), snbeta(3), csbeta(3)
	character(10) ::current_time
	! one dimenssional wigner 6j array    

	!$ integer  OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS
	!$ external OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS 


        ! make small approximation when small than 0.5 deg to avoid trunck
        if(CPT_ang .lt. 0.0087)then
                sn2 = 2*CPT_ang-4*CPT_ang**3/3._mcp+4*CPT_ang**5/15._mcp
                cs2 = 1-2*CPT_ang**2+2*CPT_ang**4/3._mcp
        else
                sn2 = sin(2*CPT_ang)
                cs2 = cos(2*CPT_ang)
        end if 
        cs2sq = cs2**2
        sn2sq = sn2**2
        cs4 = cs2sq-sn2sq
        sn4 = 2*sn2*cs2
 

	!default satisfy the even triangles conditions
	read(perm_inds(1:1),'(I2)') id1
	read(perm_inds(2:2),'(I2)') id2
	read(perm_inds(3:3),'(I2)') id3
	
 

	!print*, 'eod2,eod3', eod2,eod3
	!print*, 'id1,id2,id3', id1,id2,id3

	!if(perm_lb1==0)then
	!    print*,'out loop:lb1,lb2,lb3', perm_lb1,perm_lb2,perm_lb3
	!   stop 
	!end if

	thread_ix = 1
	!$ thread_ix = OMP_GET_MAX_THREADS()	

	allocate(Wig6j(0:lmax_bis+lmax_ea+5, thread_ix) )
	allocate(lly(4, thread_ix))

	Wigloop_131 = 0
	Wigloop_213 = 0
	llly=0	


	!$OMP PARALLEL DO DEFAULT(PRIVATE),  &
	!$OMP SHARED(tp, perm_lb1, perm_lb2, perm_lb3, id1, id2, id3, lmax_alpha, lmax_teb, lmax_ea), &
	!$OMP SHARED(threej_m22, threej_20, CPT_ang, cs2, sn2, CEa, Caa, CEE, llly), schedule(dynamic)      

	do l1=2, lmax_ea  

	
		! calculate 3j (l2, l1, perm_lb2;0, -2, 2), l2 as free index 
		! be careful that  l2 can be less than 2 in this 3j symbol, we should pass those l2 less than 2.


		botl = max(2, abs(perm_lb1-l1))		
		topl =  perm_lb1+l1
		nl = topl- botl+1
		Wigloop_131(botl:topl ) = threej_20(id1, l1, botl:topl)
		! calculate 3j (l1, l3, perm_lb1;0, -2, 2), l3 is free index
		
		botl = max(0, abs(perm_lb3-l1))
		topl = perm_lb3+l1
		nl= topl- botl+1
		Wigloop_213(botl:topl ) = threej_m22(id3, l1, botl:topl)
		! calculate 3j (l2, l1, perm_lb3;0, -2, 2), l3 is free index	


		l2start = max(2,abs(perm_lb3-l1) )
                if(tp==1)then
		        l2end = min(lmax_ea, perm_lb3+l1)
                else
                        l2end = min(lmax_alpha, perm_lb3+l1)
                end if 
		! satisfying even or odd conditions, l2 only dependends on l1, perm_lb2.
		! add eod2, for transformation from B to E.	

		wigloop_322 = 0
		Wigloop_232 = 0
		Wig6j=0
		lly=0

		do l2 = l2start, l2end
		! l2 for index of alphaalpha power spectrum.	

			!print*,'inner loop: lb1,lb2,lb3', perm_lb1,perm_lb2,perm_lb3

			

			topl = perm_lb2+l2
			botl = max(0, abs(perm_lb2-l2))
			nl = topl- botl+1
			wigloop_322(botl:) = threej_m22(id2, l2, botl:)
			!3j(l3, l2, lb2; 0, -2, -2)
			
			botl = max(2, abs(perm_lb2-l2))
			nl = topl- botl+1
			Wigloop_232(botl:) = threej_20(id2, l2, botl:)	
			!3j(l2, l3, lb2; 0, -2, -2)
			

			botl_6j = max( abs(l1-perm_lb1), abs(l2-perm_lb2))
			topl_6j = min(l1+perm_lb1, l2+perm_lb2)
			nl_6j = topl_6j- botl_6j+1

			!$   thread_ix = OMP_GET_THREAD_NUM()+1

			call DRC6J (dble(l2), dble(perm_lb2), dble(perm_lb3), dble(perm_lb1), dble(l1), dble(botl_6j), dble(topl_6j), &
			Wig6j(botl_6j:, thread_ix), nl_6j, IER)
			! calculate wigner 6j {perm_lb1, perm_lb2, perm_lb3;l2, l1, l3} = {l3, l2, perm_lb2;perm_lb3, perm_lb1, l1},&
			!	 l3 as free index

			
			botl_6j = max(2, botl_6j)
			! enphonsize l3 is great tha 2, but this is not necessary in 6j.
						
								
			l3start = max(2, botl_6j) 
			if(tp==1)then
				l3end = min(topl_6j, lmax_ea)
			else
				l3end = min(topl_6j, lmax_teb)
			end if 

			beta(3) = 0.5_mcp*(l1+l2+perm_lb3)*pi-2*CPT_ang
                        
                        ! cos(0.5_mcp*(l1+l2+perm_lb3)*pi)
                        fac_cs3 = (1-mod(l1+l2+perm_lb3, 2))*(1- mod(l1+l2+perm_lb3, 4))
                        ! sin(0.5_mcp*(l1+l2+perm_lb3)*pi)
                        fac_sn3 = (1-mod(l1+l2+perm_lb3-1, 2))*(1- mod(l1+l2+perm_lb3-1, 4))

                        csbeta(3) = fac_cs3*cs2+fac_sn3*sn2
                        snbeta(3) = fac_sn3*cs2-fac_cs3*sn2
			
                        ly=0
			do l3 = l3start, l3end
			! l3  is for EE power spectrum indices, cut is lmax_teb

				!angles  
				beta(1) = 0.5_mcp*(l1+l3+perm_lb1)*pi-2*CPT_ang
				! cos(0.5_mcp*(l1+l2+perm_lb3)*pi)
                                fac_cs1 = (1-mod(l1+l3+perm_lb1, 2))*(1- mod(l1+l3+perm_lb1, 4))
                                ! sin(0.5_mcp*(l1+l2+perm_lb3)*pi)
                                fac_sn1 = (1-mod(l1+l3+perm_lb1-1, 2))*(1- mod(l1+l3+perm_lb1-1, 4))

                                csbeta(1) = fac_cs1*cs2+fac_sn1*sn2
                                snbeta(1) = fac_sn1*cs2-fac_cs1*sn2
                                
                                
                                beta(2) = 0.5_mcp*(l2+l3+perm_lb2)*pi-2*CPT_ang
				! cos(0.5_mcp*(l1+l2+perm_lb3)*pi)
                                fac_cs2 = (1-mod(l2+l3+perm_lb2, 2))*(1- mod(l2+l3+perm_lb2, 4))
                                ! sin(0.5_mcp*(l1+l2+perm_lb3)*pi)
                                fac_sn2 = (1-mod(l2+l3+perm_lb2-1, 2))*(1- mod(l2+l3+perm_lb2-1, 4))

                                csbeta(2) = fac_cs2*cs2+fac_sn2*sn2
                                snbeta(2) = fac_sn2*cs2-fac_cs2*sn2
				
				
				!angles coefficients of BBB, EEE
				coe_bbb = csbeta(1)*csbeta(2)*csbeta(3)
				coe_eee = snbeta(1)*snbeta(2)*snbeta(3)

				!angles coefficients of EBB
				coe_ebb(1) =  snbeta(1)*csbeta(2)*csbeta(3)
				coe_ebb(2) =  csbeta(1)*snbeta(2)*csbeta(3)
				coe_ebb(3) =  csbeta(1)*csbeta(2)*snbeta(3)

 				!angles coefficients of EEB
				coe_eeb(1) =  snbeta(1)*snbeta(2)*csbeta(3)
				coe_eeb(2) =  csbeta(1)*snbeta(2)*snbeta(3)
				coe_eeb(3) =  snbeta(1)*csbeta(2)*snbeta(3)

				!type 1 , or part 1 for 2 permutations, type 2 for 5 permutations
				if (tp==1)then
					u1 =  Wig6j(l3, thread_ix)*Wigloop_131(l3)*wigloop_322(l3)*CEa(l3)*(2*l3+1)/4./pi
				else
					u1 =  Wig6j(l3, thread_ix)*Wigloop_131(l3)*wigloop_232(l3)*CEE(l3)*(1-2*mod(l3, 2)) *(2*l3+1)/4./pi
				end if 
				!(2l3+1)/4/pi is from I^{0-22} devided by 3j symbols, the sqrt((2lb1+1)(2lb2+1)(2lb3+1)/4pi)&
				!  factor is elimited with geometrical factor.  

				ly(1, l3) = coe_bbb*u1
				ly(4, l3) = coe_eee*u1

				if (tp==1)then	
					ly(2, l3) = coe_ebb(1)*u1
					ly(3, l3) = coe_eeb(1)*u1
				else
					! orders, 123, 132, 231, 213, 321, 312
					if (perm_inds=='123')then
						ly(2, l3) =  coe_ebb(1)*u1
						ly(3, l3) =  coe_eeb(1)*u1
					else if(perm_inds=='132')then
						ly(2, l3) =  coe_ebb(1)*u1
						ly(3, l3) =  coe_eeb(3)*u1
					else if(perm_inds=='231')then
						ly(2, l3) =  coe_ebb(3)*u1
						ly(3, l3) =  coe_eeb(3)*u1
					else if(perm_inds=='213')then
						ly(2, l3) =  coe_ebb(2)*u1
						ly(3, l3) =  coe_eeb(1)*u1
					else if(perm_inds=='321')then
						ly(2, l3) =  coe_ebb(3)*u1
						ly(3, l3) =  coe_eeb(2)*u1
					else
						ly(2, l3) =  coe_ebb(2)*u1
						ly(3, l3) =  coe_eeb(2)*u1
					end if 
				end if 	
			end do 
			
			!$   thread_ix = OMP_GET_THREAD_NUM()+1

			do i=1, 4
				if (tp==1)then
					lly(i, thread_ix) = lly(i, thread_ix) + sum(ly(i, :))*Wigloop_213(l2)*CEa(l2)*(2*l2+1)
				else 
					lly(i, thread_ix) = lly(i, thread_ix)+ sum(ly(i, :))*Wigloop_213(l2)*Caa(l2)*(1-2*mod(l2+perm_lb1+perm_lb3, 2))*(2*l2+1)  
				end if 
			end do 

		end do	
		
		do j=1, 4
			if (tp==1)then
				llly(j, l1) = sum(lly(j, :))*CEa(l1)*(1-2*mod(perm_lb1+perm_lb2+perm_lb3, 2))*(2*l1+1)
			else 
				llly(j, l1) = sum(lly(j, :))*CEa(l1)*(2*l1+1) 
			end if 
		end do 

	end do 
	!$OMP END PARALLEL DO



	deallocate( Wig6j)
	deallocate(lly)

	do i=1, 4
		! no 8*e^{-6corr0}*h^{-1}(l1, l2, l3) factor 
		out_polar(i) = sum(llly(i, :))
	end do 

	!call date_and_time(time = current_time)
	!print *, 'one triple loop calculation of polar bispectrum has been done at time',  current_time

	end subroutine triple_loop_polar 


