    ! command
    ! f2py -c -m cpt_transform cpt_transform.f90 --opt='-O3' --fcompiler=intelem --f90flags="-qopenmp -D__OPENMP" --f77flags="-qopenmp -D__OPENMP" -liomp5
    !add cpt transformation inherited from camb/lensing.f90 by me
    
    !Program test
    !implicit none

    !integer, parameter :: lmax_teb = 2500
    !integer lmax_ea
    !real(8)  CPT_ang, corr0, delta_CPT_ang(6)
    !real(8)  CEE(2:lmax_teb), CBB(2:lmax_teb), CTE(2:lmax_teb)
    !real(8) :: CEE_new(2:lmax_teb), CBB_new(2:lmax_teb), CTE_new(2:lmax_teb), CTB_new(2:lmax_teb), CEB_new(2:lmax_teb)



    !CPT_ang = 1.0
    !corr0 = 0.02
    !delta_CPT_ang = (/1.d-8, 1.d-8, 1.d-8, 1.d-8,1.d-8, 1.d-8/)
    !!lmax_teb = 2500
    !lmax_ea = 2*lmax_teb

    !CTE(2:lmax_teb) = 0
    !CEE(2:lmax_teb) = 0
    !CBB(2:lmax_teb) = 0


    !CALL cpt_transform( CPT_ang, corr0, delta_CPT_ang, CEE, CBB, CTE, CEE_new, CBB_new, CTE_new, CTB_new, CEB_new, lmax_teb, lmax_ea)

    !PRINT *, 'CEE_new', CEE_new

    !end PROGRAM test     
    !---------------------------------------------------------------------------------------------------------------------------------------
    SUBROUTINE temp_bispectrum( eod, CPT_ang, CEE, CTE, Caa, CEa, CTa, bispectrum_temp, lmax_alpha, lmax_ea, lstart, lmax_bis, lmax_teb, lmax_computed_cl)
    
    implicit none
    
    integer, parameter :: mcp=KIND(1.d0)
    real(mcp), parameter :: pi = 3.14159265358979323846264338328d0
    
    INTEGER, intent(in) :: lmax_alpha,  lmax_teb, lmax_ea, lstart, lmax_bis, lmax_computed_cl
    real(mcp), intent(in) :: CPT_ang
    character(4), intent(in) :: eod! even or odd
    real(mcp), intent(in) :: CEE(lmax_teb),  CTE(lmax_teb), Caa(lmax_ea), CEa(lmax_ea), CTa(lmax_ea)
    real(mcp), intent(out) :: bispectrum_temp(5, lmax_bis)
    ! orders, TTB, TTE, TBB, TEB, TEE
    integer :: i, j, l, l1,  l2, npoints, ntime
    integer :: lb1, lb2, lb3
    real(mcp) AccuracyLevel, AccuracyBoost
    real(mcp) theta, dtheta,llp1,fac, fac1,fac2,fac3
    real(mcp), allocatable :: P(:)
    real(mcp)  dP(lmax_teb)
    real(mcp) sinth,halfsinth, x
    real(mcp) lfacs(lmax_teb), lfacs2(lmax_teb), lrootfacs(lmax_teb)
    real(mcp) d_11(lmax_teb),d_m11(lmax_teb), d_22(lmax_teb),d_2m2(lmax_teb),d_20(lmax_teb)
    real(mcp), allocatable :: Wig3j_000(:,:), Wig3j_112(:,:), Wig3j_121(:, :), Wig3j_211(:, :)
    real(mcp), allocatable :: Wig3j_123(:,:), Wig3j_231(:,:), Wig3j_321(:,:)
    real(mcp), allocatable ::  wig3_202(:, :), wig2_202(:,:), wig3_022(:,:), wig2_022(:, :), wig3_000(:, :), wig2_000(:, :)
    integer :: topl, botl, nl, IER, l2start, l2end
    real(mcp) geom(lmax_bis), ifac_123(lmax_bis), ifac_231(lmax_bis), ifac_321(lmax_bis)
    ! precalculated 3j symbol in specific value.
    real(mcp) xtemp(4)
    ! stands for 4 general parts of bispectrum 
    real(mcp)  y_tbb(4), yy_tbb(4), yyy_tbb(4)
    ! approximation for TBB
    real(mcp), allocatable :: xtbb(:,:), xteb(:,:), xtee(:,:)
    ! part of bispectrum_temp

    integer thread_ix, interp_fac, apodize_point_width
    real(mcp) pmm, pmmp1
    logical :: short_integral_range
    real(mcp) range_fac
    real(mcp) theta_cut(lmax_teb) 
    
    real(mcp) corr0_CPT_ang, corr_CPT_ang
    real(mcp) cs2, sn2, cs4, sn4, corr0, exp4_corr0, exp8_corr
    real(mcp)  Caa_corr0(lmax_ea), Caa_corr(lmax_ea)
    real(mcp) px, mx, px1, px2, px3, px4, ppx1, ppx2, pppx1, pppx2        
    real(mcp) pwd, mwd, pcd, mcd, pdd, mdd
    character(10) :: current_time

        ! a check of this cpt calculation module when all delta_CPT_ang is zero, then the new cls has a linear relation with original  
        ! here we write the old cls in file
        ! only used in action=4
   
    !$ integer  OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS
    !$ external OMP_GET_THREAD_NUM, OMP_GET_MAX_THREADS 
    
    
    !lmax_ea = 2700

    !cs2 = cos(2*CPT_ang)
    !sn2 = sin(2*CPT_ang)
    !cs4 = cos(4*CPT_ang)
    !sn4 = sin(4*CPT_ang)
    if(CPT_ang .lt. 0.0087)then
            ! five order approximation
            !sn2 = 2*CPT_ang-4*CPT_ang**3/3._mcp+4*CPT_ang**5/15._mcp
            !cs2 = 1-2*CPT_ang**2+2*CPT_ang**4/3._mcp
            ! changed in 2020/4/27, lower order to decrease oscillations
            sn2 = 2*CPT_ang-4*CPT_ang**3/3._mcp
            cs2 = 1-2*CPT_ang**2
    else
            sn2 = sin(2*CPT_ang)
            cs2 = cos(2*CPT_ang)
    end if
    cs4 = cs2**2-sn2**2
    sn4 = 2*sn2*cs2
	!if (lmax_bis.gt.50)then
	!	lstart = lmax_bis-50
	!else
	!	lstart=3
	!end if 

    !lstart=3
    !default from 3, we can change from other values

    ! here use  cl^alpha, not Dl^alpha
    ! variance, lcut at lmax_alpha
    do l=2, lmax_alpha
        Caa_corr0(l) = 0.25_mcp*( 2*real(l) + 1 )/pi*Caa(l)
 
    end do    
     
    corr0 =  sum(Caa_corr0(2:lmax_alpha))
    
    ! small order approximation for variance smaller than 1e-4 rad^2
    if(corr0 .le. 1e-4)then
        exp4_corr0  = 1-4*corr0
    else
        exp4_corr0  = exp(-4*corr0)
    end if 
    !lmax_eall = max(lmax_ea, lmax_teb, lmax_bis)
    thread_ix = 1
    !$ thread_ix = OMP_GET_MAX_THREADS()
    
    allocate(Wig3j_000(0:2*lmax_bis+5, thread_ix))
    allocate(Wig3j_112(0:2*lmax_bis+5, thread_ix))
    allocate(Wig3j_121(0:2*lmax_bis+5, thread_ix))
    allocate(Wig3j_211(0:2*lmax_bis+5, thread_ix))
    allocate(Wig3j_123(0:2*lmax_bis+5, thread_ix))
    allocate(Wig3j_321(0:2*lmax_bis+5, thread_ix))
    allocate(Wig3j_231(0:2*lmax_bis+5, thread_ix))

    Wig3j_000= 0
    Wig3j_112= 0
    Wig3j_121= 0
    Wig3j_211= 0
    Wig3j_123= 0
    Wig3j_321= 0
    Wig3j_231= 0
    
    geom = 0
    ifac_123= 0
    ifac_231=0
    ifac_321=0

    bispectrum_temp=0

    !$OMP PARALLEL DO DEFAULT(PRIVATE),  &
    !$OMP SHARED( eod, lstart, geom, ifac_123, ifac_321, ifac_231, Wig3j_123, Wig3j_321, Wig3j_231), &
	!$OMP SHARED( lmax_ea, lmax_bis,  Wig3j_000,  Wig3j_112, Wig3j_121, Wig3j_211), schedule(dynamic)      
    
    ! precalculate 3j symbols
    do l = lstart, lmax_bis

        !$   thread_ix = OMP_GET_THREAD_NUM()+1
        
        lb1=4
        lb2=l
        if(eod=='even')then
                lb3=l+4
        else
                lb3 = l+3
        end if 
        topl =  lb1+lb3
        if(eod=='even')then
            botl = max(0, abs(lb1-lb3))
            nl = topl- botl + 1
            !3j(lb1,lb2,lb3; 0,0,0)
            call DRC3JJ(dble(lb1), dble(lb3), dble(0), dble(0), dble(botl), dble(topl), Wig3j_000(botl:, thread_ix ), &
            nl, IER)
            geom(lb2 ) = Wig3j_000(lb2, thread_ix )
        else
            botl = max(1, abs(lb1-lb3))
            nl = topl- botl + 1
            !3j(lb1,lb2,lb3;1,1,-2)
            call DRC3JJ(dble(lb3), dble(lb1), dble(-2), dble(1),dble(botl), dble(topl), Wig3j_112(botl:,thread_ix ), &
            nl, IER)

            botl = max(2, abs(lb1-lb3))
            nl = topl- botl + 1
            !3j(lb1,lb2,lb3;1,-2,1)
            call DRC3JJ(dble(lb3), dble(lb1), dble(1), dble(1), dble(botl), dble(topl), Wig3j_121(botl:,thread_ix ), &
            nl, IER)

            botl = max(1, abs(lb1-lb3))
            nl = topl- botl + 1
            !3j(lb1,lb2,lb3;-2,1,1)
            call DRC3JJ(dble(lb3), dble(lb1), dble(1), dble(-2), dble(botl), dble(topl), Wig3j_211(botl:,thread_ix ), &
            nl, IER)

            geom(lb2 ) = (Wig3j_112(lb2, thread_ix ) + Wig3j_121(lb2,thread_ix )+Wig3j_211(lb2,thread_ix ))/3.

        end if 
        
        !print *, 'geom', geom(lb2 )
        topl =  lb1+lb3

        botl = max(2, abs(lb1-lb3))
        nl = topl- botl + 1
        !3j(lb1,lb2, lb3; 0, -2, 2)
        call DRC3JJ(dble(lb3), dble(lb1), dble(2), dble(0), dble(botl), dble(topl), Wig3j_123(botl:, thread_ix), &
         nl, IER)
         ifac_123(lb2 ) = Wig3j_123(lb2, thread_ix )

         botl = max(2, abs(lb1-lb3))
         nl = topl- botl + 1
        !3j(lb3,  lb2,  lb1; 0,-2,  2)
         call DRC3JJ(dble(lb1), dble(lb3), dble(2), dble(0), dble(botl), dble(topl), Wig3j_321(botl:, thread_ix), &
         nl, IER)
         ifac_321(lb2) = Wig3j_321(lb2, thread_ix )

         botl = max(0, abs(lb1-lb3))
         nl = topl- botl + 1
        !3j(lb2, lb3, lb1; 0, -2, 2)
         call DRC3JJ(dble(lb3), dble(lb1), dble(-2), dble(2), dble(botl), dble(topl), Wig3j_231(botl:, thread_ix), &
         nl, IER)
         ifac_231(lb2 ) = Wig3j_231(lb2, thread_ix )


    end do     
    !$OMP END PARALLEL DO
    


     !print *, 'geom', geom
    !bispectrum TTB, TTE
    do lb2 = lstart, lmax_bis
        lb1 = 4
        if(eod=='even')then
                lb3=lb2+4
        else
                lb3 = lb2+3
        end if 
        
        if(corr0 .le. 1e-4)then
        	fac = 2*(1-2*corr0)/geom(lb2)
        else
                fac = 2*exp(-2*corr0)/geom(lb2)
        end if 

       	if (lb2 .le. lmax_ea)then
		if (eod=='even')then
                    !TTB
		    bispectrum_temp(1, lb2) = cs2*(ifac_231(lb2)*CTE(lb1)*CTa(lb2)+ ifac_123(lb2)*CTE(lb2)*CTa(lb1) )*fac
		    !TTE
		    bispectrum_temp(2, lb2) = -sn2*(ifac_231(lb2)*CTE(lb1)*CTa(lb2)+ ifac_123(lb2)*CTE(lb2)*CTa(lb1) )*fac
		else
                    !TTB
		    bispectrum_temp(1, lb2) = sn2*(ifac_231(lb2)*CTE(lb1)*CTa(lb2)- ifac_123(lb2)*CTE(lb2)*CTa(lb1) )*fac
		    !TTE
		    bispectrum_temp(2, lb2) = cs2*(ifac_231(lb2)*CTE(lb1)*CTa(lb2)- ifac_123(lb2)*CTE(lb2)*CTa(lb1) )*fac
		end if 
	end if 
    end do 

    !print *, 'bispectrum_temp(2)', bispectrum_temp(2, :)
    deallocate(Wig3j_000)
    deallocate(Wig3j_112)
    deallocate(Wig3j_121)
    deallocate(Wig3j_211)
    deallocate(Wig3j_123)
    deallocate(Wig3j_321)
    deallocate(Wig3j_231)
    

    print *, 'wigner 3j completed!'
    
    !!----------------------------------------------------------------------------------------------------------------------------------
    
    ! default, the upper limit 
    !lmax_computed_cl = lmax_computed_cl
    AccuracyLevel = 2._mcp
        
    ! if I choose the ntimes number too large, the time will become large, here I choose the same as lensing.f90, ntime = 2
    ntime = 2
    npoints = lmax_computed_cl * ntime *AccuracyLevel
    short_integral_range = .false.!.not. this%accurate_BB 
    dtheta = pi / npoints
    if (lmax_computed_cl > 3500) dtheta=dtheta/1.3_mcp
    apodize_point_width = nint(0.003 / dtheta)
    npoints = int(pi/dtheta)
    if (short_integral_range) then
        range_fac= max(1._mcp,32/AccuracyLevel ) !fraction of range to integrate
        npoints = int(npoints /range_fac)
            !OK for TT, EE, TE but inaccurate for low l BB
            !this induces high frequency ringing on very small scales
            !which is then mitigated by the apodization below
    else
        range_fac=1
    end if

    AccuracyBoost = 2._mcp
    interp_fac = max(1,min(nint(10._mcp/AccuracyBoost),int(range_fac*2)-1))
        
    do l=2,lmax_teb
        lfacs(l) = real(l*(l+1),mcp)
        lfacs2(l) = real((l+2)*(l-1),mcp)
        lrootfacs(l) = sqrt(lfacs(l)*lfacs2(l))
    end do
    do l=2,lmax_teb
        theta_cut(l) = 0.244949_mcp/sqrt(3._mcp*lfacs(l) - 8._mcp)
    end do
  
    
    thread_ix = 1
    !$ thread_ix = OMP_GET_MAX_THREADS()
    
    allocate(xtbb(lmax_bis, thread_ix))
    allocate(xteb(lmax_bis, thread_ix))
    allocate(xtee(lmax_bis, thread_ix))
    xtbb=0
    xteb=0
    xtee=0
    
    


    ! only applied for EE,BB,EB, which need integration 

    !uncomment second line for PGF90 workaround
    !$OMP PARALLEL DO DEFAULT(PRIVATE),  &
    !OMP PRIVATE(P,dP,d11,dm11,d22,d2m2,d_20,corrcontribs,ddcontribs),&
    !$OMP SHARED(eod, lstart, lfacs,lfacs2,lrootfacs, lmax_alpha, lmax_teb, lmax_ea, lmax_bis, theta_cut), &
    !$OMP SHARED( Caa,CEa, CTa, CEE, CTE, cs4, sn4, corr0, exp4_corr0, geom, ifac_123, ifac_231, ifac_321),&
    !$OMP SHARED(dtheta,npoints,interp_fac,short_integral_range,apodize_point_width, xtbb, xteb, xtee)        

        
    do  i=1, npoints-1

        !print *, 'point ', i
        
        theta = i * dtheta
        x = cos(theta)
        sinth = sin(theta)
        halfsinth = sinth/2

        pmm=1
        pmmp1=x
        
        ! legendre polynomal
        allocate(P(max(lmax_teb, lmax_ea)))
        
        d_11(1) = cos(theta/2)**2
        d_m11(1) = sin(theta/2)**2
        P(1) = x
        d_22(1)=0
        d_2m2(1)=0
        d_20(1)=0
        

        ! calculate the wigner small d matrix values.
        do l=2,lmax_teb

            P(l)= ((2*l-1)* x *pmmp1 - (l-1)*Pmm)/ l
            dP(l) = l*(pmmp1-x*P(l))/sinth**2
            Pmm=pmmp1
            pmmp1=P(l)
            llp1 = lfacs(l)

            fac1 = (1-x)
            fac2 = (1+x)
            fac = fac1/fac2

            d_11(l) =  fac1*dP(l)/llp1 + P(l)
            d_m11(l) = fac2*dP(l)/llp1 - P(l)

            d_22(l) = ( ((4*x-8)/fac2 + llp1)*P(l) &
                + 4*fac*( fac2 + (x - 2)/llp1)*dP(l) )/ lfacs2(l)

            !For small theta use Taylor expansion for better stability (thanks Pavel Motloch)
            if (theta > theta_cut(l)) then
                d_2m2(l) = ( (llp1- (4*x+8)/fac1) *P(l) &
                    +4/fac*( -fac1 + (x+2)/llp1) *dP(l) )/lfacs2(l)
            else
                d_2m2(l) = lfacs(l)*lfacs2(l)*theta**4*(1._mcp/384._mcp &
                    - (3._mcp*lfacs(l) - 8._mcp)/23040._mcp*theta**2)
            endif

            d_20(l) = (2*x*dP(l) - llp1*P(l) ) / lrootfacs(l)

        end do    
        if (lmax_teb<lmax_ea) then
            do l= lmax_teb +1, lmax_ea
                P(l)= ((2*l-1)* x *pmmp1 - (l-1)*Pmm)/ l
                Pmm=pmmp1
                pmmp1=P(l)  
                !print *, 'x, l, P(l, x)',x, l, P(l)           
            end do 
        end if 


        !calculate C^\alpha(\beta)
        do l=2, lmax_alpha
            Caa_corr(l) =  0.25_mcp*( 2*real(l) + 1 )/pi*Caa(l)*P(l)
        end do 
        corr_CPT_ang = sum(Caa_corr(2:lmax_alpha))
        
        deallocate(P)

        !print *, 'loop of one point', i

        !$   thread_ix = OMP_GET_THREAD_NUM()+1   

        !print *, 'thread_ix', thread_ix
        
        ! start calculate bispectrum_temp
        ! TBB, TEB, TEE have four parts, xtbb, xteb, xtee
        


        do lb2 = lstart, lmax_bis
            ! specify lb1+1 = lb2= lb3-1, so lb2 is even number to satisfy lb1+lb2+lb3=even selection rule	
            ! usually select lmax_bis = 2000
            !lb2 = l
            lb1 = 4
            if(eod=='even')then
                        lb3 = lb2+4
                else
                        lb3 = lb2+3
                end if 
            ! precalculate the 3j 6j arrays

            px1=0
            px2=0
            px3=0
            px4=0
            ! four module parts 
            do l = 2, lmax_teb
		if (l .le. lmax_ea)then
        	    px1 = px1 +  d_20(l)*(2*l+1.)/2.*CEa(l)
		    px4 = px4 +  (2*l+1.)*d_20(l)*CEa(l)
		end if 

                px2 = px2 + d_2m2(l)*(2*l+1.)/2.*CEE(l)
                ! this part is pure odd
                px3 = px3 + d_22(l)*(2*l+1.)/2.*CEE(l)
    
            end do 

            px4 = 0.5_mcp*px4**2/pi
            ! the double summation can be transformed as square of single summation

            ! integration factor
            if(corr0.le. 1e-4)then
                 fac3 = sinth*(1-4*corr_CPT_ang+8*corr_CPT_ang**2)
                 exp8_corr = 1+8*corr_CPT_ang+ 32*corr_CPT_ang**2
            else
                 fac3 = sinth*exp(-4*corr_CPT_ang)
                 exp8_corr = exp(8*corr_CPT_ang)
             end if 
            ! four parts
            xtemp = 0
            ! four module functions
            xtemp( 1) =   px1*CTE(lb1)*fac3
            if (lb1.le. lmax_ea)then
		    xtemp( 2) =   px2*ifac_123(lb2)*CTa(lb1)*fac3 
		    xtemp( 3) =   px3*ifac_123(lb2)*CTa(lb1)*exp8_corr*fac3 
		    xtemp( 4) =   px4*ifac_123(lb2)*CTa(lb1)*fac3
	    end if 
            !theta factors were put in earlier (already in corr)
            
            !if (i>npoints-apodize_point_width*3) &
            !    xtemp=xtemp*exp(-(i-npoints+apodize_point_width*3)**2/real(2*apodize_point_width**2))

            mwd=ifac_321(lb2)*d_20(lb3) - ifac_231(lb2)*d_20(lb2)
            pwd=ifac_321(lb2)*d_20(lb3) + ifac_231(lb2)*d_20(lb2)
            mcd = d_2m2(lb3) - d_2m2(lb2)
            pcd = d_2m2(lb3) + d_2m2(lb2)
            mdd = d_22(lb3) - d_22(lb2)
            pdd = d_22(lb3) + d_22(lb2)

            ! parity even
            if(eod=='even')then
                ! tbb
                xtbb(lb2, thread_ix) =xtbb(lb2, thread_ix) + xtemp( 1)*sn4*pwd &
                + xtemp( 2)*sn4*pcd -  xtemp( 4)*sn4*pcd

                !teb
                xteb(lb2, thread_ix) = xteb(lb2, thread_ix) + xtemp( 1)*(cs4*pwd - exp8_corr*mwd) &
                +xtemp( 2)*cs4*pcd - xtemp( 3)*mdd- xtemp( 4)*(cs4*pcd + exp8_corr*mdd)

                !tee
                xtee(lb2, thread_ix) = xtee(lb2, thread_ix) - xtemp( 1)*sn4*pwd &
                - xtemp( 2)*sn4*pcd + xtemp( 4)*sn4*pcd                

           else! parity odd
                !tbb
                xtbb(lb2, thread_ix) = xtbb(lb2, thread_ix) + xtemp( 1)*(cs4-exp8_corr)*mwd &
                - xtemp( 2)*cs4*mcd + xtemp( 3)*mdd -  xtemp( 4)*(cs4*mcd+exp8_corr*mdd)
                !teb
                xteb(lb2, thread_ix) = xteb(lb2, thread_ix) - xtemp( 1)*sn4*mwd &
                +xtemp( 2)*sn4*mcd- xtemp(4)*sn4*mcd
                
                !tee
                xtee(lb2, thread_ix) = xtee(lb2, thread_ix) - xtemp( 1)*(cs4+exp8_corr)*mwd &
                + xtemp( 2)*cs4*mcd + xtemp( 3)*mdd+ xtemp( 4)*(-cs4*mcd + exp8_corr*mdd)   

           end if  

           
        end do        
        
        !print *, 'one point is completed!'

        
    end do
    !$OMP END PARALLEL DO
    
    !print *, 'xtbb', xtbb(:, thread_ix)



    do lb2 = lstart, lmax_bis
        fac = exp4_corr0/geom(lb2)*dtheta
        if(eod=='even')then
                bispectrum_temp(3, lb2) =   sum(xtbb(lb2,:))*fac
        end if 
        bispectrum_temp(4, lb2) =   sum(xteb(lb2,:))*fac
        bispectrum_temp(5, lb2) =   sum(xtee(lb2,:))*fac
    
    end do       

    deallocate(xtbb)
    deallocate(xteb)
    deallocate(xtee)

    print *, 'integration completed!'

    if(eod/='even')then
            thread_ix = 1
            !$ thread_ix = OMP_GET_MAX_THREADS()

            
            allocate(wig3_202(0:lmax_bis+lmax_teb+5, thread_ix) )
            allocate(wig2_202(0:lmax_bis+lmax_teb+5, thread_ix) )
            allocate(wig3_022(0:lmax_bis+lmax_teb+5, thread_ix ))
            allocate(wig2_022(0:lmax_bis+lmax_teb+5, thread_ix ))
            allocate(wig3_000(0:lmax_bis+lmax_teb+5, thread_ix ))
            allocate(wig2_000(0:lmax_bis+lmax_teb+5, thread_ix ))
             wig3_202= 0
             wig2_202 = 0
             wig3_000= 0
             wig2_000 = 0
             wig3_022 = 0
             wig2_022 = 0
               
            print *, 'tbb approximation'

            !$OMP PARALLEL DO DEFAULT(PRIVATE),  &
            !$OMP SHARED( bispectrum_temp, lstart, lmax_bis, lmax_alpha, lmax_teb, eod ,exp4_corr0, geom, lmax_ea), &
            !$OMP SHARED(wig3_202, wig3_022, wig2_202, wig2_022, wig3_000, wig2_000), &
            !$OMP SHARED(cs4, ifac_123, ifac_231, ifac_321, CTa,  CEa, Caa, CEE, CTE), schedule(dynamic) 
            do lb2 = lstart, lmax_bis
                fac = exp4_corr0/geom(lb2)
                
               ! make appproximation for TBB, added in 2020/4/28
               lb1= 4
               lb3 = lb2+3
              
              
            
                ppx1=0
                ppx2= 0
                yyy_tbb = 0
                yy_tbb= 0
                do l1 = 2, lmax_teb
                        thread_ix = 1
                       !$   thread_ix = OMP_GET_THREAD_NUM()+1

                        yy_tbb = 0
                        topl = l1+lb2
                        ! 3j(l1,l2,lb2; 2,0,-2)
                        botl = max(0, abs(l1-lb2))
                        nl = topl- botl+1
                        call DRC3JJ(dble(lb2), dble(l1), dble(-2), dble(2), dble(botl), dble(topl),wig2_202(botl:, thread_ix), nl, IER)

                      topl = l1+lb3
                        ! 3j(l1,l2,lb3; 2,0,-2)
                        botl = max(0, abs(l1-lb3))
                        nl = topl- botl+1
                       call DRC3JJ(dble(lb3), dble(l1), dble(-2), dble(2), dble(botl), dble(topl), wig3_202(botl:, thread_ix), nl, IER)
                        
                       topl = l1+lb2
                        ! 3j(l1,l2,lb2; 0,-2, 2)
                        botl = max(2, abs(l1-lb2))
                        nl = topl- botl+1
                        call DRC3JJ(dble(lb2), dble(l1), dble(2), dble(0), dble(botl), dble(topl),wig2_022(botl:, thread_ix), nl, IER)

                      topl = l1+lb3
                        ! 3j(l1,l2,lb3; 0,-2, 2)
                        botl = max(2, abs(l1-lb3))
                        nl = topl- botl+1
                       call DRC3JJ(dble(lb3), dble(l1), dble(2), dble(0), dble(botl), dble(topl), wig3_022(botl:, thread_ix), nl, IER)
                        
                       topl = l1+lb2
                        ! 3j(l1,l2,lb2; 0,0,0)
                        botl = max(0, abs(l1-lb2))
                        nl = topl- botl+1
                        call DRC3JJ(dble(lb2), dble(l1), dble(0), dble(0), dble(botl), dble(topl),wig2_000(botl:, thread_ix), nl, IER)

                      topl = l1+lb3
                        ! 3j(l1,l2,lb3; 0,0,0)
                        botl = max(0, abs(l1-lb3))
                        nl = topl- botl+1
                       call DRC3JJ(dble(lb3), dble(l1), dble(0), dble(0), dble(botl), dble(topl), wig3_000(botl:, thread_ix), nl, IER)
                               

                         px1 = 0 
                         px2=0       
                         y_tbb=0
                         l2start = max(2,abs(lb3-l1), abs(lb2-l1) )  
                         l2end = min(lmax_ea, lb3+l1, lb2+l1) 
                         
                         do l2 =l2start , l2end
                                px1= px1 + CEa(l2)*(cs4 + (1-2*mod(l1+l2+lb3, 2)) )*wig3_202(l2, thread_ix)*wig3_022(l2, thread_ix)*(2*l2+1.)/4./pi
                                px2= px2 + CEa(l2)*(cs4 + (1-2*mod(l1+l2+lb2, 2)) )*wig2_202(l2, thread_ix)*wig2_022(l2, thread_ix)*(2*l2+1.)/4./pi
                         
                                if(l2 .le. lmax_alpha)then
                                        y_tbb(1) = y_tbb(1) +  Caa(l2)*wig3_000(l2, thread_ix)*wig3_202(l2, thread_ix)*(2*l2+1)/4./pi
                                        y_tbb(2) = y_tbb(2) + Caa(l2)*(1- (1-2*mod(l1+l2+lb3, 2))*cs4 )*wig3_202(l2, thread_ix)**2*(2*l2+1)/4./pi
                                        y_tbb(3) = y_tbb(3)+ Caa(l2)*wig2_000(l2, thread_ix)*wig2_202(l2, thread_ix)*(2*l2+1)/4./pi
                                        y_tbb(4) = y_tbb(4) + Caa(l2)*(1- (1-2*mod(l1+l2+lb2, 2))*cs4 )*wig2_202(l2, thread_ix)**2*(2*l2+1)/4./pi

                                end if 
                         end do 
                         
                         if(l1.le. lmax_ea)then
                                 ppx1 = ppx1 +  px1*CEa(l1)*(2*l1+1)       
                                 ppx2 = ppx2 +  px2*CEa(l1)*(2*l1+1)       
                         end if 
                        
                         yy_tbb(1) = yy_tbb(1) +  y_tbb(1)*CEa(l1)*(2*l1+1)
                         yy_tbb(2) = yy_tbb(2)  +  y_tbb(2)*CEE(l1)*(2*l1+1)
                         yy_tbb(3) = yy_tbb(3) +  y_tbb(3)*CEa(l1)*(2*l1+1)
                         yy_tbb(4) = yy_tbb(4) +  y_tbb(4)*CEE(l1)*(2*l1+1)
              end do 

              pppx1  = ppx1*4*ifac_123(lb2)*CTa(lb1)
              pppx2  = -ppx2*4*ifac_123(lb2)*CTa(lb1)

              yyy_tbb(1)  = -yy_tbb(1)*8*(cs4+1)*ifac_321(lb2)*CTE(lb1)
              yyy_tbb(2)  = -yy_tbb(2)*4*ifac_123(lb2)*CTa(lb1)
              yyy_tbb(3)  = yy_tbb(3)*8*(cs4+1)*ifac_231(lb2)*CTE(lb1)
              yyy_tbb(4)  = yy_tbb(4)*4*ifac_123(lb2)*CTa(lb1)

              if (lb3.gt. lmax_ea)then
                        bispectrum_temp(3, lb2) =( (cs4-1)*( - ifac_123(lb2)*(CTa(lb1)*CEE(lb3)- CTa(lb1)*CEE(lb2) ) ) &
                                + pppx1+pppx2+ sum(yyy_tbb(:)) )*fac
 
              else if(lb2.gt. lmax_ea)then
                         bispectrum_temp(3, lb2) =( (cs4-1)*( ifac_321(lb2)*CTE(lb1)*CEa(lb3) &
                                - ifac_123(lb2)*(CTa(lb1)*CEE(lb3)- CTa(lb1)*CEE(lb2) ) ) &
                                + pppx1+pppx2+ sum(yyy_tbb(:)) )*fac
              else
                        bispectrum_temp(3, lb2) =( (cs4-1)*( ifac_321(lb2)*CTE(lb1)*CEa(lb3) &
                                -  ifac_231(lb2)*CTE(lb1)*CEa(lb2)&
                                - ifac_123(lb2)*(CTa(lb1)*CEE(lb3)- CTa(lb1)*CEE(lb2) ) ) &
                                + pppx1+pppx2+ sum(yyy_tbb(:)) )*fac
              end if 
            end do 
           !$OMP END PARALLEL DO
           
            deallocate(wig3_202)
            deallocate(wig2_202)
            deallocate(wig3_022)
            deallocate(wig2_022)
            deallocate(wig3_000)
            deallocate(wig2_000)
         
      end if 

   
    !call date_and_time(time = current_time)

    print *, 'temperature bispectrum ', bispectrum_temp(3, lmax_bis-50:lmax_bis)

    END SUBROUTINE temp_bispectrum
    
    !---------------------------------------------------------------------------------------------------------------------------------
