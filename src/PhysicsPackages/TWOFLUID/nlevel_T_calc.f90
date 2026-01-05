PROGRAM radexpintcalc
	IMPLICIT NONE
	double precision::sol,oldsol,diff,minT,maxT,T
	double precision::nuarr,exf,dy,yn,ymax,ymin,E0y,E1y,E2y
	double precision,parameter::cli=299792458.d0 !Speed of light in m/s
	double precision,parameter::kboltz=1.38064852e-23 !Boltzmann Constant [m^2 kg s^-2 K^-1]
	double precision,parameter::h=6.62607004e-34 !Planck's constant in m2 kg s^-1
	double precision,parameter::a0bohr=5.29e-11 !Bohr's radius [m] 
	double precision,parameter::melec=9.10938356e-31 !Electron mass [kg]
	integer,parameter::n_levels=5 !number of levels
	double precision,parameter::pi=3.14159265359 !pi
    double precision::Eion(6) !Energy to ionise
    double precision::rn(6),bn(6),garr(6),An0,Bn0,dt,E0z,E1z,E2z,yhat,zhat
    double precision::zizn,ziyn,zn
    double precision::xrat(n_levels,n_levels),Enn(n_levels,n_levels),rnn(n_levels,n_levels),gauntfac(n_levels,n_levels)
    double precision::fnn(n_levels,n_levels),Ann(n_levels,n_levels),Bnn(n_levels,n_levels),G_T(n_levels+1,n_levels+1)
	integer::ii,jj,ti,nsamps
	
	! output data into a file 
    open(1, file = 'colexp.dat', status = 'replace')  
   
    nsamps=100000
    
    ymin=0.54d0/13.6d0*2.18e-18/kboltz/100000.d0 !Assuming a maximum electron temperature of 100000 K, and 0.54 eV ionisation energy
    !ymax=0.45d0+2.18e-18/kboltz/100.d0 !Minimum simulation temperature 100 K,
	ymax=108.0 !Values are practially zero after this (10^-50).

    minT=2.0 !log temperatures
    maxT=8.0

	write(1,*) nsamps+1, minT,maxT

	dt=(maxT-minT)/(nsamps)
	
	
!print*,Telec
    Eion=[13.6d0,3.4d0,1.51d0,0.85d0,0.54d0,0.0d0] !in eV
    Eion=Eion/13.6d0*2.18e-18 !Convert to joules (to be dimensionally correct)

    rn(1)=0.45d0 !Equation 31
    rn(2)=1.94d0*2.d0**(-1.57d0) !Equation 32
    rn(3)=1.94d0*3.d0**(-1.57d0) !Equation 32
    rn(4)=1.94d0*4.d0**(-1.57d0) !Equation 32
    rn(5)=1.94d0*5.d0**(-1.57d0) !Equation 32
    rn(6)=1.94d0*6.d0**(-1.57d0) !Equation 32

    bn(1)=-0.603d0 !Equation 25
    bn(2)=(1.0d0/2.d0)*(4.0d0-18.63d0/2.d0+36.24d0/(2.d0**2)-&
                    28.09d0/(2.d0**3)) !Equation 26
    bn(3)=(1.0d0/3.d0)*(4.0d0-18.63d0/3.d0+36.24d0/(3.d0**2)-&
                    28.09d0/(3.d0**3)) !Equation 26
    bn(4)=(1.0d0/4.d0)*(4.0d0-18.63d0/4.d0+36.24d0/(4.d0**2)-&
                    28.09d0/(4.d0**3)) !Equation 26
    bn(5)=(1.0d0/5.d0)*(4.0d0-18.63d0/5.d0+36.24d0/(5.d0**2)-&
                    28.09d0/(5.d0**3)) !Equation 26
    bn(6)=(1.0d0/6.d0)*(4.0d0-18.63d0/6.d0+36.24d0/(6.d0**2)-&
                    28.09d0/(6.d0**3)) !Equation 26
    
    do ii=1,n_levels; do jj=ii+1,n_levels
        xrat(ii,jj)=1.0d0-(dble(ii)/dble(jj))**2 !ratio of transition energy to ionisation energy of lower level
        Enn(ii,jj)=Eion(ii)-Eion(jj) !Difference in ionisation energies
        rnn(ii,jj)=rn(ii)*xrat(ii,jj) !Equation 30

        !Gaunt factors from table 1 and equation 4
        !Constants so can be moved for speed
        if (ii .eq. 1) then 
            gauntfac(ii,jj)=1.1330d0-0.4059d0/xrat(ii,jj)+0.07014d0/(xrat(ii,jj)**2)
        else if (ii .eq. 2) then 
            gauntfac(ii,jj)=1.0785d0-0.2319d0/xrat(ii,jj)+0.02947d0/(xrat(ii,jj)**2)
        else
            gauntfac(ii,jj)=0.9935d0+0.2328d0/dble(ii)-0.1296d0/(dble(ii)**2)&
                -(1.0d0/xrat(ii,jj))*(1.0d0/dble(ii))*(0.6282d0-0.5598d0/dble(ii)+0.5299d0/(dble(ii)**2))&
                +(1.0d0/xrat(ii,jj))**2*(1.0d0/dble(ii)**2)*(0.3887d0-1.181d0/dble(ii)+1.470d0/(dble(ii)**2))
        endif

        fnn(ii,jj)=32.0d0/3.0d0/dsqrt(3.d0)/pi*dble(ii)/dble(jj)**3/(xrat(ii,jj)**3)*gauntfac(ii,jj) !Equation 3

        Ann(ii,jj)=2.0d0*dble(ii)**2/xrat(ii,jj)*fnn(ii,jj) !Equation 11
        Bnn(ii,jj)=4.0d0*dble(ii)**4/(dble(jj)**3)/(xrat(ii,jj)**2)*(1.0d0+4.0d0/3.0d0/xrat(ii,jj)+bn(ii)/(xrat(ii,jj)**2)) !Equation 23
    enddo;enddo
    
	do ti=0,nsamps

        T=10.0**(minT+dt)
        
        do ii=1,n_levels
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !This loop does the excitation rate part
            do jj=ii+1,n_levels

		        !yn=phi_I(ni)/13.6d0*2.18e-18/kboltz/T 
		        
		        yhat=Enn(ii,jj)/kboltz/T !Equation 37

                zhat=rnn(ii,jj)+Enn(ii,jj)/kboltz/T   !Equation 38

                !call ionexpfittest(yhat,0.d0,0.001d0,0.0001d0,E0y)   !Equation 8
                call ionexpfittest(yhat,1.d0,0.001d0,0.0001d0,E1y)
                call ionexpfittest(yhat,2.d0,0.001d0,0.0001d0,E2y)
                
                call ionexpfittest(zhat,1.d0,0.001d0,0.0001d0,E1z)
                call ionexpfittest(zhat,2.d0,0.001d0,0.0001d0,E2z)
				
                !Temperature dependent part of rate coefficient for excitation Equation 36 using electron mass
                G_T(ii,jj)=dsqrt(8.d0*kboltz*T/pi/melec)*2.0d0*&
                        dble(ii)**2/xrat(ii,jj)*pi*a0bohr**2*yhat**2*&
                              (Ann(ii,jj)*((1.0d0/yhat+0.5d0)*E1y-(1.0d0/zhat+0.5d0)*E1z)+&
                              (Bnn(ii,jj)-Ann(ii,jj)*dlog(2.d0*dble(ii)**2/xrat(ii,jj)))*(E2y/yhat-E2z/zhat))

            enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !This is for the ionisation part
            yn=Eion(ii)/kboltz/T
            zn=rn(ii)+Eion(ii)/kboltz/T

            call ionexpfittest(yhat,0.d0,0.001d0,0.0001d0,E0y)   !Equation 8
            call ionexpfittest(yhat,1.d0,0.001d0,0.0001d0,E1y)
            call ionexpfittest(yhat,2.d0,0.001d0,0.0001d0,E2y)
            
            call ionexpfittest(zhat,0.d0,0.001d0,0.0001d0,E0z)
            call ionexpfittest(zhat,1.d0,0.001d0,0.0001d0,E1z)
            call ionexpfittest(zhat,2.d0,0.001d0,0.0001d0,E2z)
            
            ziyn=E0y-2.0d0*E1y+E2y !Equation 42
            zizn=E0z-2.0d0*E1z+E2z !Equation 42
            
            
            if (ii .eq. 1) then 
                    garr(1)= 1.1330d0
                    garr(2)=-0.4059d0
                    garr(3)= 0.07014d0
                else if (ii .eq. 2) then 
                    garr(1)= 1.0785d0
                    garr(2)=-0.2319d0
                    garr(3)= 0.02947d0
                else
                    garr(1)=0.9935d0+0.2328d0/dble(ii)-0.1296d0/(dble(ii)**2)
                    garr(2)=(-1.0d0/dble(ii))*(0.6282d0-0.5598d0/dble(ii)+0.5299d0/(dble(ii)**2))
                    garr(3)=(1.0d0/dble(ii))**2*(0.3887d0-1.181d0/dble(ii)+1.470d0/(dble(ii)**2))
            endif
                
            !Equation 20
            An0=32.0d0/3.0d0/dsqrt(3.d0)/pi*dble(ii)+garr(1)/3.0d0+garr(2)/4.0d0+garr(3)/5.0d0
            Bn0=2.0d0/3.0d0*dble(ii)**2*(5.0d0+bn(ii)) !Equation 24

            !Ionisation coefficients Equation 35 using electron mass
            G_T(ii,n_levels+1)=dsqrt(8.d0*kboltz*T/pi/melec)&
                *2.0d0*dble(ii)**2*pi*a0bohr**2*yn**2*&
                (An0*(E1y/yn-E1z/zn)+&
                (Bn0-An0*dlog(2.d0*dble(ii)**2))*(ziyn-zizn))
        enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !some writting routine
        write(1,*) minT+dt, G_T(1,2), G_T(1,3), G_T(1,4), G_T(1,5), G_T(1,6), &
                                      G_T(2,3), G_T(2,4), G_T(2,5), G_T(2,6), &
                                                G_T(3,4), G_T(3,5), G_T(3,6), &
                                                          G_T(4,5), G_T(4,6), &
                                                                    G_T(5,6)
	 enddo
	 
	 close(1) 

END PROGRAM radexpintcalc
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine ionexpfittest(zhat,istar,stepsize,tol,sol)
    !numerical integration of the exponential
    !This needs work
    double precision,intent(in)::zhat,istar,stepsize,tol
    double precision,intent(out)::sol
    double precision::a,b,fa,fb,sol0,dif,solt

    a=1.d0
    b=a+stepsize
    fa=exp(-zhat*a)*a**(-istar)
    fb=exp(-zhat*b)*b**(-istar)
    sol0=0.5d0*(b-a)*(fa+fb)
    sol=sol0

    dif=1.d0

    do while (dif .gt. tol)
        a=b
        b=a+stepsize
        fa=exp(-zhat*a)*a**(-istar)
        fb=exp(-zhat*b)*b**(-istar)
        solt=0.5d0*(b-a)*(fa+fb)
        sol=sol+solt
        dif=solt/sol0
    enddo

  end subroutine ionexpfittest
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



