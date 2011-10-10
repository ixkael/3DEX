MODULE cosmotools

    USE healpix_types
    IMPLICIT NONE
    
    REAL(DP), PARAMETER :: c = 299792458_dp  	! Speed of light
    REAL(DP), PARAMETER :: tcmb = 2726000_dp 	! CMB temperature [micro Kelvin]

    CONTAINS  
    
    ! ---------------------------------------------------------------------------------------!    
    
    
    SUBROUTINE cosmo_z2s( h_in, omega_m_in, omega_l_in, omega_b_in, wa_in, w0_in, z, nbpts_in, sk )   ! need in cosmo : h, omega_b, omega_m, omega_l
	    
	INTEGER(I4B) :: nbpts_in, nbpts
	REAL(DP) :: h, omega_m, omega_l, omega_b, omega_k, wa, w0, sqrtk
	REAL(DP) :: h_in, omega_m_in, omega_l_in, omega_b_in, wa_in, w0_in
	COMMON h, omega_m, omega_l, omega_b, omega_k, wa, w0, sqrtk, nbpts
	
	REAL(DP) :: omega_dm, omega, gamma, h0, rh, k, r0
	
	REAL(DP), DIMENSION(1:nbpts_in) :: z, chi, sk
		
	h = h_in
	omega_m = omega_m_in
	omega_l = omega_l_in
	omega_b = omega_b_in
	wa = wa_in
	w0 = w0_in
	nbpts = nbpts_in
	
	h0 = 100.0_dp * h
	rh = c/(1000*h0)
	
	omega_dm = omega_m - omega_b 	! Dark matter density
	omega = omega_m + omega_l      	! Total density
	omega_k = 1.0_dp - omega           		! Curvature
	gamma = omega_m * h * &				! Sugiyama (1995, APJS, 100, 281)
	& exp(-omega_b*(1.0_dp + sqrt(2.0_dp*h)/omega_m)) 
	
	IF( omega > 1.001 ) THEN 			! closed
	    k = 1.0_dp                         
	    sqrtk = sqrt(omega-1.0_dp)   	! useful factor
	ENDIF                       			! with  a=r/r0
	IF( omega > 0.999 .and. omega < 1.001 ) THEN
	    k = 0.0_dp                   		! flat
	    sqrtk = 1.0_dp
	ENDIF
	IF( omega < 0.999 ) THEN        	! open
	    k = -1.0_dp  
	    sqrtk = sqrt(1.0_dp - omega)       
	ENDIF
	
	r0 = rh / sqrtk  

	! comoving distance in units of R_0  [1] 
	print*,"Computing chi..."
	call chi_z(z,chi)
	
	! sk(chi): radial comoving distance in units of R_0  [1],
	IF( omega > 1.001 ) THEN 
  	    sk = chi
  	ENDIF                       			
	IF( omega > 0.999 .and. omega < 1.001 ) THEN 
	    sk = sin(chi)
	ENDIF
	IF( omega < 0.999 ) THEN 
  	    sk = sinh(chi)
  	ENDIF
  	
  	sk = abs(sk) * r0
  	
  	RETURN

    END SUBROUTINE cosmo_z2s
    
    
    ! ---------------------------------------------------------------------------------------!
    
    
    SUBROUTINE chi_z(zin, chi)

	INTEGER(I4B) :: nbpts
	REAL(DP) :: h, omega_m, omega_l, omega_b, omega_k, wa, w0, sqrtk
	COMMON h, omega_m, omega_l, omega_b, omega_k, wa, w0, sqrtk, nbpts
	
    	INTEGER(I4B) :: i
    	REAL(DP), DIMENSION(1:nbpts) :: zin,chi
    	
    	!real, external :: chi_a_int
    	!real, external :: fz
    	
    	real, parameter :: a = 0.0E+00
    	real abserr
    	real b
    	real, parameter :: epsabs = 0.0E+00
    	real, parameter :: epsrel = 0.001E+00
    	integer ier
    	integer, parameter :: key = 6
    	integer neval
    	real result
    
    	DO i = 1, nbpts
    	
    	    b = zin(i)

       	    call qag ( fz, a, b, epsabs, epsrel, key, result, abserr, neval, ier )
       	           	    	
       	    chi(i) = result
    
        ENDDO
        	    
    END SUBROUTINE chi_z
    
    	    
    ! ---------------------------------------------------------------------------------------!
    
    	
    REAL FUNCTION fz(z)

	INTEGER(I4B) :: nbpts
	REAL(DP) :: h, omega_m, omega_l, omega_b, omega_k, wa, w0, sqrtk
	COMMON h, omega_m, omega_l, omega_b, omega_k, wa, w0, sqrtk, nbpts
	
	REAL(DP) :: z
		
	fz = ( omega_m*((1+z)**3) + omega_k*((1+z)**2) + omega_l )**(-0.5_dp) 
	
	RETURN

    END FUNCTION fz
    	    
    
    ! ---------------------------------------------------------------------------------------!
    
    
END MODULE cosmotools
