pro clkn_diff_contour, clkn_list1, clkn_list2, CARR=clkn_arr, LARR=l_arr, KARR=kln_arr, NLMAX=nlmax, NNMAX=nnmax, NLEVELS=nlevels, NPOINTS=npoints, KLOG=klog, ELIM=elim, NSIDE=nside,RMAX=rmax
;Written by Boris Leistedt, October 2011
;Routine used to output power spectrum error plot,
; as in Leistedt, Rassat, Refregier and Stark (2011).

if keyword_set(NLMAX) then lmax=nlmax else lmax=min( [max(clkn_list1(*,0)),max(clkn_list2(*,0))] )
if keyword_set(NLMAX) then nmax=nnmax else nmax=min( [max(clkn_list1(*,1)),max(clkn_list2(*,1))] )

clkn_arr = dblarr(lmax+1,nmax)
l_arr = dblarr(lmax+1,nmax)
kln_arr = dblarr(lmax+1,nmax)

if keyword_set(NPOINTS) then np=npoints else np=nmax
if keyword_set(ELIM) then lim=elim else lim=0.001
if keyword_set(NLEVELS) then nlev=nlevels else nlev=100
 
for l=0, lmax do begin
    indl1 = where(clkn_list1(*,0) eq l) 
    indl2 = where(clkn_list2(*,0) eq l) 
    for n=1, nmax do begin
        cl1 = clkn_list1(indl1,*)
        cl2 = clkn_list2(indl2,*)
        indn1 = where( cl1(*,1) eq n)
        indn2 = where( cl2(*,1) eq n)
        v1 = cl1(indn1,3) ;clkn1_sm(n-1)
        v2 = cl2(indn2,3) ;clkn2_sm(n-1)
        l_arr(l,n-1) = l
        if keyword_set(KLOG) then kln_arr(l,n-1) = clkn_list1(indn1,2) else kln_arr(l,n-1) = n
        if ( clkn_list1(indn1,2) - clkn_list2(indn2,2) < 0.000000001 ) then begin
            print, "PROBLEM : kln"
            print,n, l
            print, clkn_list1(indn1,2),clkn_list2(indn2,2)
         endif
        if( (v1 gt 0) and (v2 gt 0) ) then begin
            clkn_arr(l,n-1) = (v2 / v1) - 1
            if( abs(clkn_arr(l,n-1)) gt lim ) then begin
               print, " "
               print, l, n
               print, v1, v2
               print, clkn_arr(l,n-1)           ;clkn_arr(l,n-1) = 1e-16 ; 
            endif  
            if( clkn_arr(l,n-1) eq 0 ) then clkn_arr(l,n-1) = 1e-16
        endif else begin
            clkn_arr(l,n-1) = 1e-16
            print, "PROBLEM : clkn"
            print, l, n
            print, kln_arr(l,n-1)
            print, v1, v2, abs( v2 / v1 - 1 )  
            print," "
        endelse
    endfor
 endfor



cmin = -lim ;max( [ min(clkn_arr), -lim ] )
cmax = lim ;min( [ max(clkn_arr), lim ] )
clim = lim ; max( [cmax, abs(cmin)] )
levels = [ xgen(-clim, clim, NPOINTS=nlevel) , clim ] ;, LOGPLOT=2)
print, "Levels between ", min(levels), max(levels)


indic='n'
if keyword_set(KLOG) then indic='k'

metatitle='Relative error for C(l,'+STRCOMPRESS(indic)+') for nside='+STRCOMPRESS(nside, /REMOVE_ALL)+' and Rmax='+STRCOMPRESS(rmax, /REMOVE_ALL)

filename='err_'+indic+'_01_'+string(nside)+'_'+string(rmax)+'.ps'
filename=STRCOMPRESS(filename, /REMOVE_ALL)

PS_Start, FILENAME=filename, FONT=0

if keyword_set(KLOG) then begin
   contour, clkn_arr, l_arr, kln_arr, /ylo, xrange=[0,lmax],xtickname=[' ','10', '20', '30', '40'],/xstyle, yrange=[0.001, 0.2], /ystyle, levels=levels, c_labels=levels, color=1, background=0, title=metatitle,charsize=1.55,xcharsize=1.3,ycharsize=1.15,ytitle='k',xtitle='l',/fill
endif else begin
   contour, clkn_arr, l_arr, kln_arr, xrange=[0,lmax],xtickname=[' ','10', '20', '30', '40'],/xstyle, yrange=[1, nmax],/ystyle, levels=levels, c_labels=levels, color=1, background=0, title=metatitle,charsize=1.55,xcharsize=1.3,ycharsize=1.3,ytitle='n',xtitle='l',/fill
endelse

colorbar, ncolors=15, divisions=2,range=[-clim,clim],format='(F6.3)',color=1,bottom=2, charsize=1.3,position=[0.07, 0.136, 0.11, 0.29];[0.98, 0.031, 0.992, 0.15]

PS_End, /PNG

;plot,kln_arr(0,*),clkn_arr(0,*), /xlo, xrange=[min(kln_arr),max(kln_arr)],/xstyle, /ylo, yrange=[min(clkn_arr), max(clkn_arr)],/ystyle, charsize=3
;print,clkn_arr(0,*),kln_arr(0,*)
;for l=1, nlmax do begin
   ;print,kln_arr(l,*),clkn_arr(l,*)
;   oplot,kln_arr(l,*),clkn_arr(l,*)
;end

end
