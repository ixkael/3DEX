#-------------------------------------------------------------
# Bourne/Korn shell functions required to configure 3DEX packages
# ------------------------------------------------------------
#
#=====================================
#=========== General usage ===========
#=======================3dex==============
#   checkDir: search for installation directories and create them
#             is necessary
#   echoLn:
#   findFITSLib: search for FITSIO library
#   findFFTW: search for FFTW library
#   fullPath: convert relative to absolute directory names
#
#
#-------------
checkDir () {
    l=""
    for d in $*; do
	[ ! -d $d ] && l="$l $d"
    done
    if [ "x$l" != "x" ]; then
	echo "Warning: The following directories could not be found:"
	for d in $l; do
	    echo "$d"
	done
	echoLn "Should I attempt to create these directories (Y|n)? "
	read answer
	if [ "x$answer" != "xn"  -a  "x$answer" != "xN"  ]; then
	    for d in $l; do
		mkdir $d 1>${DEVNULL} 2>&1
		if [ $? -gt 0 ]; then
		    echo "Error: Could not create directory $d"
		    crashAndBurn
		fi
	    done
	else
	    echo "Create installation directories first."
	    crashAndBurn
	fi
    fi
}    
#-------------
echoLn () {
#     if [ "${OS}" = "Linux" -o "${OS}" = "Darwin" ]; then
# 	echo -n "$*"
#     else
# 	echo "$*\c"
#     fi
    ${PRINTF} "$*"
}
#-------------
findFITSLib () {
    for dir in $1 /usr/lib /usr/lib64 /usr/local/lib /usr/local/lib64 /usr/local/lib/cfitsio /usr/local/lib64/cftisio /usr/local/src/cfitsio ${HOME}/lib ${HOME}/lib64 ./src/cxx/${DMPX_TARGET}/lib/ ; do
	if [ -r "${dir}/lib${LIBFITS}.a" ] ; then
	    FITSDIR=$dir
	    break
	fi	    
    done
}
#-------------
findFITSInclude () {
    for dir in $* /usr/include /usr/local/include /usr/local/src/cfitsio ${HOME}/include ${HOME}/include64 ./src/cxx/${DMPX_TARGET}/include/ ; do
	if [ -r "${dir}/fitsio.h" ] ; then
	    FITSINC=$dir
	    break
	fi
    done
}

#-------------
fullPath () {
    t='TEMP=`cd $TEMP; pwd`'
    for d in $*; do
	eval `echo $t | sed 's/TEMP/'$d'/g'`
    done
}
#-------------
askMake () {
if [ "${MAKESET}" = 0 ] ; then
    echoLn "Enter make command ($MAKE): "
    read answer
    [ "x$answer" != "x" ] && MAKE=$answer    
    MAKESET=1
fi
}
#-------------
goodBye () {
    echo     
    if [ -s Makefile -a ${edited_makefile} -eq 1 ] ; then
	echo
	echo "You can run \"(GNU)make\" to build all the packages configured so far,"
	echo
    fi
    echo "Good Bye !"
    echo
    exit 0
}
#-------------
crashAndBurn () {
    echo
    echo "Something went wrong ..."
    echo "Quitting configuration script !"
    echo
    exit -1
}
#=====================================
#=========== F90 pakage ===========
#=====================================
#
#   setF90Defaults: set default values of variables
#   sun_modules : test weither the Sun compiler creates modules ending with .M or .mod
#   ifc_modules : test weither the IFC compiler creates .d or .mod (version7) modules
#   GuessCompiler: tries to guess compiler from operating system
#   askFFT: ask user for his choice of fft, find fftw library
#   askOpenMP: ask user for compilation of OpenMP source files
#   countUnderScore: match trailing underscores with fftw
#   IdentifyCompiler : identify Non native f90 compiler
#   add64bitF90Flags: add 64 bit flags to F90 (and C) compiler
#   askUserF90:  ask user the f90 compiler command
#   askUserMisc:  ask user to confirm or change various defaults
#   editF90Makefile: create makefile from template
#   makeProfile: create profile
#   installProfile: modify user's shell profile if agreed
#   showDefaultDirs: show default directories
#   updateDirs: update those directories
#   showActualDirs: show actual directories
#
#-------------
setF90Defaults () {
    FC="f90"
    #HEALPIX=""
    CC="cc"
    FFLAGS="-I\$(F90_INCDIR)"
    CFLAGS="-O"
    LDFLAGS="-L\$(F90_LIBDIR) -l3dex"
    PIXFLAGS="-I\$(HEALPIX)/include -L\$(HEALPIX)/lib -L\$(FITSDIR) -lhealpix -lhpxgif -l\$(LIBFITS)"
    F90_BINDIR="./bin"
    F90_INCDIR="./include"
    F90_LIBDIR="./lib"
    F90_OUTDIR="./out"
    F90_PIXFLAGS=""
    DIRSUFF=""
    MOD="mod"
    #OS=`uname -s`
    FFTSRC="healpix_fft"
    FFTLD=" "
    ADDUS=" "
    LIBFFTW="dfftw"
    FPP="-D"
    PARALL=""
    PRFLAGS=""
    AR="ar rv"
    FTYPE=""
    PPFLAGS=""
    FF64=""
    CF64=""
    AR64=""
    PGLIBSDEF="-L/usr/local/pgplot -lpgplot -L/usr/X11R6/lib -lX11"
    WLRPATH="" # to add a directory to the (linker) runtime library search path

    echo "you seem to be running $OS"

    case $OS in
	AIX)
	    FC="xlf90_r";;
	Linux)
	    FC="";;
	Darwin)
	    FC="";;
	SUPER-UX)
	    FC="f90";;
    esac

    FCNAME="$OS Native compiler"
}

# -----------------------------------------------------------------

sun_modules () {
tmpfile=to_be_removed
${CAT} > ${tmpfile}.f90 << EOF
   module ${tmpfile}
       integer :: i
   end module ${tmpfile}
EOF
   $FC -c ${tmpfile}.f90 -o ${tmpfile}.o

   if test -s ${tmpfile}.M  ; then
       MOD="M"
   else
       MOD="mod"
   fi

   ${RM} ${tmpfile}.*
}

# -----------------------------------------------------------------
ifc_modules () {
tmpfile=to_be_removed
${CAT} > ${tmpfile}.f90 << EOF
   module ${tmpfile}
       integer :: i
   end module ${tmpfile}
EOF
   $FC -c ${tmpfile}.f90 -o ${tmpfile}.o 2> ${DEVNULL}

   if test -s ${tmpfile}.d  ; then
    # version 5 and 6 of ifc
	echo "This version of ifc is no longer supported"
	echo "use a more recent version (7 or higher)"
	crashAndBurn
#         IFCMOD="d"
#         IFCINC="-cl,\$(DMPX)/include/list.pcl"
#         IFCVERSION="ifcold"
   else
       IFCMOD="mod"
#         IFCINC="-I\$(DMPX)/include"
       IFCINC="-I\$(F90_INCDIR)"
#         IFCVERSION="ifnew"
       IFCVERSION="ifcnew"
   fi

   ${RM}  ${tmpfile}.* 
   ${RM}  TO_BE_REMOVED.*
}

#----------
checkF90Fitsio () {
    cfitsiolib=$1
    sanity=`${NM} ${cfitsiolib} 2> ${DEVNULL} | ${GREP} read | ${GREP} T | ${WC} -l` # make sure that nm, grep and wc are understood
    if [ $sanity -gt 0 ] ; then
	check=`${NM} ${cfitsiolib} 2> ${DEVNULL} | ${GREP} ftgkey | ${GREP} T | ${WC} -l` # count ftgkey definition
	if [ $check -eq 0 ] ; then
	    echo 
	    echo "*WARNING*: the cfitsio library ${cfitsiolib}"
	    echo "does not seem to include the Fortran interface;"
	    echo "this will prevent compilation of the 3DEX F90 package."
	    echo 
	    echo "When installing cfitsio, make sure that a Fortran compiler is known"
	    echo "to the cfitsio configure script."
	    echo 
	fi
    fi

}
# ----------------
checkF90FitsioLink () {
# check that F90 routines can link with F90-fitsio wrappers
# requires compilation of F90 code
    tmpfile=to_be_removed
    # write simple program to link with fitsio
cat > ${tmpfile}.f90 << EOF
    program needs_fitsio
	character(len=6) :: string='abcdef'
	call ftupch(string)
    end program needs_fitsio
EOF
    # compile and link
    ${FC} ${FFLAGS}  ${tmpfile}.f90 -o ${tmpfile}.x -L${FITSDIR} -l${LIBFITS} #${WLRPATH}

    # test
    if [ ! -s ${tmpfile}.x ]; then
	echo
	echo "F90 codes do not link correctly with ${FITSDIR}/lib${LIBFITS}.a"
	echo "Check that in the cfitsio library:"
	echo " - the Fortran wrappers were correctly compiled, and"
	echo " - the library (C routines and F90 wrappers) was compiled "
	echo "   with a number of bits compatible with ${FC} ${FFLAGS}"
	crashAndBurn
    fi

    # clean up
    ${RM} ${tmpfile}.*
    

}
 
# -----------------------------------------------------------------

GuessCompiler () {
    case $OS in
	AIX)
	    IdentifyCompiler;;
# 	    FTYPE="xlf"
#	    FFLAGS="$FFLAGS -qsuffix=f=f90:cpp=F90"
#	    OFLAGS="-O"
#	    CFLAGS="$CFLAGS $FPP""RS6000"
#	    FPP="-WF,-D"
#	    PRFLAGS="-qsmp=omp"
#	    AR="ar -rsv" # archive with index table
#	    FF64="-q64"
#	    CF64="-q64"
#	    AR64="-X64";;
	SunOS)
	    sun_modules
	    FFLAGS=`echo $FFLAGS | ${SED} "s/-I/-M/g"`
	    LDFLAGS="$LDFLAGS -lm -lnsl -lsocket"
	    OFLAGS="-fast";;
	IRIX*)
	    OS="IRIX"
	    LDFLAGS="$LDFLAGS -lm"
	    OFLAGS="-O"
	    PRFLAGS="-mp";;
	Linux)
	    AR="ar -rsv" # archive with index table
  	    OFLAGS="-O"
	    IdentifyCompiler;;
	Darwin)
  	    OFLAGS="-O"
	    AR="libtool -static -s -o"  # archive with index table
	    IdentifyCompiler;;
	OSF*)
	    OS="OSF"
	    OFLAGS="-O5 -fast"
	    PRFLAGS="-omp";;
	SUPER-UX)
	    FC="f90"
	    FFLAGS="$FFLAGS"
	    OFLAGS="-C vopt"
	    CFLAGS="-C vopt"
	    FPP="-D"
	    PRFLAGS="-P openmp";;
	CYGWIN*)
	    OFLAGS="-O"
	    IdentifyCompiler;;
	*)
	    echo "\"$OS\" is not supported yet"
	    crashAndBurn;;
    esac
}

# -----------------------------------------------------------------

askOpenMP () {
    OpenMP="0"
    echo "Do you want to use :"
    echo " 0) the standard serial implementation ?"
    echo " 1) the parallel implementation (slightly slower in single CPU usage with some compilers)"
    echoLn "Enter choice                                      ($OpenMP): "
    read answer
    [ "x$answer" != "x" ] && OpenMP="$answer"
    if [ $OpenMP = 1 ] ; then
	if [ "x$PRFLAGS" != "x" ] ; then
	    # update FFLAGS
	    FFLAGS="$FFLAGS $PRFLAGS"
##	    PARALL="_omp" # no need for a different source file
	else
	    echo "3DEX+OpenMP not tested for  \"$FCNAME\" under \"$OS\" "
	    echo "Contact healpix@jpl.nasa.gov if you already used OpenMP in this configuration."
	    echo "Will perform serial implementation instead"
	    #crashAndBurn
	fi 
    fi
}

# -----------------------------------------------------------------

countUnderScore () {
tmpfile=to_be_removed
${CAT} > ${tmpfile}.f90 << EOF
    subroutine sub1()
      return
    end subroutine sub1
EOF
 case $FTYPE in
  xlf)
    $FC -qsuffix=f=f90 -c ${tmpfile}.f90 -o ${tmpfile}.o  2>&1 ${DEVNULL} ;;
  *)
    $FC -c ${tmpfile}.f90 -o ${tmpfile}.o  2>&1 ${DEVNULL} ;;
 esac

    stwo=`${NM} ${tmpfile}.o | ${GREP} sub1__ | ${WC} -l`
    sone=`${NM} ${tmpfile}.o | ${GREP} sub1_  | ${WC} -l`
    ltwo=`${NM} $lib | ${GREP} fftw_f77_one__ | ${WC} -l`
    lone=`${NM} $lib | ${GREP} fftw_f77_one_  | ${WC} -l`

    if [ $ltwo != 0 ] ; then
      if [ $stwo != 0 ] ; then
        ADDUS="$FPP""ADD0US"
      elif [ $sone != 0 ] ; then
        ADDUS="$FPP""ADD1US"
      else
        ADDUS="$FPP""ADD2US"
      fi
    elif [ $lone != 0 ] ; then
      if [ $stwo != 0 ] ; then
        echo "uncompatible trailing underscores"
        crashAndBurn
      elif [ $sone != 0 ] ; then
        ADDUS="$FPP""ADD0US"
      else
        ADDUS="$FPP""ADD1US"
      fi
    else
      if [ $stwo != 0 ] ; then
        echo "uncompatible trailing underscores"
	crashAndBurn
      elif [ $sone != 0 ] ; then
        echo "uncompatible trailing underscores"
	crashAndBurn
      else
        ADDUS="$FPP""ADD0US"
      fi
    fi

#      echo $ADDUS
   ${RM}  ${tmpfile}.*


}
# -----------------------------------------------------------------

IdentifyCompiler () {
# For Linux and Darwin
# Lahey and Fujitsu still have to be tested
        nima=`$FC -V 2>&1 | ${GREP} -i imagine1 | ${WC} -l`
        nnag=`$FC -V 2>&1 | ${GREP} -i nagware  | ${WC} -l`
        nifc=`$FC -V 2>&1 | ${GREP} -i intel    | ${WC} -l`
        npgf=`$FC -V 2>&1 | ${GREP} -i portland | ${WC} -l`
	nlah=`$FC --version 2>&1 | ${GREP} -i lahey | ${WC} -l`
	nfuj=`$FC -V 2>&1 | ${GREP} -i fujitsu | ${WC} -l`
#	nvas=`$FC | ${GREP} -i sierra | ${WC} -l`
#  	nxlf=`man $FC | ${HEAD} -10 | ${GREP} XL | ${WC} -l`
	nxlf=`$FC --help 2>&1 | ${HEAD} -15 | ${GREP} XL | ${WC} -l`
#	nxlf=`$FC -qversion 2>&1 | ${GREP} XL | ${WC} -l` # to be tested
	nabs=`$FC -V 2>&1 | ${GREP} 'Pro Fortran' | ${WC} -l`
	ng95=`$FC -dumpversion 2>&1 | ${GREP} 'g95' | ${WC} -l`
#	ngfortran=`$FC -dumpversion 2>&1 | ${GREP} 'GNU Fortran' | ${GREP} 'GCC' | ${WC} -l` # corrected 2008-11-17
#	ngfortran=`$FC -dumpversion 2>&1 | ${GREP} 'GNU Fortran' | ${WC} -l` # corrected 2009-10-12
	ngfortran=`$FC --version 2>&1 | ${GREP} 'GNU Fortran' | ${WC} -l`
	npath=`$FC -v 2>&1 | ${GREP} -i ekopath | ${WC} -l`
        if [ $nima != 0 ] ; then
                FCNAME="Imagine F compiler"
                FFLAGS="$FFLAGS -DNAG -w -dusty -mismatch_all"
		echo "$FCNAME is not supported yet"
		crashAndBurn
        elif [ $nnag != 0 ] ; then
                FCNAME="NAGWare compiler"
		PPFLAGS="-fpp"
# very sloppy compiler flags: no longer needed
#                FFLAGS="$FFLAGS -DNAG -w -dusty -mismatch_all"
# compiler flags for very thorough checking. use for debugging
#                FFLAGS="$FFLAGS -DNAG -strict95 -g -gline -C=all -u -colour"
# standard flags
                FFLAGS="$FFLAGS -DNAG -strict95"
		FI8FLAG="-double" # change default INTEGER and FLOAT to 64 bits
        elif [ $nifc != 0 ] ; then 
		ifc_modules
                FCNAME="Intel Fortran Compiler"
#                 FFLAGS="$IFCINC -Vaxlib -cm -w -vec_report0" # June 2007
                FFLAGS="$IFCINC -cm -w -vec_report0 -sox"
		MOD="$IFCMOD"
		FTYPE="$IFCVERSION"
#  		OFLAGS="-O3 -tpp7 -xW" # pentium 4
  		OFLAGS="-O3"
#		OFLAGS="-O3 -axiMKW" # generates optimized code for each Pentium platform
# 		PRFLAGS="-openmp" # Open MP enabled
		PRFLAGS="-openmp -openmp_report0" # Open MP enabled # June 2007
		FI8FLAG="-i8" # change default INTEGER to 64 bits
##		FI8FLAG="-integer-size 64" # change default INTEGER to 64 bits
		[ $OS = "Linux" ] && WLRPATH="-Wl,-R"
        elif [ $npgf != 0 ] ; then
                FCNAME="Portland Group Compiler"
		PRFLAGS="-mp" # Open MP enabled, to be tested
        elif [ $nlah != 0 ] ; then
                FCNAME="Lahey/Fujitsu Compiler"
#  		FFLAGS="$FFLAGS --nap --nchk --npca --ntrace --tpp --trap dio"
		FFLAGS="$FFLAGS --nap --nchk --npca --ntrace --tpp --trap" # (on trial version)
        elif [ $nfuj != 0 ] ; then
                FCNAME="Fujitsu Compiler"
  		FFLAGS="$FFLAGS -Am -X9 -static"
	elif [ $nxlf != 0 ] ; then
	    FTYPE="xlf"
	    if [ "$OS" = "AIX" ] ; then
		FC="xlf90_r"
		FCNAME="IBM XL Fortran"
		FFLAGS="$FFLAGS -qsuffix=f=f90:cpp=F90"
		OFLAGS="-O"
		#####CC="gcc"
		CFLAGS="$CFLAGS -DRS6000" #### 
		FPP="-WF,-D"
		PRFLAGS="-qsmp=omp" # Open MP enabled
		AR="ar -rsv" # archive with index table
		FF64="-q64"
		CF64="-q64"
		AR64="-X64"
	    else
		FC="xlf90"
		FCNAME="IBM XL Fortran for Mac OS"
		FFLAGS="$FFLAGS -qfree=f90 -qsuffix=f=f90:cpp=F90"
		OFLAGS="-O"
		CC="gcc"
		CFLAGS="$CFLAGS -DRS6000" #### 
		#### FPP="-WF,-D"
		PRFLAGS="-qsmp=omp" # Open MP enabled
	    fi
	    FI8FLAG="-qintsize=8" # change default INTEGER to 64 bits
	elif [ $nabs != 0 ] ; then
	        FCNAME="Absoft Pro Compiler"
		FFLAGS=`echo $FFLAGS | ${SED} "s/-I/-p/g"`
		FFLAGS="$FFLAGS -YEXT_NAMES=LCS -YEXT_SFX=_ -q"
		OFLAGS="-O3 -cpu:host"
		LDFLAGS="$LDFLAGS -lU77"
		CC="gcc"
	elif [ $ng95 != 0 ] ; then
	        FCNAME="g95 compiler"
		FFLAGS="$FFLAGS -DGFORTRAN"
		OFLAGS="-O3"
		CC="gcc"
		FI8FLAG="-i8" # change default INTEGER to 64 bits
	elif [ $ngfortran != 0 ] ; then
	        FCNAME="gfortran compiler"
		FFLAGS="$FFLAGS -DGFORTRAN -fno-second-underscore"
		OFLAGS="-O3"
		PRFLAGS="-fopenmp" # Open MP enabled
		CC="gcc"
		FI8FLAG="-fdefault-integer-8" # change default INTEGER to 64 bits
		[ $OS = "Linux" ] && WLRPATH="-Wl,-R"
	elif [ $npath != 0 ] ; then
	        FCNAME="PathScale EKOPath compiler"
		FFLAGS="$FFLAGS"
		OFLAGS="-O"
		CC="pathcc"	    
		PRFLAGS="-mp" # Open MP enabled
		FI8FLAG="-i8" # change default INTEGER to 64 bits
		#FI8FLAG="-default64" # change default INTEGER and FLOAT to 64 bits
        else
	    nvas=`$FC | ${GREP} -i sierra | ${WC} -l`
            if [ $nvas != 0 ] ; then
                FCNAME="Pacific/Sierra Compiler"
		echo "$FCNAME is not supported"
		crashAndBurn
	    else
                echo "$FC: Unknown compiler"
                crashAndBurn
	    fi
        fi
}

# -----------------------------------------------------------------
add64bitF90Flags () {

    if [ "x$FF64$CF64$AR64" != "x" ]; then
	echo "Do you want to make a 64 bit compilation ? [y/N]"
	read answer
	if [ "x$answer" = "xy" -o "x$answer" = "xY" ]; then
	    FFLAGS="$FFLAGS $FF64"
	    CFLAGS="$CFLAGS $CF64"
	    AR="$AR $AR64"
	fi
    fi
}
# -----------------------------------------------------------------
countF90Bits () {
    tmpfile=to_be_removed
${CAT} > ${tmpfile}.f90 <<EOF
program test
end program test
EOF

     $FC $FFLAGS ${tmpfile}.f90 -o ${tmpfile} 1>${DEVNULL} 2>&1
     f90_64=`${FILE} ${tmpfile} | ${GREP} 64 | ${WC} -l`
     ${RM}  ${tmpfile}*
}
# -----------------------------------------------------------------
countCBits () {
    tmpfile=to_be_removed
${CAT} > ${tmpfile}.c <<EOF
int main(){
}
EOF

    $CC $CFLAGS ${tmpfile}.c -o ${tmpfile} 1>${DEVNULL} 2>&1
    c_64=`${FILE} ${tmpfile} | ${GREP} 64 | ${WC} -l`
    ${RM}  ${tmpfile}*
}
# -----------------------------------------------------------------
checkF90Compilation () {
    # check that F90 compiler actually work
    # requires compilation and execution of F90 code
    tmpfile=./to_be_removed
${CAT} > ${tmpfile}.f90 <<EOF
program test
    print*,'hello'
end program test
EOF
    canrun=0
    cancompile=0
    $FC $FFLAGS ${tmpfile}.f90 -o ${tmpfile}  1>${DEVNULL} 2>&1
    [ -s ${tmpfile} ] && cancompile=1
    if [ -x ${tmpfile} ] ; then
	canrun=`${tmpfile} | grep hello | ${WC} -l`
    fi
    ${RM} ${tmpfile}*

    if [ $cancompile -eq 0 ]; then
	echo
	echo "  ERROR: Compilation with "
	echo "${FC} ${FFLAGS}"
	echo "currently fails."
	echo
	echo "Please check that this compiler is supported by your system"
	crashAndBurn
    fi
    if [ $canrun -eq 0 ]; then
	echo
	echo "  WARNING: Currently the codes compiled with "
	echo "${FC} ${FFLAGS}"
	echo "can not be executed."
	echo "Most likely, some compiler related dynamic libraries are not found."
	echo "(Check the LD_LIBRAY_PATH variable and read the compiler documentation.)"
	echo 
	echo "That can affect negatively the result of this configuration script."
	echo
	echo
    fi
}
# -----------------------------------------------------------------
checkF90LongLong () {
    # check that F90 support 8 byte integers
    # requires compilation and execution of F90 code
    tmpfile=./to_be_removed


    # conduct real test
${CAT} > ${tmpfile}.f90 <<EOF
program test
    if (selected_int_kind(16) > selected_int_kind(9)) print*,'OK'
end program test
EOF
    longlong=0
    $FC $FFLAGS ${tmpfile}.f90 -o ${tmpfile}  1>${DEVNULL} 2>&1
    if [ -x ${tmpfile} ] ; then
	longlong=`${tmpfile} | grep OK | ${WC} -l`
    fi
    ${RM} ${tmpfile}*

}
# -----------------------------------------------------------------

askUserF90 () {
    echoLn "enter name of your F90 compiler ($FC): "
    read answer
    [ "x$answer" != "x" ] && FC="$answer"
    echoLn "Healpix library ---"
    echoLn "Enter Healpix directory ($HEALPIX): "
    read answer
    [ "x$answer" != "x" ] && HEALPIX="$answer"
}

# -----------------------------------------------------------------

showDefaultDirs () {
    echo " compiled 3DEX products will be:"
    echo "F90_BINDIR =  ${F90_BINDIR}[suffix]"
    echo "F90_INCDIR =  ${F90_INCDIR}[suffix]"
    echo "F90_LIBDIR =  ${F90_LIBDIR}[suffix]"
    echo "F90_OUTDIR =  ${F90_OUTDIR}[suffix]"
#    echo " and the Makefile will be copied into Makefile[suffix]"
}

updateDirs () {
    F90_BINDIR=${F90_BINDIR}$DIRSUFF
    F90_INCDIR=${F90_INCDIR}$DIRSUFF
    F90_LIBDIR=${F90_LIBDIR}$DIRSUFF
    F90_OUTDIR=${F90_OUTDIR}$DIRSUFF
}

showActualDirs () {
    echo " compiled 3DEX products will be:"
    echo "F90_BINDIR =  ${F90_BINDIR}"
    echo "F90_INCDIR =  ${F90_INCDIR}"
    echo "F90_LIBDIR =  ${F90_LIBDIR}"
    echo "F90_OUTDIR =  ${F90_OUTDIR}"
}
# -----------------------------------------------------------------
askUserMisc () {
    echo "  Note: your Fortran compiler is $FCNAME"

    echoLn " "

    add64bitF90Flags

    showDefaultDirs
    echoLn "enter suffix for directories ($DIRSUFF): "
    read answer
    [ "x$answer" != "x" ] && DIRSUFF="$answer"
    updateDirs
    showActualDirs

    checkDir $F90_BINDIR $F90_INCDIR $F90_LIBDIR $F90_OUTDIR
    fullPath F90_BINDIR F90_INCDIR F90_LIBDIR F90_OUTDIR

    echoLn " "

    echoLn "enter compilation flags for $FC compiler ($FFLAGS): "
    read answer
    [ "x$answer" != "x" ] && FFLAGS="$answer"

    echoLn "enter optimisation flags for $FC compiler ($OFLAGS): "
    read answer
    [ "x$answer" != "x" ] && OFLAGS="$answer"

    checkF90Compilation

    checkF90LongLong
    if [ ${longlong} = 0 ] ; then
	echo "Your compiler does not seem to support 8-byte integers"
	echo "The compilation flag ${FPP}NO64BITS will be added to prevent their usage."
	FFLAGS="${FFLAGS} ${FPP}NO64BITS"
    fi

    FFLAGS="$OFLAGS $FFLAGS"
    echo "  Fortran code will be compiled with $FC $FFLAGS"


    echoLn "enter name of your C compiler ($CC): "
    read answer
    [ "x$answer" != "x" ] && CC="$answer"

    echoLn "enter compilation/optimisation flags for C compiler ($CFLAGS): "
    read answer
    [ "x$answer" != "x" ] && CFLAGS="$answer"

    countF90Bits
    countCBits
    if [ $c_64 != $f90_64 ] ; then
	echo "Warning: "
	if [ $f90_64 != 0 ] ; then
	    echo "F90 compiler generates 64 bit code, "
	    echo "while C compiler generates 32 bit code"
	else
	    echo "F90 compiler generates 32 bit code, "
	    echo "while C compiler generates 64 bit code"
	fi
	echoLn "you may want to change the C compilations options ($CFLAGS): "
	read answer
	[ "x$answer" != "x" ] && CFLAGS="$answer"
	echo "you also may have to recompile cfitsio with the correct options to ensure that its C and Fortran codes are consistent with each other and with 3DEX"
    fi
    echo "  C subroutines will be compiled with $CC $CFLAGS"

    echoLn "enter command for library archiving ($AR): "
    read answer
    [ "x$answer" != "x" ] && AR="$answer"

    echoLn "enter full name of cfitsio library (lib${LIBFITS}.a): "
    read answer
    [ "x$answer" != "x" ] && LIBFITS=`${BASENAME} $answer ".a" | ${SED} "s/^lib//"`

    findFITSLib $LIBDIR
    echoLn "enter location of cfitsio library ($FITSDIR): "
    read answer
    [ "x$answer" != "x" ] && FITSDIR=$answer
    fullPath FITSDIR

    lib="${FITSDIR}/lib${LIBFITS}.a"
    if [ ! -r $lib ]; then
	echo
	echo "error: fits library $lib not found"
	echo
	crashAndBurn
    fi

    # add option on where to search runtime libraries, on compilers supporting it
    if [ "x$WLRPATH" != "x" ] ; then
	WLRPATH="${WLRPATH}\$(FITSDIR)"
	LDFLAGS="${LDFLAGS} ${WLRPATH}"
    fi

    checkF90Fitsio ${lib}
    checkF90FitsioLink

}


# -----------------------------------------------------------------

#  checkNAG () {
#      [ "x`$FC -V 2>&1 | ${GREP} NAG`" != "x" ] && \
#  	FFLAGS="$FFLAGS -DLINUX -w -dusty -mismatch_all"
#  }

# -----------------------------------------------------------------

editF90Makefile () {

    echoLn "Editing top Makefile for F90 ..."
#    [ -r Makefile ] && mv Makefile Makefile.bak



    mv -f Makefile Makefile_tmp
    ${CAT} Makefile_tmp |\
	${SED} "s|^HEALPIX.*$|HEALPIX	= $HEALPIX|" |\
	${SED} "s|^F90_FC.*$|F90_FC	= $FC|" |\
	${SED} "s|^F90_FFLAGS.*$|F90_FFLAGS	= $FFLAGS|" |\
	${SED} "s|^F90_LDFLAGS.*$|F90_LDFLAGS	= $LDFLAGS|" |\
	${SED} "s|^F90_PIXFLAGS.*$|F90_PIXFLAGS	= $PIXFLAGS|" |\
	${SED} "s|^F90_CC.*$|F90_CC	= $CC|" |\
	${SED} "s|^F90_CFLAGS.*$|F90_CFLAGS	= $CFLAGS|" |\
	${SED} "s|^FITSDIR.*$|FITSDIR	= $FITSDIR|" |\
	${SED} "s|^LIBFITS.*$|LIBFITS	= $LIBFITS|" |\
	${SED} "s|^F90_BINDIR.*$|F90_BINDIR	= $F90_BINDIR|" |\
	${SED} "s|^F90_INCDIR.*$|F90_INCDIR	= $F90_INCDIR|" |\
	${SED} "s|^F90_LIBDIR.*$|F90_LIBDIR	= $F90_LIBDIR|" |\
	${SED} "s|^F90_AR.*$|F90_AR        = $AR|" |\
	${SED} "s|^F90_FFTSRC.*$|F90_FFTSRC	= $FFTSRC|" |\
	${SED} "s|^F90_ADDUS.*$|F90_ADDUS	= $ADDUS|" |\
####	${SED} "s|^F90_PARALL.*$|F90_PARALL	= $PARALL|" |\
	${SED} "s|^F90_MOD.*$|F90_MOD	= $MOD|" |\
	${SED} "s|^F90_FTYPE.*$|F90_FTYPE	= $FTYPE|" |\
	${SED} "s|^F90_PPFLAGS.*$|F90_PPFLAGS	= $PPFLAGS|" |\
	${SED} "s|^F90_OS.*$|F90_OS	= $OS|" |\
	${SED} "s|^F90_I8FLAG.*$|F90_I8FLAG  = $FI8FLAG|" |\
	${SED} "s|^ALL\(.*\) f90-void\(.*\)|ALL\1 f90-all\2|" |\
	${SED} "s|^TESTS\(.*\) f90-void\(.*\)|TESTS\1 f90-test\2|" |\
	${SED} "s|^CLEAN\(.*\) f90-void\(.*\)|CLEAN\1 f90-clean\2|" |\
	${SED} "s|^DISTCLEAN\(.*\) f90-void\(.*\)|DISTCLEAN\1 f90-distclean\2|" |\
	${SED} "s|^TIDY\(.*\) f90-void\(.*\)|TIDY\1 f90-tidy\2|" > Makefile


# 	if [ "x$DIRSUFF" != "x" ] ; then
# 	    if [ "x$DIRSUFF" != "x.in" ] ; then
# 		${CP} Makefile Makefile$DIRSUFF
# 	    fi
# 	fi

    echo " done."
    edited_makefile=1

}

# -----------------------------------------------------------------

generateConfF90File () {
	echo "Generating $DMPX_CONF_F90"

    	echo "# F90 configuration for 3DEX `date`" > $DMPX_CONF_F90
# put 3DEX variable back into F90_BINDIR
    dollar="$"
    F90_BINDIR_SHORT=`echo $F90_BINDIR | sed "s|$DMPX|{DMPX}|g"`
    F90_BINDIR_SHORT="${dollar}${F90_BINDIR_SHORT}"

    case $SHELL in
    csh|tcsh)
	${CAT} <<EOF >>$DMPX_CONF_F90
setenv HEXE    ${F90_BINDIR_SHORT}
setenv PATH    \${HEXE}:\${PATH}
EOF
	;;
    sh|ksh|bash)
	${CAT} <<EOF >>$DMPX_CONF_F90
HEXE=${F90_BINDIR_SHORT}
PATH="\${HEXE}:\${PATH}"
export HEXE PATH
EOF
	;;
    *)
	echo "Shell $SHELL not supported yet."
	${RM}  $DMPX_CONF_F90
	;;
    esac
}

# -----------------------------------------------------------------

offerF90Compilation () {
    echo "F90 Configuration finished."
    echo "You can run \"(GNU)make\" to build the package,"
    echo "        and \"(GNU)make test\" to test it."
    echoLn "You can also choose to build the package right now from here (Y|n): "
    read answer
    if [ "x$answer" != "xn" -a  "x$answer" != "xN" ]; then
    # find out make command
	askMake
    # make compilation
	${MAKE}   || crashAndBurn
	${MAKE} test || crashAndBurn
    #
	echo
	echo
	echo "F90 package installed !"
	echo
    fi
}

# -----------------------------------------------------------------

f90_config () {
    DMPX_CONF_F90=$1
    setF90Defaults
    askUserF90
    GuessCompiler
    askUserMisc
    askOpenMP
    #makeProfile
    generateConfF90File
    editF90Makefile
    [ $NOPROFILEYET = 1 ] && installProfile
#    offerF90Compilation
}

#=====================================
#=========== Check Configuration ===========
#=====================================
checkConfFiles () {

    echo "Currently, the configuration files created are :"
    echo "__________________________________________________________________"
    for conffile in ${DMPX_CONF_DIR}/*; do
	echo "${conffile} : "
	cat ${conffile}
	echo
    done
    echo "__________________________________________________________________"

}

#=====================================
#=========== Top package ===========
#=====================================

#-------------
installProfile () {
    # will modity user's configuration file to invoke 3DEX configuration

    case $SHELL in
    sh|ksh|bash)
	prof="${HOME}/.profile"
	comd="[ -r ${DMPX_CONF_MAIN} ] && . ${DMPX_CONF_MAIN}";;
    csh)
	prof="${HOME}/.cshrc"
	comd="if ( -e ${DMPX_CONF_MAIN} ) source ${DMPX_CONF_MAIN}";;
    tcsh)
	prof="${HOME}/.tcshrc"
	[ ! -r $prof -a -r "${HOME}/.cshrc" ] && prof="${HOME}/.cshrc"
	comd="if ( -e ${DMPX_CONF_MAIN} ) source ${DMPX_CONF_MAIN}";;
    *) ;;
    esac
    [ ! -r $prof ] && touch $prof
    # do not do edition if it was previously done
    check=`${GREP} ${DMPX_CONF_MAIN} $prof | ${WC} -l`
    if [ $check -eq 0 ]; then
	${CAT} <<EOF

The following line should be inserted into your home shell profile ($prof):

  $comd

 Where the file ${DMPX_CONF_MAIN} contains:
EOF
${CAT} ${DMPX_CONF_MAIN}

	echo ""
	echoLn "Do you want this modification to be done (y|N)? "
	read answer
	if [ "x$answer" = "xy" -o "x$answer" = "xY" ]; then
	    ${CP} $prof ${prof}".save"
	    echo "" >> $prof
	    echo "# modifications by 3DEXautoconf ${DMPXVERSION}" >> $prof
	    echo $comd >> $prof
	    echo "Modification done and previous shell profile saved"
	    echo "as ${prof}.save."
	fi
    else
	echo "Your home shell profile ($prof)"
	echo "has already been edited."
    fi
    NOPROFILEYET=0
}
#-------------
makeTopConf(){

    mkdir -p ${DMPX_CONF_DIR}

    case $SHELL in
    sh|ksh|bash)
	cat<<EOF >| ${DMPX_CONF_MAIN}
# configuration for 3DEX $DMPXVERSION
DMPX=${DMPX} ; export DMPX 
DMPX_CONF_DIR=${DMPX_CONF_DIR}
if [ -r ${DMPX_CONF_F90} ] ; then . ${DMPX_CONF_F90} ; fi
EOF
    echo ' ' ;;
    csh|tcsh)
	cat<<EOF >| ${DMPX_CONF_MAIN}
# configuration for 3DEX $DMPXVERSION
setenv DMPX $DMPX
setenv DMPX_CONF_DIR ${DMPX_CONF_DIR}
if ( -e ${DMPX_CONF_F90} ) source ${DMPX_CONF_F90}
EOF
    echo ' '  ;;
    *) ;;
    esac
      
}
#-------------
readyTopMakefile () {

    # backup name
    sdate=`date +%s`
    [ "x${sdate}" = "x" ] && sdate="1"
    mkbk=Makefile_bk${sdate}
    [ -s ${mkbk} ] && mkbk="${mkbk}a"

    if [ -s Makefile ] ; then
	${CP} -f Makefile ${mkbk}
    else
	if [ ! -r Makefile.in ] ; then
	    echo "top makefile template (Makefile.in) was not found. Can not proceed."
	    crashAndBurn
	fi
	${CP} -f Makefile.in Makefile
    fi

}
#-------------
restartFromScratch () {

    echo "Removing Main Makefile"
    ${RM} Makefile
    echo "Removing configuration files in " ${DMPX_CONF_DIR}
    for hfile in ${DMPX_CONF_MAIN} ${DMPX_CONF_F90}; do
	eval thisfile=${hfile}
        ${RM} ${thisfile}
    done
    echo "Removing configuration directory: " ${DMPX_CONF_DIR}
    ${RMDIR} ${DMPX_CONF_DIR}

}
#-------------
setTopDefaults() {

    AWK="awk"
    BASENAME="basename"
    CAT="cat"
    CP="cp"
    DEVNULL="/dev/null"
    DIRNAME="dirname"
    FILE="file"
    GREP="grep"
    HEAD="head" # introduced 2008-11-21
    LS="ls"
    MAKE="make"
    NM="nm"
    PRINTF="printf"
    PWD="pwd"
    RM="/bin/rm -f"
    RMDIR="rmdir"
    SED="sed"
    WC="wc"
    OS=`uname -s`

    NOPROFILEYET=1
    SHELL=`${BASENAME} ${SHELL-/bin/sh}`

    MAKESET=0

    LIBFITS="cfitsio"
    FITSDIR="/usr/local/lib"

    edited_makefile=0


    DMPX_VERSION=`echo $DMPXVERSION | ${SED} "s|\.|_|g"`
    DMPX_CONF_DIR_HOME=${HOME}/.3DEX/${DMPX_VERSION}_${OS}
    DMPX_CONF_DIR_INPLACE=${HOME}/.3DEX/confdir/${DMPX_VERSION}_${OS}

}

setConfDir () {

    case $SHELL in
    sh|ksh|bash)
        suffix=sh;;
    csh|tcsh)
	suffix=csh;;
    *) ;;
    esac

    DMPX_CONF_MAIN=$DMPX_CONF_DIR/config
    DMPX_CONF_F90=\${DMPX_CONF_DIR}/f90.${suffix}
    
}

#----------------------------------------
