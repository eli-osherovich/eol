subroutine mcsrch(n,x,f,g,s,stp,ftol,xtol,maxfev,info,nfev, gtol, stpmin, stpmax, wa)
  implicit none
  
  integer, intent(IN)                      :: n
  double precision, intent(IN OUT)         :: x(n)
  double precision, intent(IN OUT)         :: f
  double precision, intent(IN)             :: g(n)
  double precision, intent(IN)             :: s(n)
  double precision, intent(IN OUT)         :: stp
  double precision, intent(IN)             :: ftol
  double precision, intent(OUT)            :: xtol
  integer, intent(IN OUT)                  :: maxfev
  integer, intent(OUT)                     :: info
  integer, intent(OUT)                     :: nfev
  double precision, intent(in)            :: gtol, stpmin, stpmax
  double precision, intent(OUT)            :: wa(n)
  

  
  
  !                     SUBROUTINE MCSRCH
  
  !     A slight modification of the subroutine CSRCH of More' and Thuente.
  !     The changes are to allow reverse communication, and do not affect
  !     the performance of the routine.

  !     THE PURPOSE OF MCSRCH IS TO FIND A STEP WHICH SATISFIES
  !     A SUFFICIENT DECREASE CONDITION AND A CURVATURE CONDITION.

  !     AT EACH STAGE THE SUBROUTINE UPDATES AN INTERVAL OF
  !     UNCERTAINTY WITH ENDPOINTS STX AND STY. THE INTERVAL OF
  !     UNCERTAINTY IS INITIALLY CHOSEN SO THAT IT CONTAINS A
  !     MINIMIZER OF THE MODIFIED FUNCTION

  !          F(X+STP*S) - F(X) - FTOL*STP*(GRADF(X)'S).

  !     IF A STEP IS OBTAINED FOR WHICH THE MODIFIED FUNCTION
  !     HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE DERIVATIVE,
  !     THEN THE INTERVAL OF UNCERTAINTY IS CHOSEN SO THAT IT
  !     CONTAINS A MINIMIZER OF F(X+STP*S).

  !     THE ALGORITHM IS DESIGNED TO FIND A STEP WHICH SATISFIES
  !     THE SUFFICIENT DECREASE CONDITION

  !           F(X+STP*S) .LE. F(X) + FTOL*STP*(GRADF(X)'S),

  !     AND THE CURVATURE CONDITION

  !           ABS(GRADF(X+STP*S)'S)) .LE. GTOL*ABS(GRADF(X)'S).

  !     IF FTOL IS LESS THAN GTOL AND IF, FOR EXAMPLE, THE FUNCTION
  !     IS BOUNDED BELOW, THEN THERE IS ALWAYS A STEP WHICH SATISFIES
  !     BOTH CONDITIONS. IF NO STEP CAN BE FOUND WHICH SATISFIES BOTH
  !     CONDITIONS, THEN THE ALGORITHM USUALLY STOPS WHEN ROUNDING
  !     ERRORS PREVENT FURTHER PROGRESS. IN THIS CASE STP ONLY
  !     SATISFIES THE SUFFICIENT DECREASE CONDITION.

  !     THE SUBROUTINE STATEMENT IS

  !        SUBROUTINE MCSRCH(N,X,F,G,S,STP,FTOL,XTOL, MAXFEV,INFO,NFEV,WA)
  !     WHERE

  !       N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER
  !         OF VARIABLES.

  !       X IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE
  !         BASE POINT FOR THE LINE SEARCH. ON OUTPUT IT CONTAINS
  !         X + STP*S.

  !       F IS A VARIABLE. ON INPUT IT MUST CONTAIN THE VALUE OF F
  !         AT X. ON OUTPUT IT CONTAINS THE VALUE OF F AT X + STP*S.

  !       G IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE
  !         GRADIENT OF F AT X. ON OUTPUT IT CONTAINS THE GRADIENT
  !         OF F AT X + STP*S.

  !       S IS AN INPUT ARRAY OF LENGTH N WHICH SPECIFIES THE
  !         SEARCH DIRECTION.

  !       STP IS A NONNEGATIVE VARIABLE. ON INPUT STP CONTAINS AN
  !         INITIAL ESTIMATE OF A SATISFACTORY STEP. ON OUTPUT
  !         STP CONTAINS THE FINAL ESTIMATE.

  !       FTOL AND GTOL ARE NONNEGATIVE INPUT VARIABLES. (In this reverse
  !         communication implementation GTOL is defined in a COMMON
  !         statement.) TERMINATION OCCURS WHEN THE SUFFICIENT DECREASE
  !         CONDITION AND THE DIRECTIONAL DERIVATIVE CONDITION ARE
  !         SATISFIED.

  !       XTOL IS A NONNEGATIVE INPUT VARIABLE. TERMINATION OCCURS
  !         WHEN THE RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
  !         IS AT MOST XTOL.

  !       STPMIN AND STPMAX ARE NONNEGATIVE INPUT VARIABLES WHICH
  !         SPECIFY LOWER AND UPPER BOUNDS FOR THE STEP. (In this reverse
  !         communication implementatin they are defined in a COMMON
  !         statement).

  !       MAXFEV IS A POSITIVE INTEGER INPUT VARIABLE. TERMINATION
  !         OCCURS WHEN THE NUMBER OF CALLS TO FCN IS AT LEAST
  !         MAXFEV BY THE END OF AN ITERATION.

  !       INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:

  !         INFO = 0  IMPROPER INPUT PARAMETERS.

  !         INFO =-1  A RETURN IS MADE TO COMPUTE THE FUNCTION AND GRADIENT.

  !         INFO = 1  THE SUFFICIENT DECREASE CONDITION AND THE
  !                   DIRECTIONAL DERIVATIVE CONDITION HOLD.

  !         INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY
  !                   IS AT MOST XTOL.

  !         INFO = 3  NUMBER OF CALLS TO FCN HAS REACHED MAXFEV.

  !         INFO = 4  THE STEP IS AT THE LOWER BOUND STPMIN.

  !         INFO = 5  THE STEP IS AT THE UPPER BOUND STPMAX.

  !         INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS.
  !                   THERE MAY NOT BE A STEP WHICH SATISFIES THE
  !                   SUFFICIENT DECREASE AND CURVATURE CONDITIONS.
  !                   TOLERANCES MAY BE TOO SMALL.

  !       NFEV IS AN INTEGER OUTPUT VARIABLE SET TO THE NUMBER OF
  !         CALLS TO FCN.

  !       WA IS A WORK ARRAY OF LENGTH N.

  !     SUBPROGRAMS CALLED

  !       MCSTEP

  !       FORTRAN-SUPPLIED...ABS,MAX,MIN

  !     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
  !     JORGE J. MORE', DAVID J. THUENTE

  !     **********



  INTEGER :: infoc,j
  LOGICAL :: brackt,stage1
  DOUBLE PRECISION :: dg,dgm,dginit,dgtest,dgx,dgxm,dgy,dgym,  &
       finit,ftest1,fm,fx,fxm,fy,fym,p5,p66,stx,sty,  &
       stmin,stmax,width,width1,xtrapf,zero
  DATA p5,p66,xtrapf,zero /0.5D0,0.66D0,4.0D0,0.0D0/

  IF(info == -1) GO TO 45
  infoc = 1

  !     CHECK THE INPUT PARAMETERS FOR ERRORS.

  IF (n <= 0 .OR. stp <= zero .OR. ftol < zero .OR.  &
       gtol < zero .OR. xtol < zero .OR. stpmin < zero  &
       .OR. stpmax < stpmin .OR. maxfev <= 0) RETURN

  !     COMPUTE THE INITIAL GRADIENT IN THE SEARCH DIRECTION
  !     AND CHECK THAT S IS A DESCENT DIRECTION.

  dginit = zero
  do  j = 1, n
     dginit = dginit + g(j)*s(j)
  end do
  if (dginit >= zero) then
     info = 7
     return
  end if

  !     INITIALIZE LOCAL VARIABLES.

  brackt = .false.
  stage1 = .true.
  nfev = 0
  finit = f
  dgtest = ftol*dginit
  width = stpmax - stpmin
  width1 = width/p5
  DO  j = 1, n
     wa(j) = x(j)
  END DO

  !     THE VARIABLES STX, FX, DGX CONTAIN THE VALUES OF THE STEP,
  !     FUNCTION, AND DIRECTIONAL DERIVATIVE AT THE BEST STEP.
  !     THE VARIABLES STY, FY, DGY CONTAIN THE VALUE OF THE STEP,
  !     FUNCTION, AND DERIVATIVE AT THE OTHER ENDPOINT OF
  !     THE INTERVAL OF UNCERTAINTY.
  !     THE VARIABLES STP, F, DG CONTAIN THE VALUES OF THE STEP,
  !     FUNCTION, AND DERIVATIVE AT THE CURRENT STEP.

  stx = zero
  fx = finit
  dgx = dginit
  sty = zero
  fy = finit
  dgy = dginit

  !     START OF ITERATION.

30 CONTINUE

  !        SET THE MINIMUM AND MAXIMUM STEPS TO CORRESPOND
  !        TO THE PRESENT INTERVAL OF UNCERTAINTY.

  IF (brackt) THEN
     stmin = MIN(stx,sty)
     stmax = MAX(stx,sty)
  ELSE
     stmin = stx
     stmax = stp + xtrapf*(stp - stx)
  END IF

  !        FORCE THE STEP TO BE WITHIN THE BOUNDS STPMAX AND STPMIN.

  stp = MAX(stp,stpmin)
  stp = MIN(stp,stpmax)

  !        IF AN UNUSUAL TERMINATION IS TO OCCUR THEN LET
  !        STP BE THE LOWEST POINT OBTAINED SO FAR.

  IF ((brackt .AND. (stp <= stmin .OR. stp >= stmax))  &
       .OR. nfev >= maxfev-1 .OR. infoc == 0  &
       .OR. (brackt .AND. stmax-stmin <= xtol*stmax)) stp = stx

  !        EVALUATE THE FUNCTION AND GRADIENT AT STP
  !        AND COMPUTE THE DIRECTIONAL DERIVATIVE.
  !        We return to main program to obtain F and G.

  DO  j = 1, n
     x(j) = wa(j) + stp*s(j)
  END DO
  info=-1
  RETURN

45 info=0
  nfev = nfev + 1
  dg = zero
  DO  j = 1, n
     dg = dg + g(j)*s(j)
  END DO
  ftest1 = finit + stp*dgtest

  !        TEST FOR CONVERGENCE.

  IF ((brackt .AND. (stp <= stmin .OR. stp >= stmax)) .OR. infoc == 0) info = 6
  IF (stp == stpmax .AND. f <= ftest1 .AND. dg <= dgtest) info = 5
  IF (stp == stpmin .AND. (f > ftest1 .OR. dg >= dgtest)) info = 4
  IF (nfev >= maxfev) info = 3
  IF (brackt .AND. stmax-stmin <= xtol*stmax) info = 2
  IF (f <= ftest1 .AND. ABS(dg) <= gtol*(-dginit)) info = 1

  !        CHECK FOR TERMINATION.

  IF (info /= 0) RETURN

  !        IN THE FIRST STAGE WE SEEK A STEP FOR WHICH THE MODIFIED
  !        FUNCTION HAS A NONPOSITIVE VALUE AND NONNEGATIVE DERIVATIVE.

  IF (stage1 .AND. f <= ftest1 .AND.  &
       dg >= MIN(ftol,gtol)*dginit) stage1 = .false.

  !        A MODIFIED FUNCTION IS USED TO PREDICT THE STEP ONLY IF
  !        WE HAVE NOT OBTAINED A STEP FOR WHICH THE MODIFIED
  !        FUNCTION HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE
  !        DERIVATIVE, AND IF A LOWER FUNCTION VALUE HAS BEEN
  !        OBTAINED BUT THE DECREASE IS NOT SUFFICIENT.

  IF (stage1 .AND. f <= fx .AND. f > ftest1) THEN

     !           DEFINE THE MODIFIED FUNCTION AND DERIVATIVE VALUES.

     fm = f - stp*dgtest
     fxm = fx - stx*dgtest
     fym = fy - sty*dgtest
     dgm = dg - dgtest
     dgxm = dgx - dgtest
     dgym = dgy - dgtest

     !           CALL CSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
     !           AND TO COMPUTE THE NEW STEP.

     CALL mcstep(stx,fxm,dgxm,sty,fym,dgym,stp,fm,dgm, brackt,stmin,stmax,infoc)

     !           RESET THE FUNCTION AND GRADIENT VALUES FOR F.

     fx = fxm + stx*dgtest
     fy = fym + sty*dgtest
     dgx = dgxm + dgtest
     dgy = dgym + dgtest
  ELSE

     !           CALL MCSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY
     !           AND TO COMPUTE THE NEW STEP.

     CALL mcstep(stx,fx,dgx,sty,fy,dgy,stp,f,dg, brackt,stmin,stmax,infoc)
  END IF

  !        FORCE A SUFFICIENT DECREASE IN THE SIZE OF THE
  !        INTERVAL OF UNCERTAINTY.

  IF (brackt) THEN
     IF (ABS(sty-stx) >= p66*width1) stp = stx + p5*(sty - stx)
     width1 = width
     width = ABS(sty-stx)
  END IF

  !        END OF ITERATION.

  GO TO 30

  !     LAST LINE OF SUBROUTINE MCSRCH.

END SUBROUTINE mcsrch

subroutine mcstep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt, stpmin,stpmax,info)

  implicit none

  double precision, intent(IN OUT)         :: stx
  double precision, intent(IN OUT)         :: fx
  double precision, intent(IN OUT)         :: dx
  double precision, intent(OUT)            :: sty
  double precision, intent(OUT)            :: fy
  double precision, intent(IN OUT)         :: dy
  double precision, intent(IN OUT)         :: stp
  double precision, intent(IN)             :: fp
  double precision, intent(IN)             :: dp
  logical, intent(OUT)                     :: brackt
  double precision, intent(IN)             :: stpmin
  double precision, intent(IN)             :: stpmax
  integer, intent(OUT)                     :: info


  logical :: bound

  !     SUBROUTINE MCSTEP

  !     THE PURPOSE OF MCSTEP IS TO COMPUTE A SAFEGUARDED STEP FOR
  !     A LINESEARCH AND TO UPDATE AN INTERVAL OF UNCERTAINTY FOR
  !     A MINIMIZER OF THE FUNCTION.

  !     THE PARAMETER STX CONTAINS THE STEP WITH THE LEAST FUNCTION
  !     VALUE. THE PARAMETER STP CONTAINS THE CURRENT STEP. IT IS
  !     ASSUMED THAT THE DERIVATIVE AT STX IS NEGATIVE IN THE
  !     DIRECTION OF THE STEP. IF BRACKT IS SET TRUE THEN A
  !     MINIMIZER HAS BEEN BRACKETED IN AN INTERVAL OF UNCERTAINTY
  !     WITH ENDPOINTS STX AND STY.

  !     THE SUBROUTINE STATEMENT IS

  !       SUBROUTINE MCSTEP(STX,FX,DX,STY,FY,DY,STP,FP,DP,BRACKT,
  !                        STPMIN,STPMAX,INFO)

  !     WHERE

  !       STX, FX, AND DX ARE VARIABLES WHICH SPECIFY THE STEP,
  !         THE FUNCTION, AND THE DERIVATIVE AT THE BEST STEP OBTAINED
  !         SO FAR. THE DERIVATIVE MUST BE NEGATIVE IN THE DIRECTION
  !         OF THE STEP, THAT IS, DX AND STP-STX MUST HAVE OPPOSITE
  !         SIGNS. ON OUTPUT THESE PARAMETERS ARE UPDATED APPROPRIATELY.

  !       STY, FY, AND DY ARE VARIABLES WHICH SPECIFY THE STEP,
  !         THE FUNCTION, AND THE DERIVATIVE AT THE OTHER ENDPOINT OF
  !         THE INTERVAL OF UNCERTAINTY. ON OUTPUT THESE PARAMETERS ARE
  !         UPDATED APPROPRIATELY.

  !       STP, FP, AND DP ARE VARIABLES WHICH SPECIFY THE STEP,
  !         THE FUNCTION, AND THE DERIVATIVE AT THE CURRENT STEP.
  !         IF BRACKT IS SET TRUE THEN ON INPUT STP MUST BE
  !         BETWEEN STX AND STY. ON OUTPUT STP IS SET TO THE NEW STEP.

  !       BRACKT IS A LOGICAL VARIABLE WHICH SPECIFIES IF A MINIMIZER
  !         HAS BEEN BRACKETED. IF THE MINIMIZER HAS NOT BEEN BRACKETED
  !         THEN ON INPUT BRACKT MUST BE SET FALSE. IF THE MINIMIZER
  !         IS BRACKETED THEN ON OUTPUT BRACKT IS SET TRUE.

  !       STPMIN AND STPMAX ARE INPUT VARIABLES WHICH SPECIFY LOWER
  !         AND UPPER BOUNDS FOR THE STEP.

  !       INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:
  !         IF INFO = 1,2,3,4,5, THEN THE STEP HAS BEEN COMPUTED
  !         ACCORDING TO ONE OF THE FIVE CASES BELOW. OTHERWISE
  !         INFO = 0, AND THIS INDICATES IMPROPER INPUT PARAMETERS.

  !     SUBPROGRAMS CALLED

  !       FORTRAN-SUPPLIED ... ABS,MAX,MIN,SQRT

  !     ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983
  !     JORGE J. MORE', DAVID J. THUENTE

  DOUBLE PRECISION :: gamma,p,q,r,s,sgnd,stpc,stpf,stpq,theta

  info = 0

  !     CHECK THE INPUT PARAMETERS FOR ERRORS.

  IF ((brackt .AND. (stp <= MIN(stx,sty) .OR. stp >= MAX(stx,sty))) .OR.  &
       dx*(stp-stx) >= 0.0 .OR. stpmax < stpmin) RETURN

  !     DETERMINE IF THE DERIVATIVES HAVE OPPOSITE SIGN.

  sgnd = dp*(dx/ABS(dx))

  !     FIRST CASE. A HIGHER FUNCTION VALUE.
  !     THE MINIMUM IS BRACKETED. IF THE CUBIC STEP IS CLOSER
  !     TO STX THAN THE QUADRATIC STEP, THE CUBIC STEP IS TAKEN,
  !     ELSE THE AVERAGE OF THE CUBIC AND QUADRATIC STEPS IS TAKEN.

  IF (fp > fx) THEN
     info = 1
     bound = .true.
     theta = 3*(fx - fp)/(stp - stx) + dx + dp
     s = MAX(ABS(theta),ABS(dx),ABS(dp))
     gamma = s*SQRT((theta/s)**2 - (dx/s)*(dp/s))
     IF (stp < stx) gamma = -gamma
     p = (gamma - dx) + theta
     q = ((gamma - dx) + gamma) + dp
     r = p/q
     stpc = stx + r*(stp - stx)
     stpq = stx + ((dx/((fx-fp)/(stp-stx)+dx))/2)*(stp - stx)
     IF (ABS(stpc-stx) < ABS(stpq-stx)) THEN
        stpf = stpc
     ELSE
        stpf = stpc + (stpq - stpc)/2
     END IF
     brackt = .true.

     !     SECOND CASE. A LOWER FUNCTION VALUE AND DERIVATIVES OF
     !     OPPOSITE SIGN. THE MINIMUM IS BRACKETED. IF THE CUBIC
     !     STEP IS CLOSER TO STX THAN THE QUADRATIC (SECANT) STEP,
     !     THE CUBIC STEP IS TAKEN, ELSE THE QUADRATIC STEP IS TAKEN.

  ELSE IF (sgnd < 0.0) THEN
     info = 2
     bound = .false.
     theta = 3*(fx - fp)/(stp - stx) + dx + dp
     s = MAX(ABS(theta),ABS(dx),ABS(dp))
     gamma = s*SQRT((theta/s)**2 - (dx/s)*(dp/s))
     IF (stp > stx) gamma = -gamma
     p = (gamma - dp) + theta
     q = ((gamma - dp) + gamma) + dx
     r = p/q
     stpc = stp + r*(stx - stp)
     stpq = stp + (dp/(dp-dx))*(stx - stp)
     IF (ABS(stpc-stp) > ABS(stpq-stp)) THEN
        stpf = stpc
     ELSE
        stpf = stpq
     END IF
     brackt = .true.

     !     THIRD CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
     !     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DECREASES.
     !     THE CUBIC STEP IS ONLY USED IF THE CUBIC TENDS TO INFINITY
     !     IN THE DIRECTION OF THE STEP OR IF THE MINIMUM OF THE CUBIC
     !     IS BEYOND STP. OTHERWISE THE CUBIC STEP IS DEFINED TO BE
     !     EITHER STPMIN OR STPMAX. THE QUADRATIC (SECANT) STEP IS ALSO
     !     COMPUTED AND IF THE MINIMUM IS BRACKETED THEN THE THE STEP
     !     CLOSEST TO STX IS TAKEN, ELSE THE STEP FARTHEST AWAY IS TAKEN.

  ELSE IF (ABS(dp) < ABS(dx)) THEN
     info = 3
     bound = .true.
     theta = 3*(fx - fp)/(stp - stx) + dx + dp
     s = MAX(ABS(theta),ABS(dx),ABS(dp))

     !        THE CASE GAMMA = 0 ONLY ARISES IF THE CUBIC DOES NOT TEND
     !        TO INFINITY IN THE DIRECTION OF THE STEP.

     gamma = s*SQRT(MAX(0.0D0,(theta/s)**2 - (dx/s)*(dp/s)))
     IF (stp > stx) gamma = -gamma
     p = (gamma - dp) + theta
     q = (gamma + (dx - dp)) + gamma
     r = p/q
     IF (r < 0.0 .AND. gamma /= 0.0) THEN
        stpc = stp + r*(stx - stp)
     ELSE IF (stp > stx) THEN
        stpc = stpmax
     ELSE
        stpc = stpmin
     END IF
     stpq = stp + (dp/(dp-dx))*(stx - stp)
     IF (brackt) THEN
        IF (ABS(stp-stpc) < ABS(stp-stpq)) THEN
           stpf = stpc
        ELSE
           stpf = stpq
        END IF
     ELSE
        IF (ABS(stp-stpc) > ABS(stp-stpq)) THEN
           stpf = stpc
        ELSE
           stpf = stpq
        END IF
     END IF

     !     FOURTH CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE
     !     SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DOES
     !     NOT DECREASE. IF THE MINIMUM IS NOT BRACKETED, THE STEP
     !     IS EITHER STPMIN OR STPMAX, ELSE THE CUBIC STEP IS TAKEN.

  ELSE
     info = 4
     bound = .false.
     IF (brackt) THEN
        theta = 3*(fp - fy)/(sty - stp) + dy + dp
        s = MAX(ABS(theta),ABS(dy),ABS(dp))
        gamma = s*SQRT((theta/s)**2 - (dy/s)*(dp/s))
        IF (stp > sty) gamma = -gamma
        p = (gamma - dp) + theta
        q = ((gamma - dp) + gamma) + dy
        r = p/q
        stpc = stp + r*(sty - stp)
        stpf = stpc
     ELSE IF (stp > stx) THEN
        stpf = stpmax
     ELSE
        stpf = stpmin
     END IF
  END IF

  !     UPDATE THE INTERVAL OF UNCERTAINTY. THIS UPDATE DOES NOT
  !     DEPEND ON THE NEW STEP OR THE CASE ANALYSIS ABOVE.

  IF (fp > fx) THEN
     sty = stp
     fy = fp
     dy = dp
  ELSE
     IF (sgnd < 0.0) THEN
        sty = stx
        fy = fx
        dy = dx
     END IF
     stx = stp
     fx = fp
     dx = dp
  END IF

  !     COMPUTE THE NEW STEP AND SAFEGUARD IT.

  stpf = MIN(stpmax,stpf)
  stpf = MAX(stpmin,stpf)
  stp = stpf
  IF (brackt .AND. bound) THEN
     IF (sty > stx) THEN
        stp = MIN(stx+0.66*(sty-stx),stp)
     ELSE
        stp = MAX(stx+0.66*(sty-stx),stp)
     END IF
  END IF
  RETURN

  !     LAST LINE OF SUBROUTINE MCSTEP.

END SUBROUTINE mcstep

