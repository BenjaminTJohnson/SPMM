
!** domain based melting model (rather than sparse)
!** update: read in any melting file step, and start from there.
PROGRAM domain_melt4
  IMPLICIT NONE
  !** parameters
  INTEGER,PARAMETER :: FID = 13, OUTFID = 17
  
  !** main storage variables
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: domain
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: pos, oTCM
  !** A liquid dipole's move depends only on the cells within two of it, so it is looked at again
  !** only after one of those cells has changed ('dirty').  'who' gives the dipole in each cell.
  INTEGER, DIMENSION(:,:,:), ALLOCATABLE :: who
  LOGICAL(1), DIMENSION(:), ALLOCATABLE :: dirty

  !** local variables
  INTEGER, DIMENSION(3) :: cpos, tpos, xpos, TCM, LCM
  INTEGER :: OUTPUT,i,j,k,PDIRS(27,3),XMAX,YMAX,ZMAX,n,NDIP,m,mx,niter, ctype, ttype, mm, nn, holes(27), ii
  INTEGER :: NNicecnt, NNicecntmax, NNliqcnt, icecnt, liqcnt, NNholcnt, meltcnt, movecnt, moveiter

  INTEGER :: xNNicecnt(27), xNNliqcnt(27), xNNholcnt(27), xtype, xoffset, yoffset, zoffset

  REAL(8) :: rx, fliq, fliqold, tmp, c2TCM, gwt

  !** Output and melt-step control, optionally set in group &melt of the file spmm.nml.
  !**   targets     melted mass fractions at which a frame is written, ascending (0 and 1 are always written)
  !**   dfmax       largest fraction of the particle melted in one step while below f_late
  !**   dfmax_late  the same, from f_late onwards
  !** Each step melts just enough to reach the next target, or the cap if that is smaller.
  INTEGER, PARAMETER :: maxtarg = 100
  REAL(8) :: targets(maxtarg), dfmax = 0.005d0, dfmax_late = 0.02d0, f_late = 0.20d0, nextt
  INTEGER :: seed = 12345, ntarg, itarget, nelig, need, left, ios
  LOGICAL :: ex
  NAMELIST /melt/ targets, dfmax, dfmax_late, f_late, seed

  !** mass conservation: ice (916.7 kg/m3 at 0C) melts to a smaller volume of water (999.84 kg/m3),
  !** so each melted ice dipole yields rho_ratio water dipoles.  The shortfall is accumulated in
  !** 'debt'; each time it reaches one whole dipole the melting dipole is removed instead of
  !** becoming liquid.  Removed dipoles are flagged with pos(n,1) = -1 and are not written out.
  REAL(8), PARAMETER :: rho_ratio = 916.7d0/999.84d0
  REAL(8) :: debt
  INTEGER :: nliq, ngone   !** liquid dipoles present, dipoles removed on melting

  !** connected bodies: dipoles (ice and liquid together) joined through any of their 26 neighbours.
  !** Labelled whenever a frame is written; body 1 is the largest, so detached droplets are body >= 2.
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: cq          !** flood-fill queue of positions
  INTEGER, DIMENSION(:), ALLOCATABLE :: body, bsize, border, brank
  INTEGER :: nbody
  CHARACTER(len=64) :: outfn, inpfn
  
    
  OUTPUT = 1
  !** target output melt fractions: 1% steps to 10%, 2.5% steps to 20%, then coarser
  targets = -1.0d0
  INQUIRE(file='spmm.nml', exist=ex)
  IF (ex) THEN
     OPEN(FID, file='spmm.nml', status='old')
     READ(FID, nml=melt, iostat=ios)
     CLOSE(FID)
     IF (ios > 0) STOP 'error reading group &melt in spmm.nml'
  END IF
  IF (targets(1) < 0.0d0) targets(1:24) = (/0.005d0, 0.01d0, 0.02d0, 0.03d0, 0.04d0, 0.05d0, 0.06d0, 0.07d0, &
       0.08d0, 0.09d0, 0.10d0, 0.125d0, 0.15d0, 0.175d0, 0.20d0, 0.25d0, 0.30d0, 0.40d0, 0.50d0, 0.60d0, &
       0.70d0, 0.80d0, 0.90d0, 0.95d0/)
  ntarg = COUNT(targets > 0.0d0 .AND. targets < 1.0d0)
  DO n = 2,ntarg
     IF (targets(n) <= targets(n-1)) STOP 'targets must be ascending, between 0 and 1'
  END DO
  itarget = 1
  CALL seed_rng(seed)

  !** compute all permutations of directions from a given point (to find NN)
  n=0
  DO i= -1,1
     DO j = -1,1
        DO k = -1,1
           n = n + 1
           PDIRS(n,1:3) = (/i,j,k/)
           !        print *, n,(/i,j,k/)
        END DO
     END DO
  END DO
  
  !** open the original shape.dat (kuo style) to get dimensions and numbe of dipoles.
  OPEN(unit=FID, file='shape.dat', status='old')
  READ(FID,*)  !** skip header line
  NDIP = 0
  XMAX = 0
  YMAX = 0
  ZMAX = 0
  DO 
     READ(FID,*,END=99) i,j,k,m
     NDIP = NDIP + 1
     IF (i > XMAX) XMAX = i
     IF (j > YMAX) YMAX = j
     IF (k > ZMAX) ZMAX = k
  END DO
99 CLOSE (FID)  !** we only need the range and number of dipoles
  
  
  !** increase the domain a bit
  XMAX = XMAX + 5
  YMAX = YMAX + 5
  ZMAX = ZMAX + 5 !** NOTE: for flat shapes, the domain needs enough space to create a sphere, so, extra offsets will be needed normal to the flat dimension
    
  !** find the maximum dimension of the domain
  mx = MAXVAL((/XMAX,YMAX,ZMAX/))
  XMAX = mx  !** make a cubic domain...
  YMAX = mx
  ZMAX = mx 
  
  xoffset = NINT(mx/2.0) - NINT(XMAX/2.0)+1
  yoffset = NINT(mx/2.0) - NINT(YMAX/2.0)+1
  zoffset = NINT(mx/2.0) - NINT(ZMAX/2.0)+1
  
  ALLOCATE(domain(XMAX,YMAX,ZMAX)) !** total 3D domain
  ALLOCATE(who(XMAX,YMAX,ZMAX), dirty(NDIP))
  who = 0
  dirty = .TRUE.
  ALLOCATE(pos(NDIP,3))            !** position of liquid or ice points 
  ALLOCATE(oTCM(NDIP,3))            !** nearest TCM for the current dipole
  ALLOCATE(cq(3,NDIP), body(NDIP), bsize(NDIP), border(NDIP), brank(NDIP))  !** connected-body labelling
  
  !** read input data !CONTINUE modifications here, to read the default file if it's the first time running
  !** or to read a melt file on subsequent times running.  Maybe add an args option to read the filename? maybe not needed.
  PRINT *, 'reading data'
  oTCM      = 0
  TCM       = 0
  domain    = -99  !** all points initially -99 (uninitialized) 
  pos       = -99

  OPEN(unit=FID, file='shape.dat', status='old')
  READ(FID,*,END=100) !** skip header
  INIT1: DO n = 1,NDIP
     READ(FID,*,END=100) i,j,k,m
     domain(i+xoffset,j+yoffset,k+zoffset) = 1     !** populate ice-occupied points in domain
     cpos           = (/i+xoffset,j+yoffset,k+zoffset/)  !** local storage of position
     pos(n,1:3)     = cpos           !** store the list of dipoles
     who(cpos(1),cpos(2),cpos(3)) = n
     TCM(1:3)      =  TCM(1:3) + cpos(1:3)  !** computation of total center of mass of the system of all points
  END DO INIT1
100 CLOSE(FID)
  fliqold = 0.0d0  !** starts fresh
  TCM = NINT(DBLE(TCM)/DBLE(NDIP))   !** total center of mass for the particle
  oTCM(:,1) = TCM(1) !** intially set the center point to the total center of mass
  oTCM(:,2) = TCM(2) !** intially set the center point to the total center of mass
  oTCM(:,3) = TCM(3) !** intially set the center point to the total center of mass
  PRINT *, 'TCM:', TCM(:), oTCM(1,:), oTCM(NDIP,:)
  PRINT *, 'done reading data', NDIP

  !** output initial file in melt file format (just in case)
           
  fliq = 0.0
  nliq  = 0
  ngone = 0
  debt  = 0.0d0
  
  PRINT '(5A)', '       ', '      iter', '   liquid', '    fraction', '   target' 
  PRINT '(A,2I9,3G12.2)', 'output0:',  niter, 0, fliq, 0.0
  IF (OUTPUT == 1) CALL write_frame()
            


  
  !** domain is now populated with the occupied points
  
  NNicecntmax = 5 !** initially
  niter = 0
  
  iterloop: DO
     niter = niter + 1
!     PRINT *, 'iteration #', niter
     !** loop over the whole domain, and check for melting criteria
     icecnt      = 0
     liqcnt      = 0
     meltcnt     = 0
     movecnt     = 0
     holes       = 0
     tmp         = 0.0d0
!     fliq        = 0.0d0
     !** MELTING BLOCK
     IF (movecnt == 0) THEN !** everything is pretty much done moving, we can melt again  .or. fliq > 0.0
        !** how many dipoles to melt this step: enough to reach the next target, at most the cap
        nelig = count_eligible()
        nextt = 1.0d0
        IF (itarget <= ntarg) nextt = targets(itarget)
        need  = MAX(1, NINT(MERGE(dfmax, dfmax_late, fliq < f_late)*NDIP))
        need  = MAX(1, MIN(need, CEILING((nextt - fliq)*NDIP - 1.0d-6)))
        left  = nelig
        meltloop: DO n = 1,NDIP
           IF (pos(n,1) < 0) CYCLE meltloop  !** dipole removed by mass conservation
           i = pos(n,1)
           j = pos(n,2)
           k = pos(n,3)
           ctype = domain(i,j,k)  !** current dipole type
           cpos(1:3) = (/i,j,k/) 
           !** deal with liquid points ** 
           IF (ctype == -1) THEN !** it's a liquid point, no need to process furhter
              liqcnt = liqcnt + 1 !** update liquid point counter for current iteration
              cycle meltloop
           END IF
           
           !** deal with ice points ** 
           
           IF (ctype == 1) THEN   !** it's an ice point
              icecnt = icecnt + 1
           END IF
           
           !** check for melting criteria
           !** 1. Check the nearest neighbors
           NNicecnt = 0
           NNliqcnt = 0
           NNholcnt = 0
           LCM      = 0
           NNloop1: DO m = 1,27
              IF (m == 14) CYCLE NNloop1  !** skip self
              tpos = cpos + pdirs(m,1:3)  !** target position for NN
              IF(ANY(tpos(:) <= 0) .OR. tpos(1) > XMAX .OR. tpos(2) > YMAX .OR. tpos(3) > ZMAX) THEN
                 CYCLE NNloop1 !** out of bounds
              END IF
              ttype = domain(tpos(1),tpos(2),tpos(3)) !** type of neighbor
              IF (ttype == 1) THEN  !** the neighbor is an icepoint
                 NNicecnt = NNicecnt + 1
                 LCM = LCM + tpos !** create a local center of mass 
              ELSEIF (ttype == -1) THEN !** the neighbor is a liquid point 
                 NNliqcnt = NNliqcnt + 1 
              ELSEIF (ttype == 0 .OR. ttype == -99) THEN
                 NNholcnt = NNholcnt + 1 !** count the holes of NNs 
                 domain(tpos(1),tpos(2),tpos(3)) = 0       !** it's a hole that has at least one ice or liquid neighbor 
              END IF
           END DO NNloop1
!!$           IF (NNicecnt > 0) THEN 
       !       oTCM(n,:) = NINT(DBLE(LCM)/DBLE(NNicecnt))
!              c2TCM = SQRT(DBLE((tpos(1)-TCM(1))**2)+DBLE((tpos(2)-TCM(2))**2) + DBLE((tpos(3)-TCM(3))**2))
!!$           ELSE
!!$              oTCM(n,:) = TCM(:)
!!$           END IF
!!$           !** 2. check the melt criteria
           IF (ctype == 1 .AND. NNicecnt <= NNicecntmax) THEN 
              CALL RANDOM_NUMBER(rx)
              !** pick 'need' of the 'left' eligible dipoles at random, each equally likely
              left = left - 1
              IF (need > 0 .AND. rx*MAX(left + 1, 1) < need) THEN
                 need = need - 1
                 !** 3. melt the current particle (conserving mass, see rho_ratio)
                 debt = debt + (1.0d0 - rho_ratio)
                 IF (debt >= 1.0d0) THEN
                    debt          = debt - 1.0d0
                    domain(i,j,k) = 0   !** its volume is lost to the density change: leave a hole
                    pos(n,1:3)    = -1  !** flag the dipole as removed
                    who(i,j,k)    = 0
                    ngone         = ngone + 1
                 ELSE
                    domain(i,j,k) = -1  !** switch to melted
                    nliq          = nliq + 1
                 END IF
                 ctype         = -1
                 meltcnt       = meltcnt + 1 !** only keep track of melting events
                 CALL mark_dirty(cpos)
              END IF
           END IF
        END DO meltloop
        IF (nelig == 0) THEN !** nothing was eligible to melt during this iteration
           NNicecntmax = NNicecntmax + 1  !** increase the melt threshhold 
        END IF
     END IF
     
     moveiter = 0
     
     moveiterloop: DO
        moveiter = moveiter + 1
        movecnt = 0
        moveloop: DO n = 1,NDIP
           !** 4. determine movement ability
           IF (.NOT. dirty(n)) CYCLE moveloop  !** nothing near it has changed since it last stayed put
           dirty(n) = .FALSE.
           IF (pos(n,1) < 0) CYCLE moveloop  !** dipole removed by mass conservation
           i = pos(n,1)
           j = pos(n,2)
           k = pos(n,3)
           ctype = domain(i,j,k)  !** current dipole type
           cpos(1:3) = (/i,j,k/) 
           
           !** deal with liquid points ** 
           LCM = TCM - cpos
           oTCM(n,:) = cpos + 0.25*LCM

           IF (ctype == -1) THEN !** it's a liquid point
              c2TCM = DBLE((cpos(1)-oTCM(n,1))**2)+DBLE((cpos(2)-oTCM(n,2))**2) + DBLE((cpos(3)-oTCM(n,3))**2)
           ELSE
              CYCLE moveloop  !** only liquid point can even theoretically move
           END IF

           !** check for movement criteria
           !** 5. Check the nearest neighbors (this can't be done in the melt loop, it wouldn't update properly!
           !**    unless... I integrate incremental movement into the melt loop itself. hmm.) 
           NNicecnt = 0
           NNliqcnt = 0
           NNholcnt = 0
           xNNliqcnt(:) = 0
           xNNicecnt(:) = 0
           xNNholcnt(:) = 0
           mx = 0
           
           NNloop2: DO m = 1,27
              IF (m == 14) CYCLE NNloop2  !** skip self
              tpos = cpos + pdirs(m,1:3)  !** target position for NN
              IF(ANY(tpos(:) <= 0) .OR. tpos(1) > XMAX .OR. tpos(2) > YMAX .OR. tpos(3) > ZMAX) THEN
                 CYCLE NNloop2 !** out of bounds
              END IF
              ttype = domain(tpos(1),tpos(2),tpos(3)) !** type of neighbor
              IF (ttype == 1) THEN  !** the neighbor is an icepoint
                 NNicecnt = NNicecnt + 1
                 !** the current point is unable to move if there are any ice points as neighbors!
                 CYCLE moveloop 
              ELSEIF (ttype == -1) THEN !** the neighbor is a liquid point 
                 NNliqcnt = NNliqcnt + 1 
                 !** we don't care about liquid point neighbors too much for movement, we want to find a hole
              ELSEIF (ttype == 0 .OR. ttype == -99) THEN
                 NNholcnt = NNholcnt + 1 !** count the hole neighbors of the current liquid point

                 domain(tpos(1),tpos(2),tpos(3)) = 0       !** update the domain knowledge for the target position (as a hole)
                 ttype = 0 
                 !** need an orphan check here
                 tmp = DBLE((tpos(1)-oTCM(n,1))**2)+DBLE((tpos(2)-oTCM(n,2))**2) + DBLE((tpos(3)-oTCM(n,3))**2) !**linear distance to dipoles
                 !** is this hole closer to the total center of mass than the current point itself? 
                 IF (tmp < c2TCM) THEN
                    
                    !** if so, it's a candidate hole for the free liquid point to move into
                    !** check the hole's neighbors to find a preferable hole to move into.  always prefer a hole next to ice. 
                    NNloop3: DO mm = 1,27
                       IF (mm == 14) CYCLE NNloop3  !** skip self
                       xpos = tpos + pdirs(mm,1:3)  !** target position for NN
                       IF(ANY(xpos(:) <= 0) .OR. &
                            xpos(1) > XMAX  .OR. &    
                            xpos(2) > YMAX  .OR. &
                            xpos(3) > ZMAX  .OR. &
                            ALL(xpos == cpos) ) THEN  !** skip cpos
                          CYCLE NNloop3 !** out of bounds
                       END IF
                       xtype = domain(xpos(1),xpos(2),xpos(3)) !** type of neighbor of the candidate hole
                       IF (xtype == 1) THEN  !** the neighbor of the candidate hole is an icepoint
                          xNNicecnt(m) = xNNicecnt(m) + 1  !** sum up all of the ice neighbors for each candidate holes
                          !** there's at least one ice neighbor to the hole under consideration
                          !** should we make a rule that the point up for moving will always
                          !** move to any ice-neighbored hole?  
                       ELSEIF (xtype == -1) THEN !** the neighbor is a liquid point
                          xNNliqcnt(m) = xNNliqcnt(m) + 1  !** sum up all of the liquid neighbors for each of the candidate holes
                       ELSEIF (xtype == 0 .OR. xtype == -99) THEN
                          xNNholcnt(m) = xNNholcnt(m) + 1 !** count the holes of NNs
                          !** don't really care about the hole neighbors of a hole, just here for completeness.
                       END IF
                    END DO NNloop3
                    
                 END IF
              END IF
              
           END DO NNloop2
!!$           IF (NNholcnt >= 26) THEN
!!$              PRINT *, 'orphan', n,m, cpos(1:3), NNicecnt, NNliqcnt, NNholcnt,MAXVAL(xNNicecnt(:)), maxval(xNNliqcnt(:))
!!$              !** special section to deal with orphans
!!$           END IF
           
           !** note: might need to randomly pick mx here, check for biases
           IF (MAXVAL(xNNicecnt(:)) > 0) THEN  !** of holes, were any of the holes' neighbors ice points?
              mx = MAXVAL(MAXLOC(xNNicecnt(:)))  !** find the hole with the most ice neighbors first
           ELSEIF (MAXVAL(xNNliqcnt(:)) > 0 .AND. MAXVAL(xNNicecnt(:)) == 0) THEN 
              mx = MAXVAL(MAXLOC(xNNliqcnt(:)))  !** in the absence of any holes having ice neighbors, pick the one with most liquid neighbors
           ELSE
              mx = 0
           END IF
           !** now we have the best target hole (mx) for the current point to move into 
           IF (mx > 0) THEN 
              tpos  = cpos + pdirs(mx,1:3) !** new position
              ttype = domain(tpos(1),tpos(2),tpos(3))  !** type of the candidate position
              IF (ttype /= 0) THEN
                 !** it should be zero (a hole), toss an error 
                 STOP 'ttype /= 0'
              END IF
              !** also update now empty hole after movement. 
              domain(tpos(1),tpos(2),tpos(3)) = -1  !** water domain moves into tpos
              domain(cpos(1),cpos(2),cpos(3)) =  0  !** empty cpos
              pos(n,1:3) = tpos(1:3)                !** update pos array
              who(tpos(1),tpos(2),tpos(3)) = n
              who(cpos(1),cpos(2),cpos(3)) = 0
              CALL mark_dirty(cpos)
              CALL mark_dirty(tpos)
              movecnt = movecnt + 1                 !** increment the move counter
           END IF
        END DO moveloop
        
        !        print '(2I5,2I7,3I6,G10.2)', niter, moveiter, icecnt, liqcnt, meltcnt, movecnt, NNicecntmax, fliq
        IF (movecnt < 20) EXIT moveiterloop
        
     END DO moveiterloop
     
     !** trigger output when specific melt criteria are met
     !** melted mass fraction: every original ice dipole carries the same mass, whether it is now
     !** a liquid dipole or was removed.  (Liquid volume fraction is nliq/(NDIP-ngone).)
     fliq = DBLE(nliq + ngone)/DBLE(NDIP)

     !** write a frame when the next target has been reached (to within half a dipole)
     IF (itarget <= ntarg) THEN
        IF (fliq >= targets(itarget) - 0.5d0/NDIP .AND. icecnt > 0) THEN
           PRINT '(A,2I9,2F12.5)', 'output1:',  niter, nliq, fliq, targets(itarget)
           DO WHILE (itarget <= ntarg)
              IF (fliq < targets(itarget) - 0.5d0/NDIP) EXIT
              itarget = itarget + 1
           END DO
           fliqold = fliq
           IF (OUTPUT == 1) CALL write_frame()
        END IF
     END IF
     

     !** final output
     IF (icecnt == 0 .AND. meltcnt == 0 .AND. movecnt == 0) THEN
        fliq = DBLE(nliq + ngone)/DBLE(NDIP)    !** melted mass fraction
        PRINT '(A,2I9,2G12.2)', 'output2:',  niter, nliq, fliq, 1.0
        PRINT '(A,I9,A,I9,A,F8.5)', 'mass conservation: dipoles removed ', ngone, ' of ', NDIP, &
             '  final/initial volume ', DBLE(NDIP - ngone)/DBLE(NDIP)
        IF (OUTPUT == 1) CALL write_frame()
        
        EXIT iterloop !** bust out of iteration loop once all points melted
     END IF
     !** need iterloop breaking criteria (i.e., when particle is completely melted
  END DO iterloop
  
CONTAINS

  SUBROUTINE seed_rng(sd)
    INTEGER, INTENT(in) :: sd
    INTEGER :: ns, is
    INTEGER, ALLOCATABLE :: sv(:)
    CALL RANDOM_SEED(size=ns)
    ALLOCATE(sv(ns))
    DO is = 1,ns
       sv(is) = sd + 37*is
    END DO
    CALL RANDOM_SEED(put=sv)
  END SUBROUTINE seed_rng

  !** number of ice dipoles that meet the melting criterion now
  INTEGER FUNCTION count_eligible()
    INTEGER :: nd, mb, tp(3), nice
    count_eligible = 0
    DO nd = 1,NDIP
       IF (pos(nd,1) < 0) CYCLE
       IF (domain(pos(nd,1),pos(nd,2),pos(nd,3)) /= 1) CYCLE
       nice = 0
       DO mb = 1,27
          IF (mb == 14) CYCLE
          tp = pos(nd,1:3) + PDIRS(mb,1:3)
          IF (ANY(tp(:) <= 0) .OR. tp(1) > XMAX .OR. tp(2) > YMAX .OR. tp(3) > ZMAX) CYCLE
          IF (domain(tp(1),tp(2),tp(3)) == 1) nice = nice + 1
       END DO
       IF (nice <= NNicecntmax) count_eligible = count_eligible + 1
    END DO
  END FUNCTION count_eligible

  !** the cell at c has changed: every dipole within two cells of it must be looked at again
  SUBROUTINE mark_dirty(c)
    INTEGER, INTENT(in) :: c(3)
    INTEGER :: a, b, d, w
    DO d = MAX(1,c(3)-2), MIN(ZMAX,c(3)+2)
       DO b = MAX(1,c(2)-2), MIN(YMAX,c(2)+2)
          DO a = MAX(1,c(1)-2), MIN(XMAX,c(1)+2)
             w = who(a,b,d)
             IF (w > 0) dirty(w) = .TRUE.
          END DO
       END DO
    END DO
  END SUBROUTINE mark_dirty

  !** write the current state for every dipole still present:
  !**   index, x, y, z, type (1 = ice, -1 = liquid), body number (1 = largest connected body)
  SUBROUTINE write_frame()
    INTEGER :: nd, nout, top(3)
    CALL label_bodies()
    WRITE(outfn,'(''domain5_f'',I6.6,''.out'')') NINT(fliq*100000.0d0)
    OPEN(OUTFID,file=outfn,status='replace')
    nout = 0
    DO nd = 1,NDIP
       IF (pos(nd,1) < 0) CYCLE  !** removed by mass conservation
       nout = nout + 1
       WRITE(OUTFID, '(6I9)') nout, pos(nd,1:3), domain(pos(nd,1),pos(nd,2),pos(nd,3)), body(nd)
    END DO
    CLOSE(OUTFID)
    top = 0                                     !** sizes of the three largest bodies (0 if fewer)
    DO nd = 1,MIN(3,nbody)
       top(nd) = bsize(border(nd))
    END DO
    PRINT '(A,I8,A,3I9,A,I9)', 'bodies: ', nbody, '  largest: ', top, '  detached dipoles: ', nout - top(1)
  END SUBROUTINE write_frame

  !** label connected bodies by flood fill.  While a body is being filled its cells in 'domain' hold
  !** 100 + 2*body + (1 if liquid), so no second domain-sized array is needed; they are restored to
  !** 1 / -1 before returning.  Sets nbody, bsize, border (body ids, largest first) and body(n).
  SUBROUTINE label_bodies()
    INTEGER :: nd, qh, qt, mb, v, cp(3), tp(3)
    nbody = 0
    DO nd = 1,NDIP
       IF (pos(nd,1) < 0) CYCLE
       v = domain(pos(nd,1),pos(nd,2),pos(nd,3))
       IF (v >= 100) CYCLE                      !** already part of a labelled body
       nbody = nbody + 1
       bsize(nbody) = 0
       qh = 1
       qt = 1
       cq(:,1) = pos(nd,1:3)
       domain(pos(nd,1),pos(nd,2),pos(nd,3)) = 100 + 2*nbody + MERGE(1,0,v == -1)
       DO WHILE (qh <= qt)
          cp = cq(:,qh)
          qh = qh + 1
          bsize(nbody) = bsize(nbody) + 1
          DO mb = 1,27
             IF (mb == 14) CYCLE                !** skip self
             tp = cp + PDIRS(mb,1:3)
             IF (ANY(tp(:) <= 0) .OR. tp(1) > XMAX .OR. tp(2) > YMAX .OR. tp(3) > ZMAX) CYCLE
             v = domain(tp(1),tp(2),tp(3))
             IF (v == 1 .OR. v == -1) THEN
                domain(tp(1),tp(2),tp(3)) = 100 + 2*nbody + MERGE(1,0,v == -1)
                qt = qt + 1
                cq(:,qt) = tp
             END IF
          END DO
       END DO
    END DO
    !** rank bodies by size, largest first
    DO nd = 1,nbody
       border(nd) = nd
    END DO
    IF (nbody > 1) CALL sort_bodies(1, nbody)
    DO nd = 1,nbody
       brank(border(nd)) = nd
    END DO
    !** per-dipole body number, and restore the domain
    DO nd = 1,NDIP
       IF (pos(nd,1) < 0) THEN
          body(nd) = 0
          CYCLE
       END IF
       v = domain(pos(nd,1),pos(nd,2),pos(nd,3)) - 100
       body(nd) = brank(v/2)
       domain(pos(nd,1),pos(nd,2),pos(nd,3)) = MERGE(-1,1,MOD(v,2) == 1)
    END DO
  END SUBROUTINE label_bodies

  !** quicksort of border(l:h) by descending body size
  RECURSIVE SUBROUTINE sort_bodies(l, h)
    INTEGER, INTENT(in) :: l, h
    INTEGER :: il, ih, pv, tmpb
    IF (l >= h) RETURN
    pv = bsize(border((l + h)/2))
    il = l
    ih = h
    DO
       DO WHILE (bsize(border(il)) > pv)
          il = il + 1
       END DO
       DO WHILE (bsize(border(ih)) < pv)
          ih = ih - 1
       END DO
       IF (il <= ih) THEN
          tmpb = border(il)
          border(il) = border(ih)
          border(ih) = tmpb
          il = il + 1
          ih = ih - 1
       END IF
       IF (il > ih) EXIT
    END DO
    CALL sort_bodies(l, ih)
    CALL sort_bodies(il, h)
  END SUBROUTINE sort_bodies
  
  INTEGER FUNCTION sgn(ns)
    INTEGER :: ns
    IF (ABS(ns) > 0) THEN
       sgn = ns/ABS(ns)
    ELSE
       sgn = 0
    END IF
  END FUNCTION sgn

  
END PROGRAM domain_melt4
