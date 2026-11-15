!** SPMM prototype: bond-count (pair-energy) liquid movement rule.
!**
!** Every site has a cohesion score  S(s) = sum over kernel offsets o of w(|o|) * a(type(s+o)),
!** with a = a_i for ice, a_w for water, 0 for empty.  A water point at s may move to an empty
!** site h when  S(h) - a_w*w(|h-s|) - S(s) > margin  (its own contribution to S(h) removed).
!** With a symmetric kernel each accepted move lowers the pair energy
!**   E = -1/2 * sum_n a(type_n) * S(pos_n),   so relaxation cannot oscillate.
!**
!** The melt rule is the one in domain_melt5.f90 (26-neighbour ice count threshold).
!** Input: shape file with nheader header lines, then rows "i j k m"; m == 2 marks water, else ice.
program spmm_bond
  implicit none
  integer, parameter :: dp = kind(1.0d0), i1 = selected_int_kind(2)

  !** run-time parameters (namelist file spmm.nml, optional)
  character(len=256) :: shape_file = 'shape.dat', prefix = 'bond'
  real(dp) :: a_i = 1.0_dp, a_w = 1.0_dp, margin = 0.0_dp, rkern = 3.0_dp
  real(dp) :: p_melt = 1.0_dp, move_frac = 0.2_dp, dfmax = 0.005_dp, temp = 0.0_dp, temp_final = -1.0_dp
  !** frames are written at the melted fractions in targets (ascending).  Each melt step melts just
  !** enough to reach the next target, at most dfmax of the particle (dfmax_late from f_late onwards).
  !** These can also be set in group &melt of spmm.nml, which spmm reads too.
  integer, parameter :: maxtarg = 100
  real(dp) :: targets(maxtarg) = -1.0_dp, dfmax_late = 0.02_dp, f_late = 0.20_dp, nextt
  integer  :: ntarg, need, left
  logical  :: do_melt = .true., pin_coating = .false., nonlocal = .true., check = .false.
  logical  :: reattach = .true.    !** move bodies that lose contact with the rest back into contact
  integer  :: max_sweeps = 40, max_final_sweeps = 20000, nheader = 1, seed = 12345, sweeps_hot = 0, final_hot = 0
  namelist /spmm/ shape_file, prefix, a_i, a_w, margin, rkern, p_melt, move_frac, dfmax, temp, temp_final, sweeps_hot, final_hot, &
       do_melt, pin_coating, nonlocal, check, max_sweeps, max_final_sweeps, nheader, seed, reattach, &
       targets, dfmax_late, f_late
  namelist /melt/ targets, dfmax, dfmax_late, f_late, seed

  real(dp), parameter :: eps = 1.0e-9_dp, pi = 3.141592653589793_dp
  integer(i1), allocatable :: typ(:,:,:)
  real(dp), allocatable :: S(:,:,:), Sref(:,:,:), kw(:), wgrid(:,:,:), wkey(:), hkey(:)
  integer, allocatable :: idx(:,:,:), mark(:,:,:), pos(:,:), comp(:), queue(:), koff(:,:)
  integer, allocatable :: wl(:), hpos(:,:), hcomp(:), worder(:), horder(:), cstw(:), csth(:), cnt(:)
  integer, allocatable :: bl(:), bq(:), bst(:)      !** body label, members body by body, start of each body
  integer(i1), allocatable :: btyp(:)
  !** nocc: occupied cells among the 26 neighbours of each cell, so interior dipoles are skipped cheaply.
  !** Water bodies are labelled over runs of water cells along x: rid is the run in each water cell,
  !** rend its last x, rpar the union-find parent, rnew the body number of a root run.
  integer(i1), allocatable :: nocc(:,:,:)
  integer, allocatable :: rid(:,:,:), rend(:), rpar(:), rnew(:), rpos(:,:)
  integer :: nre, nst, totre, totst
  integer :: pd(3,26), dims(3), lo(3), hi(3), mn(3), mx(3), off(3), p(3), q(3)
  integer :: N, nk, R, i, j, k, m, n1, ios, fid, hcap, sweepid, ncomp, nws, nh
  integer :: niter, nmov, nsw, totsw, totmov, nice, nliq, meltcnt, nnmax, nnice, itarget
  integer(8) :: c0, c1, crate
  integer :: nelig
  real(dp) :: rc, rr, aeq, fliq, rx, E, Eold
  logical :: ex

  fid = 13
  inquire(file='spmm.nml', exist=ex)
  if (ex) then
     open(fid, file='spmm.nml', status='old'); read(fid, nml=spmm, iostat=ios); close(fid)
     if (ios > 0) stop 'error reading group &spmm in spmm.nml'
     open(fid, file='spmm.nml', status='old'); read(fid, nml=melt, iostat=ios); close(fid)
     if (ios > 0) stop 'error reading group &melt in spmm.nml'
  end if
  call seed_rng(seed)
  call system_clock(c0, crate)

  n1 = 0
  do i = -1, 1; do j = -1, 1; do k = -1, 1
     if (i == 0 .and. j == 0 .and. k == 0) cycle
     n1 = n1 + 1; pd(:, n1) = (/i, j, k/)
  end do; end do; end do

  !** kernel: all offsets with 0 < r <= rkern, smooth compact weight
  R = ceiling(rkern); rc = rkern + 1.0_dp
  allocate(wgrid(-R:R, -R:R, -R:R)); wgrid = 0.0_dp
  nk = 0
  do i = -R, R; do j = -R, R; do k = -R, R
     rr = sqrt(real(i*i + j*j + k*k, dp))
     if (rr > 0.0_dp .and. rr <= rkern + 1.0e-9_dp) then
        nk = nk + 1; wgrid(i, j, k) = (1.0_dp - (rr/rc)**2)**2
     end if
  end do; end do; end do
  allocate(koff(3, nk), kw(nk)); nk = 0
  do i = -R, R; do j = -R, R; do k = -R, R
     if (wgrid(i, j, k) > 0.0_dp) then
        nk = nk + 1; koff(:, nk) = (/i, j, k/); kw(nk) = wgrid(i, j, k)
     end if
  end do; end do; end do

  print '(A,F8.3,A,F8.3,A,F8.3)', 'kernel weights: half-space A=', sum(kw, mask=koff(3,:) < 0), &
       ' in-plane P=', sum(kw, mask=koff(3,:) == 0), ' total=', sum(kw)
  !** pass 1: count and extents
  open(fid, file=trim(shape_file), status='old')
  do i = 1, nheader; read(fid, *); end do
  N = 0; mn = huge(1); mx = -huge(1)
  do
     read(fid, *, iostat=ios) p, m
     if (ios /= 0) exit
     N = N + 1; mn = min(mn, p); mx = max(mx, p)
  end do
  close(fid)
  if (N == 0) stop 'no dipoles read'

  !** centred domain with room for the equal-volume sphere and the kernel halo
  aeq = (3.0_dp*N/(4.0_dp*pi))**(1.0_dp/3.0_dp)
  do i = 1, 3
     dims(i) = max(mx(i) - mn(i) + 1 + ceiling(aeq), ceiling(2.0_dp*aeq) + 2) + 2*(R + 3)
     off(i)  = (dims(i) - (mx(i) - mn(i) + 1))/2 - mn(i) + 1
  end do
  lo = R + 2; hi = dims - R - 1
  allocate(typ(dims(1), dims(2), dims(3)), S(dims(1), dims(2), dims(3)))
  allocate(idx(dims(1), dims(2), dims(3)), mark(dims(1), dims(2), dims(3)))
  allocate(nocc(dims(1), dims(2), dims(3)), rid(dims(1), dims(2), dims(3)), rend(N), rpar(N), rnew(N), rpos(3, N))
  nocc = 0_i1; rid = 0
  allocate(pos(3, N), comp(N), queue(N), wl(N), wkey(N), worder(N), cstw(N + 2), csth(N + 2), cnt(N + 1))
  hcap = 8*N + 1000
  allocate(hpos(3, hcap), hcomp(hcap), hkey(hcap), horder(hcap))
  typ = 0_i1; S = 0.0_dp; idx = 0; mark = 0; sweepid = 0
  if (reattach) allocate(bl(N), bq(N), bst(N + 1), btyp(N))
  totre = 0; totst = 0

  !** pass 2: load
  open(fid, file=trim(shape_file), status='old')
  do i = 1, nheader; read(fid, *); end do
  do n1 = 1, N
     read(fid, *) p, m
     p = p + off
     if (idx(p(1), p(2), p(3)) /= 0) stop 'duplicate dipole in shape file'
     pos(:, n1) = p; idx(p(1), p(2), p(3)) = n1
     if (m == 2) then
        typ(p(1), p(2), p(3)) = -1_i1
     else
        typ(p(1), p(2), p(3)) = 1_i1
     end if
  end do
  close(fid)
  do n1 = 1, N
     call add_occ(pos(:, n1), 1_i1)
  end do
  call recompute(S)
  print '(A,I9,A,3I5,A,I4,A,F6.2)', 'N=', N, ' domain=', dims, ' kernel sites=', nk, ' aeq=', aeq
  print '(A,3F7.3,A,2L2,A,L2)', 'a_i a_w margin=', a_i, a_w, margin, ' pin nonlocal=', pin_coating, nonlocal, ' melt=', do_melt

  if (targets(1) < 0.0_dp) targets(1:24) = (/0.005_dp, 0.01_dp, 0.02_dp, 0.03_dp, 0.04_dp, 0.05_dp, 0.06_dp, 0.07_dp, &
       0.08_dp, 0.09_dp, 0.10_dp, 0.125_dp, 0.15_dp, 0.175_dp, 0.20_dp, 0.25_dp, 0.30_dp, 0.40_dp, 0.50_dp, 0.60_dp, &
       0.70_dp, 0.80_dp, 0.90_dp, 0.95_dp/)
  ntarg = count(targets > 0.0_dp .and. targets < 1.0_dp)
  do i = 2, ntarg
     if (targets(i) <= targets(i - 1)) stop 'targets must be ascending, between 0 and 1'
  end do
  itarget = 1
  totsw = 0; totmov = 0; niter = 0; nnmax = 5
  call counts(nice, nliq); fliq = real(nliq, dp)/N
  if (do_melt) then
     call output_frac(fliq)
  else
     call output_name(trim(prefix)//'_init.out')
  end if

  if (do_melt) then
     meltloop: do
        niter = niter + 1; meltcnt = 0
        !** cap the melt per step: count currently eligible ice points, then thin the draw
        nelig = 0
        do n1 = 1, N
           p = pos(:, n1)
           if (typ(p(1), p(2), p(3)) /= 1_i1) cycle
           nnice = 0
           do m = 1, 26
              q = p + pd(:, m)
              if (typ(q(1), q(2), q(3)) == 1_i1) nnice = nnice + 1
           end do
           if (nnice <= nnmax) nelig = nelig + 1
        end do
        !** how many to melt this step: enough to reach the next target, at most the cap
        nextt = 1.0_dp
        if (itarget <= ntarg) nextt = targets(itarget)
        need = max(1, nint(merge(dfmax, dfmax_late, fliq < f_late)*N))
        need = max(1, min(need, ceiling((nextt - fliq)*N - 1.0e-6_dp)))
        if (p_melt < 1.0_dp) need = max(1, min(need, nint(p_melt*nelig)))
        left = nelig
        do n1 = 1, N
           p = pos(:, n1)
           if (typ(p(1), p(2), p(3)) /= 1_i1) cycle
           nnice = 0
           do m = 1, 26
              q = p + pd(:, m)
              if (typ(q(1), q(2), q(3)) == 1_i1) nnice = nnice + 1
           end do
           if (nnice <= nnmax) then
              !** pick 'need' of the 'left' eligible points at random, each equally likely
              call random_number(rx); left = left - 1
              if (need > 0 .and. rx*max(left + 1, 1) < need) then
                 need = need - 1
                 typ(p(1), p(2), p(3)) = -1_i1
                 call add_field(p, a_w - a_i)
                 meltcnt = meltcnt + 1
              end if
           end if
        end do
        if (nelig == 0) nnmax = nnmax + 1
        call relax(max_sweeps, sweeps_hot)
        if (reattach) then
           call reattach_bodies(nre, nst); totre = totre + nre; totst = totst + nst
        end if
        call counts(nice, nliq); fliq = real(nliq, dp)/N
        if (nice == 0) exit meltloop
        if (itarget <= ntarg) then
           if (fliq >= targets(itarget) - 0.5_dp/N) then
              call output_frac(fliq)
              do while (itarget <= ntarg)
                 if (fliq < targets(itarget) - 0.5_dp/N) exit
                 itarget = itarget + 1
              end do
           end if
        end if
     end do meltloop
  end if
  if (temp_final >= 0.0_dp) temp = temp_final   !** separate noise level for the final rounding
  if (reattach) then
     call reattach_bodies(nre, nst); totre = totre + nre; totst = totst + nst
  end if
  call relax(max_final_sweeps, final_hot)
  if (reattach) then      !** pieces that separated during the final rounding: bring back, then quench again
     do i = 1, 20
        call reattach_bodies(nre, nst); totre = totre + nre; totst = totst + nst
        if (nre == 0) exit
        call relax(max_final_sweeps, 0)
     end do
     print '(A,I8,A,I8)', 'reattach: body shifts=', totre, ' blocked=', totst
  end if
  if (do_melt) then
     call output_frac(1.0_dp)
  else
     call output_name(trim(prefix)//'_final.out')
  end if
  call system_clock(c1)
  print '(A,I7,A,I8,A,I10,A,ES15.7,A,F8.1,A,L2)', 'DONE melt_iters=', niter, ' sweeps=', totsw, ' moves=', totmov, &
       ' E=', energy(), ' seconds=', real(c1 - c0, dp)/crate, ' converged=', (nmov == 0)

contains

  subroutine seed_rng(sd)
    integer, intent(in) :: sd
    integer :: ns, ii
    integer, allocatable :: sv(:)
    call random_seed(size=ns); allocate(sv(ns))
    do ii = 1, ns; sv(ii) = sd + 37*ii; end do
    call random_seed(put=sv)
  end subroutine seed_rng

  pure real(dp) function aval(t)
    integer(i1), intent(in) :: t
    if (t == 1_i1) then
       aval = a_i
    else if (t == -1_i1) then
       aval = a_w
    else
       aval = 0.0_dp
    end if
  end function aval

  subroutine add_occ(pp, d)
    integer, intent(in) :: pp(3)
    integer(i1), intent(in) :: d
    integer :: kk
    do kk = 1, 26
       nocc(pp(1) + pd(1, kk), pp(2) + pd(2, kk), pp(3) + pd(3, kk)) = &
            nocc(pp(1) + pd(1, kk), pp(2) + pd(2, kk), pp(3) + pd(3, kk)) + d
    end do
  end subroutine add_occ

  integer function rroot(r0) result(r)
    integer, intent(in) :: r0
    integer :: t, nx
    r = r0
    do while (rpar(r) /= r); r = rpar(r); end do
    t = r0
    do while (rpar(t) /= r)
       nx = rpar(t); rpar(t) = r; t = nx
    end do
  end function rroot

  subroutine add_field(pp, delta)
    integer, intent(in) :: pp(3)
    real(dp), intent(in) :: delta
    integer :: kk
    do kk = 1, nk
       S(pp(1) + koff(1, kk), pp(2) + koff(2, kk), pp(3) + koff(3, kk)) = &
            S(pp(1) + koff(1, kk), pp(2) + koff(2, kk), pp(3) + koff(3, kk)) + delta*kw(kk)
    end do
  end subroutine add_field

  subroutine recompute(F)
    real(dp), intent(out) :: F(:, :, :)
    integer :: nn, kk, pp(3)
    real(dp) :: av
    F = 0.0_dp
    do nn = 1, N
       pp = pos(:, nn); av = aval(typ(pp(1), pp(2), pp(3)))
       do kk = 1, nk
          F(pp(1) + koff(1, kk), pp(2) + koff(2, kk), pp(3) + koff(3, kk)) = &
               F(pp(1) + koff(1, kk), pp(2) + koff(2, kk), pp(3) + koff(3, kk)) + av*kw(kk)
       end do
    end do
  end subroutine recompute

  real(dp) function energy()
    integer :: nn, pp(3)
    energy = 0.0_dp
    do nn = 1, N
       pp = pos(:, nn)
       energy = energy - 0.5_dp*aval(typ(pp(1), pp(2), pp(3)))*S(pp(1), pp(2), pp(3))
    end do
  end function energy

  subroutine counts(ni, nl)
    integer, intent(out) :: ni, nl
    integer :: nn, pp(3)
    ni = 0; nl = 0
    do nn = 1, N
       pp = pos(:, nn)
       if (typ(pp(1), pp(2), pp(3)) == 1_i1) ni = ni + 1
       if (typ(pp(1), pp(2), pp(3)) == -1_i1) nl = nl + 1
    end do
  end subroutine counts

  subroutine invariants(tag)
    character(len=*), intent(in) :: tag
    integer :: nn, pp(3)
    do nn = 1, N
       pp = pos(:, nn)
       if (idx(pp(1), pp(2), pp(3)) /= nn) then
          print *, 'INVARIANT FAIL (index grid) ', tag, nn; stop 1
       end if
       if (typ(pp(1), pp(2), pp(3)) == 0_i1) then
          print *, 'INVARIANT FAIL (empty dipole) ', tag, nn; stop 1
       end if
    end do
    do nn = 1, N
       pp = pos(:, nn)
       if (nocc(pp(1), pp(2), pp(3)) /= count(typ(pp(1)-1:pp(1)+1, pp(2)-1:pp(2)+1, pp(3)-1:pp(3)+1) /= 0_i1) - 1) then
          print *, 'INVARIANT FAIL (neighbour count) ', tag, nn; stop 1
       end if
    end do
    if (count(typ /= 0_i1) /= N .or. count(idx /= 0) /= N) then
       print *, 'INVARIANT FAIL (count) ', tag, count(typ /= 0_i1), count(idx /= 0), N; stop 1
    end if
  end subroutine invariants

  subroutine relax(maxsw, nhot)
    integer, intent(in) :: maxsw, nhot
    integer :: isw
    real(dp) :: dmax
    nmov = 0
    if (temp > 0.0_dp) then      !** noisy phase: lets flat lattice faces nucleate new layers
       do isw = 1, nhot
          nmov = relax_sweep(temp); totsw = totsw + 1; totmov = totmov + nmov
       end do
       if (check) call invariants('hot')
    end if
    if (check) Eold = energy()
    do isw = 1, maxsw
       nmov = relax_sweep(0.0_dp)
       totsw = totsw + 1; totmov = totmov + nmov
       if (check) then
          E = energy()
          if (E > Eold + 1.0e-6_dp) then
             print *, 'INVARIANT FAIL (energy increased) ', Eold, E; stop 1
          end if
          Eold = E
          call invariants('sweep')
          if (mod(totsw, 50) == 0) then
             if (.not. allocated(Sref)) allocate(Sref(dims(1), dims(2), dims(3)))
             call recompute(Sref); dmax = maxval(abs(Sref - S))
             if (dmax > 1.0e-6_dp) then
                print *, 'INVARIANT FAIL (score field drift) ', dmax; stop 1
             end if
          end if
       end if
       if (nmov == 0) exit
    end do
  end subroutine relax

  !** Bodies (ice and water together, 26-connectivity) that have lost contact with everything else
  !** are shifted rigidly, without changing shape, to the nearest position that touches another
  !** body.  The largest body never moves.  nmoved counts shifts made, nstuck bodies that could
  !** not be shifted this time.
  subroutine reattach_bodies(nmoved, nstuck)
    integer, intent(out) :: nmoved, nstuck
    integer :: nn, c, mm, nb, qh, qt, nbody, main, b, pp(3), qq(3), v(3), vt(3), s, L, kk, ii, i0, i1e, try
    integer :: nfar, npart, pass, nshift, njoin
    logical :: touching
    real(dp) :: cm(3), cb(3)
    logical :: hit, ok
    nmoved = 0; nstuck = 0; nfar = 0; npart = 0
    !** a body may land on another detached body; relabel and go again until all are one
    passes: do pass = 1, 200
    nshift = 0; njoin = 0; nstuck = 0
    bl = 0; nbody = 0; qt = 0
    do nn = 1, N
       if (bl(nn) /= 0) cycle
       nbody = nbody + 1; bl(nn) = nbody; qt = qt + 1; bq(qt) = nn; qh = qt; bst(nbody) = qt
       do while (qh <= qt)
          c = bq(qh); qh = qh + 1
          do mm = 1, 26
             qq = pos(:, c) + pd(:, mm); nb = idx(qq(1), qq(2), qq(3))
             if (nb /= 0) then
                if (bl(nb) == 0) then
                   bl(nb) = nbody; qt = qt + 1; bq(qt) = nb
                end if
             end if
          end do
       end do
    end do
    bst(nbody + 1) = N + 1
    if (nbody == 1) exit passes
    main = 1
    do b = 2, nbody
       if (bst(b + 1) - bst(b) > bst(main + 1) - bst(main)) main = b
    end do
    cm = 0.0_dp
    do ii = bst(main), bst(main + 1) - 1; cm = cm + pos(:, bq(ii)); end do
    cm = cm/(bst(main + 1) - bst(main))

    do b = 1, nbody
       if (b == main) cycle
       i0 = bst(b); i1e = bst(b + 1) - 1
       touching = .false.       !** another body was set down against this one earlier in this pass
       chk: do ii = i0, i1e
          do mm = 1, 26
             qq = pos(:, bq(ii)) + pd(:, mm); nb = idx(qq(1), qq(2), qq(3))
             if (nb /= 0) then
                if (bl(nb) /= b) then
                   touching = .true.; exit chk
                end if
             end if
          end do
       end do chk
       if (touching) then
          njoin = njoin + 1; cycle
       end if
       tries: do try = 1, 400     !** a shift can be cut short by an obstacle, so repeat until the body touches
       !** breadth-first search outward through empty sites for the nearest site that touches another body;
       !** each site remembers the dipole of this body it started from, which gives the shift
       sweepid = sweepid + 1; qt = 0; qh = 1; hit = .false.; v = 0
       do ii = i0, i1e
          do mm = 1, 26
             qq = pos(:, bq(ii)) + pd(:, mm)
             if (typ(qq(1), qq(2), qq(3)) /= 0_i1 .or. any(qq < lo) .or. any(qq > hi)) cycle
             if (mark(qq(1), qq(2), qq(3)) == sweepid .or. qt >= hcap) cycle
             mark(qq(1), qq(2), qq(3)) = sweepid; qt = qt + 1; hpos(:, qt) = qq; hcomp(qt) = bq(ii)
          end do
       end do
       search: do while (qh <= qt)
          pp = hpos(:, qh); s = hcomp(qh); qh = qh + 1
          do mm = 1, 26
             qq = pp + pd(:, mm); nb = idx(qq(1), qq(2), qq(3))
             if (nb /= 0) then
                if (bl(nb) /= b) then
                   hit = .true.; v = pp - pos(:, s); exit search
                end if
             end if
          end do
          do mm = 1, 26
             qq = pp + pd(:, mm)
             if (typ(qq(1), qq(2), qq(3)) /= 0_i1 .or. any(qq < lo) .or. any(qq > hi)) cycle
             if (mark(qq(1), qq(2), qq(3)) == sweepid .or. qt >= hcap) cycle
             mark(qq(1), qq(2), qq(3)) = sweepid; qt = qt + 1; hpos(:, qt) = qq; hcomp(qt) = s
          end do
       end do search
       if (.not. hit) then      !** nothing within reach: one step toward the centre of the largest body
          cb = 0.0_dp
          do ii = i0, i1e; cb = cb + pos(:, bq(ii)); end do
          cb = cm - cb/(i1e - i0 + 1)
          if (maxval(abs(cb)) > 0.0_dp) v = nint(min(8.0_dp, maxval(abs(cb)))*cb/maxval(abs(cb)))
          nfar = nfar + 1
       end if
       !** longest part of the shift that keeps every dipole inside the domain and off other bodies
       L = maxval(abs(v)); ok = .false.; vt = 0
       do kk = L, 1, -1
          vt = nint(real(v, dp)*kk/L)
          if (all(vt == 0)) cycle
          ok = .true.
          do ii = i0, i1e
             qq = pos(:, bq(ii)) + vt
             if (any(qq < lo) .or. any(qq > hi)) then
                ok = .false.; exit
             end if
             nb = idx(qq(1), qq(2), qq(3))
             if (nb /= 0) then
                if (bl(nb) /= b) then
                   ok = .false.; exit
                end if
             end if
          end do
          if (ok) exit
       end do
       if (.not. ok) then
          nstuck = nstuck + 1; exit tries
       end if
       do ii = i0, i1e            !** lift the whole body out, then set it down shifted
          nn = bq(ii); pp = pos(:, nn); btyp(nn) = typ(pp(1), pp(2), pp(3))
          call add_field(pp, -aval(btyp(nn))); call add_occ(pp, -1_i1)
          typ(pp(1), pp(2), pp(3)) = 0_i1; idx(pp(1), pp(2), pp(3)) = 0
       end do
       do ii = i0, i1e
          nn = bq(ii); pp = pos(:, nn) + vt
          typ(pp(1), pp(2), pp(3)) = btyp(nn); idx(pp(1), pp(2), pp(3)) = nn; pos(:, nn) = pp
          call add_field(pp, aval(btyp(nn))); call add_occ(pp, 1_i1)
       end do
       nmoved = nmoved + 1; nshift = nshift + 1
       if (hit .and. all(vt == v)) exit tries      !** full shift made: the body now touches another
       npart = npart + 1
       end do tries
    end do
    if (check) print '(A,I4,A,I6,A,I6,A,I6,A,I6,A,I6)', 'reattach pass', pass, ': bodies=', nbody, ' shifts=', nshift, &
         ' joined=', njoin, ' out of search range=', nfar, ' blocked=', nstuck
    if (nshift == 0 .and. njoin == 0) exit passes
    end do passes
    if (check) call invariants('reattach')
  end subroutine reattach_bodies

  subroutine move(nn, h)
    integer, intent(in) :: nn, h(3)
    integer :: ss(3)
    ss = pos(:, nn)
    call add_field(ss, -a_w); call add_occ(ss, -1_i1)
    typ(ss(1), ss(2), ss(3)) = 0_i1; idx(ss(1), ss(2), ss(3)) = 0
    typ(h(1), h(2), h(3)) = -1_i1; idx(h(1), h(2), h(3)) = nn
    pos(:, nn) = h
    call add_field(h, a_w); call add_occ(h, 1_i1)
  end subroutine move

  !** score gain for water point nn moving to empty site h (own contribution removed)
  real(dp) function gain(nn, h)
    integer, intent(in) :: nn, h(3)
    integer :: ss(3), d(3)
    ss = pos(:, nn); d = h - ss
    gain = S(h(1), h(2), h(3)) - S(ss(1), ss(2), ss(3))
    if (all(abs(d) <= R)) gain = gain - a_w*wgrid(d(1), d(2), d(3))
  end function gain

  integer function relax_sweep(tk) result(nm)
    real(dp), intent(in) :: tk
    integer :: nn, c, mm, nb, pp(3), qq(3), ii, jj, w0, w1, h0, h1, npair, best(3)
    logical :: mobile, hashole
    integer :: nrun
    integer, parameter :: fwd(2, 4) = reshape((/1, -1, 1, 0, 1, 1, 0, 1/), (/2, 4/))
    real(dp) :: g, gbest, r1, r2
    integer(i1) :: t
    nm = 0
    !** 1. label connected water bodies (26-connectivity), numbered by their lowest dipole index.
    !**    Runs of water cells along x are found, joined to the runs they touch in the four
    !**    neighbouring columns ahead of them, and the joined sets are the bodies.
    nrun = 0
    do nn = 1, N
       pp = pos(:, nn)
       if (typ(pp(1), pp(2), pp(3)) /= -1_i1) cycle
       if (typ(pp(1) - 1, pp(2), pp(3)) == -1_i1) cycle          !** not the start of a run
       nrun = nrun + 1; rpar(nrun) = nrun; rnew(nrun) = 0; rpos(:, nrun) = pp
       ii = pp(1)
       do while (typ(ii, pp(2), pp(3)) == -1_i1)
          rid(ii, pp(2), pp(3)) = nrun; ii = ii + 1
       end do
       rend(nrun) = ii - 1
    end do
    do c = 1, nrun
       do mm = 1, 4
          qq(2) = rpos(2, c) + fwd(1, mm); qq(3) = rpos(3, c) + fwd(2, mm)
          ii = rpos(1, c) - 1
          do while (ii <= rend(c) + 1)
             if (typ(ii, qq(2), qq(3)) == -1_i1) then
                nb = rid(ii, qq(2), qq(3))
                w0 = rroot(c); w1 = rroot(nb)
                if (w0 /= w1) rpar(max(w0, w1)) = min(w0, w1)
                ii = rend(nb) + 1
             else
                ii = ii + 1
             end if
          end do
       end do
    end do
    ncomp = 0
    do nn = 1, N
       pp = pos(:, nn)
       if (typ(pp(1), pp(2), pp(3)) /= -1_i1) then
          comp(nn) = 0; cycle
       end if
       c = rroot(rid(pp(1), pp(2), pp(3)))
       if (rnew(c) == 0) then
          ncomp = ncomp + 1; rnew(c) = ncomp
       end if
       comp(nn) = rnew(c)
    end do
    if (ncomp == 0) return
    !** 2. mobile surface water points, and empty sites bordering water
    sweepid = sweepid + 1; nws = 0; nh = 0
    do nn = 1, N
       pp = pos(:, nn)
       if (typ(pp(1), pp(2), pp(3)) /= -1_i1) cycle
       if (nocc(pp(1), pp(2), pp(3)) == 26_i1) cycle             !** interior: no empty site beside it
       mobile = .true.; hashole = .false.
       do mm = 1, 26
          qq = pp + pd(:, mm); t = typ(qq(1), qq(2), qq(3))
          if (t == 1_i1 .and. pin_coating) mobile = .false.
          if (t == 0_i1) then
             hashole = .true.
             if (all(qq >= lo) .and. all(qq <= hi) .and. nh < hcap) then
                if (mark(qq(1), qq(2), qq(3)) /= sweepid) then
                   mark(qq(1), qq(2), qq(3)) = sweepid
                   nh = nh + 1; hpos(:, nh) = qq; hcomp(nh) = comp(nn)
                end if
             end if
          end if
       end do
       if (mobile .and. hashole) then
          nws = nws + 1; wl(nws) = nn
       end if
    end do
    if (nws == 0 .or. nh == 0) return
    do ii = 1, nws
       pp = pos(:, wl(ii)); call random_number(r1)
       wkey(ii) = S(pp(1), pp(2), pp(3)) + 1.0e-7_dp*r1 - tk*log(max(r1, 1.0e-12_dp))   !** ascending: least bound first
    end do

    if (.not. nonlocal) then
       !** local mode: each point may step to its best adjacent empty site
       do ii = 1, nws; worder(ii) = ii; end do
       call qsort(worder, 1, nws, wkey)
       do ii = 1, nws
          nn = wl(worder(ii)); pp = pos(:, nn); gbest = margin + eps; best = 0
          do mm = 1, 26
             qq = pp + pd(:, mm)
             if (typ(qq(1), qq(2), qq(3)) /= 0_i1) cycle
             if (any(qq < lo) .or. any(qq > hi)) cycle
             g = gain(nn, qq)
             if (g > gbest) then
                gbest = g; best = qq
             end if
          end do
          if (best(1) /= 0) then
             call move(nn, best); nm = nm + 1
          end if
       end do
       return
    end if

    !** 3. non-local mode: within each water body pair the least-bound points with the best empty sites
    do ii = 1, nh
       call random_number(r1)
       hkey(ii) = -S(hpos(1, ii), hpos(2, ii), hpos(3, ii)) + 1.0e-7_dp*r1 - tk*log(max(r1, 1.0e-12_dp))   !** ascending: best site first
    end do
    cnt(1:ncomp + 1) = 0
    do ii = 1, nws; c = comp(wl(ii)); cnt(c) = cnt(c) + 1; end do
    cstw(1) = 1
    do c = 1, ncomp; cstw(c + 1) = cstw(c) + cnt(c); end do
    cnt(1:ncomp) = 0
    do ii = 1, nws
       c = comp(wl(ii)); worder(cstw(c) + cnt(c)) = ii; cnt(c) = cnt(c) + 1
    end do
    cnt(1:ncomp + 1) = 0
    do ii = 1, nh; c = hcomp(ii); cnt(c) = cnt(c) + 1; end do
    csth(1) = 1
    do c = 1, ncomp; csth(c + 1) = csth(c) + cnt(c); end do
    cnt(1:ncomp) = 0
    do ii = 1, nh
       c = hcomp(ii); horder(csth(c) + cnt(c)) = ii; cnt(c) = cnt(c) + 1
    end do
    do c = 1, ncomp
       w0 = cstw(c); w1 = cstw(c + 1) - 1; h0 = csth(c); h1 = csth(c + 1) - 1
       if (w1 < w0 .or. h1 < h0) cycle
       call qsort(worder, w0, w1, wkey)
       call qsort(horder, h0, h1, hkey)
       npair = max(1, int(move_frac*min(w1 - w0 + 1, h1 - h0 + 1)))
       do jj = 0, npair - 1
          nn = wl(worder(w0 + jj)); qq = hpos(:, horder(h0 + jj))
          if (typ(qq(1), qq(2), qq(3)) /= 0_i1) cycle
          g = gain(nn, qq)
          if (g > margin + eps) then
             call move(nn, qq); nm = nm + 1
          else if (tk > 0.0_dp) then
             call random_number(r2)
             if (r2 < exp((g - margin)/tk)) then
                call move(nn, qq); nm = nm + 1
             end if
          else
             exit
          end if
       end do
    end do
  end function relax_sweep

  recursive subroutine qsort(a, l, h, key)
    integer, intent(inout) :: a(:)
    integer, intent(in) :: l, h
    real(dp), intent(in) :: key(:)
    integer :: ii, jj, tmp
    real(dp) :: pv
    if (l >= h) return
    pv = key(a((l + h)/2)); ii = l; jj = h
    do
       do while (key(a(ii)) < pv); ii = ii + 1; end do
       do while (key(a(jj)) > pv); jj = jj - 1; end do
       if (ii <= jj) then
          tmp = a(ii); a(ii) = a(jj); a(jj) = tmp; ii = ii + 1; jj = jj - 1
       end if
       if (ii > jj) exit
    end do
    call qsort(a, l, jj, key); call qsort(a, ii, h, key)
  end subroutine qsort

  subroutine output_frac(f)
    real(dp), intent(in) :: f
    character(len=300) :: fn
    write(fn, '(A,''_f'',I6.6,''.out'')') trim(prefix), nint(f*100000.0_dp)
    call output_name(trim(fn))
    print '(A,I7,A,F8.5,A,I8,A,I10)', 'frame iter=', niter, ' fliq=', f, ' sweeps=', totsw, ' moves=', totmov
  end subroutine output_frac

  subroutine output_name(fn)
    character(len=*), intent(in) :: fn
    integer :: nn, pp(3)
    open(17, file=fn, status='replace')
    do nn = 1, N
       pp = pos(:, nn)
       write(17, '(5I9)') nn, pp - off, int(typ(pp(1), pp(2), pp(3)))
    end do
    close(17)
  end subroutine output_name

end program spmm_bond
