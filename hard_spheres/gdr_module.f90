MODULE gdr_module

    CONTAINS

  ! Compute the radial distribution function g(r) from positions array 
  ! with (ndim, N) layout and filling the arrays bins and hist as a result.
  ! The two arrays bins and hist must be already allocated by the client code.
  ! Bins are defined as i*dr + dr/2 with i=1,...,size(bins) and filled here.
  ! The histogram is reset at each call and is not normalized at the end.
  subroutine gr_self(positions, box, dr, hist, bins)
    real(8), intent(in)       :: positions(:,:), dr
    integer(8), intent(inout) :: hist(:)
    real(8), intent(inout)    :: bins(:)
    real(8), intent(in)       :: box(:)
    real(8)                   :: distances(size(positions,2))  ! stack
    real(8)    :: dist(size(box)), dist_sq, pos(size(box)), hbox(size(box)), rdist, rmax
    integer(8) :: i, j, ii, bin, k
    rmax = dr * size(bins)
    hbox = box / 2
    hist = 0
    do i = 1, size(positions,2)
       pos = positions(:,i)
       ! Compute distances with particle i
       k = 0
       do j=i+1,size(positions,2)
          dist(:) = positions(:,j) - pos(:)
          where (abs(dist) > hbox)
             dist = dist - sign(box,dist)
          end where
          rdist = sqrt(sum(dist**2))
          if (rdist < rmax) then
             k = k+1
             distances(k) = rdist
          end if
       end do
       ! Bin distances
       do j=1,k
          bin = floor(distances(j) / dr) + 1
          hist(bin) = hist(bin) + 1
       end do
    end do
    ! Bins
    do j=1,size(bins)
       bins(j) = dr * (j-1) + dr / 2
    end do
  end subroutine gr_self

    
END MODULE gdr_module