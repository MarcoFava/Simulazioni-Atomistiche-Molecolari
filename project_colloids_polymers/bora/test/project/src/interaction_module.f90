module interaction_module
use potential_module

implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine pbc(r, box)
        double precision, intent(inout) :: r(:)
        double precision, intent(in) :: box(:)
        where (abs(r) > box/2)
          r = r - sign(box, r)
        end where
    end subroutine pbc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine full_pbc(pos, box)
        double precision, intent(inout) :: pos(:,:)
        double precision, intent(in) :: box(:)
        integer :: i
        do i = 1, size(pos,2)
            where (abs(pos(:,i)) > box/2)
                pos(:,i) = pos(:,i) - sign(box, pos(:,i))
            end where
        end do
    end subroutine full_pbc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine compute_interaction(pos, box, u)
        double precision, intent(in) :: pos(:,:), box(:)
        double precision, intent(out) :: u
        double precision :: rij(size(pos,1)), rsq, uij
        logical :: uij_is_zero
        integer :: i,j
        u = 0
        do i = 1,size(pos, 2)
        do j = i+1,size(pos, 2)
            rij = pos(:, i) - pos(:, j)
            call pbc(rij,box)
            rsq = sum(rij**2)
            call is_zero(rsq, uij_is_zero)
            if (uij_is_zero) cycle
            call potential(rsq, uij)
            u = u + uij
        end do
        end do
    end subroutine compute_interaction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine compute_single_interaction(pos, box, index, u)
        double precision, intent(in) :: pos(:,:), box(:)
        integer, intent(in) :: index
        integer :: i
        double precision, intent(out) :: u
        double precision :: rij(size(pos,1)), rsq, uij
        logical :: uij_is_zero
        u = 0
        do i = 1,size(pos, 2)
            if (index == i) cycle
            rij = pos(:, i) - pos(:, index)
            call pbc(rij,box)
            rsq = sum(rij**2)
            call is_zero(rsq, uij_is_zero)
            if (uij_is_zero) cycle
            call potential(rsq, uij)
            u = u + uij
        end do
    end subroutine compute_single_interaction

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine compute_virial(pos, box, w)
        double precision, intent(in) :: pos(:,:), box(:)
        double precision, intent(out) :: w
        double precision :: rij(size(pos,1)), rsq, wij
        logical :: uij_is_zero
        integer :: i,j
        w = 0
        do i = 1,size(pos, 2)
        do j = i+1,size(pos, 2)
            rij = pos(:, i) - pos(:, j)
            call pbc(rij,box)
            rsq = sum(rij**2)
            call is_zero(rsq, uij_is_zero)
            if (uij_is_zero) cycle
            call virial(rsq, wij)
            w = w - wij
        end do
        end do
        w = w / size(box)
    end subroutine compute_virial

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine check_overlap(pos,box,n_coll,result)
        double precision, intent(in) :: pos(:,:), box(:)
        integer, intent(in) :: n_coll
        logical, intent(out) :: result
        integer :: i,j
        double precision :: rij(size(pos,1)), rsq

        outer: do i = 1, n_coll
        do j = i+1, size(pos, 2)
            rij = pos(:, i) - pos(:, j)
            call pbc(rij,box)
            rsq = sum(rij**2)
            call overlap(rsq,i,j,result)
            if (result) exit outer
        end do
        end do outer
    end subroutine check_overlap

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine check_single_overlap(pos,box,n_coll,index,result)
        double precision, intent(in) :: pos(:,:), box(:)
        integer, intent(in) :: n_coll
        integer, intent(in) :: index
        logical, intent(out) :: result
        integer :: j, last_cycle
        double precision :: rij(size(pos,1)), rsq

        if (index>n_coll) then
            last_cycle = n_coll        ! non interacting particle
        else
            last_cycle = size(pos, 2)  ! colloid
        endif

        do j = 1,last_cycle
            if (j == index) cycle
            rij = pos(:, index) - pos(:, j)
            call pbc(rij,box)
            rsq = sum(rij**2)
            call overlap(rsq,index,j,result)
            if (result) exit
        end do
    end subroutine check_single_overlap

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    subroutine check_overlap_all(pos,box,result)
        double precision, intent(in) :: pos(:,:), box(:)
        logical, intent(out) :: result
        integer :: i,j
        double precision :: rij(size(pos,1)), rsq

        outer: do i = 1,size(pos, 2)
        do j = i+1,size(pos, 2)
            rij = pos(:, i) - pos(:, j)
            call pbc(rij,box)
            rsq = sum(rij**2)
            call overlap(rsq,i,j,result)
            if (result) exit outer
        end do
        end do outer
    end subroutine check_overlap_all

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! subroutine check_single_overlap(pos,box,index,result)
    !     double precision, intent(in) :: pos(:,:), box(:)
    !     integer, intent(in) :: index
    !     logical, intent(out) :: result
    !     integer :: j
    !     double precision :: rij(size(pos,1)), rsq

    !     do j = 1,size(pos, 2)
    !         if (j == index) cycle
    !         rij = pos(:, index) - pos(:, j)
    !         call pbc(rij,box)
    !         rsq = sum(rij**2)
    !         call overlap(rsq,index,j,result)
    !         if (result) exit
    !     end do
    ! end subroutine check_single_overlap


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! subroutine interaction(pos, box, u, w)
    !     double precision, intent(in) :: pos(:,:), box(:)
    !     double precision, intent(out) :: u, w
    !     integer :: i,j
    !     double precision :: rijsq, rij(size(pos,1)), uij, wij
    !     u = 0
    !     w = 0
    !     do i = 1,size(pos, 2)
    !       do j = i+1,size(pos, 2)
    !         rij = pos(:, i) - pos(:, j)
    !         call pbc(rij, box)
    !         rijsq = sum(rij**2)
    !         call potential(rijsq, uij, wij)
    !         u = u + uij
    !         w = w - wij
    !       end do
    !     end do
    !     w = w / size(box)
    !   end subroutine interaction

end module interaction_module