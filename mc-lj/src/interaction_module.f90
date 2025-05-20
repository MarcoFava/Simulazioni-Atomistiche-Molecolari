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