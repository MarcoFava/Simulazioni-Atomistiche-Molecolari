module potential_module

IMPLICIT NONE

double precision :: sigma_, sigsq_, epsilon_, rcut_, ucut_

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initialize(sigma, epsilon, rcut)
    double precision, intent(in) :: sigma, epsilon, rcut
    sigma_ = sigma
    sigsq_ = sigma**2
    epsilon_ = epsilon
    rcut_ = rcut
    ucut_ = 0.
    call potential(rcut**2,ucut_)
end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine is_zero(rsq, result)
    double precision, intent(in) :: rsq
    logical, intent(out) :: result
    result = rsq > rcut_**2
end subroutine is_zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine potential(rsq, u)
    double precision, intent(in) :: rsq
    double precision, intent(out) :: u
    u = 4 * epsilon_ * ((sigsq_/rsq)**6 - (sigsq_/rsq)**3) - ucut_
end subroutine potential

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine virial(rsq, w)
    double precision, intent(in) :: rsq
    double precision, intent(out) :: w
    ! w = (du/dr)*r
    w = - 24 * epsilon_ * (2*(sigsq_/rsq)**6 - (sigsq_/rsq)**3)
end subroutine virial

end module potential_module