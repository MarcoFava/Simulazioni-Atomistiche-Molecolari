module potential_module

IMPLICIT NONE

double precision :: r_col_par_sq, r_col_col_sq
logical, allocatable :: chem_spec_(:)

contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine initialize(r_part, r_coll, n_part, n_coll)
    double precision, intent(in) :: r_part, r_coll
    integer, intent(in) :: n_part, n_coll
    integer :: i

    r_col_col_sq = (2*r_coll)**2
    r_col_par_sq = (r_coll+r_part)**2

    ! chem_spec_ = (/ .False.,.False.,.True.,.True.,.True.,.True.,.True.,.True. ....../)
    chem_spec_ = (/ (.False., i=1,n_coll),(.True., i=1,n_part) /)

end subroutine initialize

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine overlap(rsq, i, j, result)
    double precision, intent(in) :: rsq
    integer, intent(in) :: i, j
    logical, intent(out) :: result

    if (chem_spec_(i) .and. chem_spec_(j)) then     ! interaction "particle-particle" is null
        result = .false. 
    elseif (chem_spec_(i) .or. chem_spec_(j)) then  ! interaction "colloid-particle"
        result = rsq < r_col_par_sq
    else                                            ! interaction "colloid-colloid"
        result = rsq < r_col_col_sq
    end if
end subroutine overlap



!! NON NEEDED SUBROUTINES:

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine is_zero(rsq, result)
    double precision, intent(in) :: rsq
    logical, intent(out) :: result
    return
end subroutine is_zero

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine potential(rsq, u)
    double precision, intent(in) :: rsq
    double precision, intent(out) :: u
    return
end subroutine potential

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine virial(rsq, w)
    double precision, intent(in) :: rsq
    double precision, intent(out) :: w
    return
end subroutine virial

end module potential_module