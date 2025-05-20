! initialize.f90
! Sets up initial configuration for MD or MC
MODULE initialize_driver

  ! Reads several variables and options from standard input using a namelist nml
  ! Leave namelist empty to accept supplied defaults

  USE, INTRINSIC :: iso_fortran_env, ONLY : input_unit, output_unit, error_unit, iostat_end, iostat_eor, &
       &                                    COMPILER_VERSION, COMPILER_OPTIONS
  
  USE config_io_module,  ONLY : write_cnf_atoms, write_cnf_mols
  USE maths_module,      ONLY : lowercase
  USE initialize_module, ONLY : allocate_arrays, deallocate_arrays, &
       &                        fcc_positions, ran_positions, ran_velocities, &
       &                        chain_positions, chain_velocities, &
       &                        n, r, e, v, w

  IMPLICIT NONE

  ! Most important variables
  INTEGER :: nc          ! Number of fcc unit cells in each coordinate direction; n = 4*nc**3
  REAL    :: temperature ! Specified temperature (used in generating velocities)
  REAL    :: inertia     ! Specified moment of inertia (for linear and nonlinear molecules)
  REAL    :: density     ! Specified number density (ignored for chain molecules)
  LOGICAL :: velocities  ! User option requiring velocities for MD (otherwise just positions)
  LOGICAL :: lattice     ! User option for lattice positions (otherwise placed randomly)
  REAL    :: length      ! Spherocylinder length (linear molecules, overlap check)
  REAL    :: bond        ! Bond length (only used for chain molecules)
  REAL    :: box         ! Box length (ignored for chain molecules)
  LOGICAL :: soft        ! Option for soft interactions (i.e. ignore overlaps)
  LOGICAL :: constraints ! Option to apply constraints on chain velocities


  CONTAINS
  
   SUBROUTINE initialize()

      CHARACTER(len=10)  :: molecules       ! Character string, input, used to specify molecule_option
      INTEGER            :: molecule_option ! User option for atoms, linear, or nonlinear molecule, or chain
      INTEGER, PARAMETER :: atom = 0, linear = 1, nonlinear = 2, chain = 3

      INTEGER            :: ioerr

      REAL, PARAMETER :: tol = 1.0e-6

      CHARACTER(len=7), PARAMETER :: filename = 'cnf.inp' ! Will be used as an input file by later simulations

      NAMELIST /nml/ nc, n, temperature, inertia, density, length, velocities, molecules, lattice, soft, constraints

      WRITE ( unit=output_unit, fmt='(a)'   ) 'initialize'
      WRITE ( unit=output_unit, fmt='(2a)'  ) 'Compiler: ', COMPILER_VERSION()
      WRITE ( unit=output_unit, fmt='(2a/)' ) 'Options:  ', COMPILER_OPTIONS()
      
      WRITE ( unit=output_unit, fmt='(a)' ) 'Sets up initial configuration file for various simulations'
      WRITE ( unit=output_unit, fmt='(a)' ) 'Options for molecules are "atom", "linear", "nonlinear", "chain"'
      WRITE ( unit=output_unit, fmt='(a)' ) 'Particle mass m=1 throughout'

      CALL RANDOM_INIT ( .FALSE., .TRUE. ) ! Initialize random number generator

      ! Set default parameters
      n           = 0        ! nc takes precedence unless n is explicitly specified
      !   nc          = 4        ! Default is N = 4*(4**3) = 256 on a fcc lattice, a small system
      temperature = 1.0      ! Should lie in the liquid region for density > 0.7 or so
      inertia     = 1.0      ! Only relevant for linear and nonlinear molecules
      !   density     = 0.75     ! Should lie in the liquid region for temperature > 0.9 or so
      length      = 0.0      ! By default, atoms are spherical
      bond        = 1.122462 ! Bond length for chain molecules
      velocities  = .FALSE.  ! By default, produce positions only, for MC simulations
      molecules   = 'atoms'  ! Options are 'atoms', 'chain', 'linear', 'nonlinear'
      lattice     = .TRUE.   ! By default, arrange atoms on a lattice with predetermined orientations
      soft        = .FALSE.  ! By default, check for overlaps when placing molecules
      constraints = .TRUE.   ! By default, constrain chain velocities relative to bonds

      ! Read run parameters from namelist
      ! Comment out, or replace, this section if you don't like namelists
      !   READ ( unit=input_unit, nml=nml, iostat=ioerr )
      !   IF ( ioerr /= 0 ) THEN
      !      WRITE ( unit=error_unit, fmt='(a,i15)') 'Error reading namelist nml from standard input', ioerr
      !      IF ( ioerr == iostat_eor ) WRITE ( unit=error_unit, fmt='(a)') 'End of record'
      !      IF ( ioerr == iostat_end ) WRITE ( unit=error_unit, fmt='(a)') 'End of file'
      !      STOP 'Error in initialize'
      !   END IF

      ! Use molecules string to deduce molecule_option
      molecules = lowercase(molecules)

      IF ( INDEX ( molecules, 'nonlinear' ) /= 0 ) THEN

         molecule_option = nonlinear
         WRITE ( unit=output_unit, fmt='(a)' ) 'Nonlinear molecules'

      ELSE IF ( INDEX ( molecules, 'linear') /= 0 ) THEN

         molecule_option = linear
         WRITE ( unit=output_unit, fmt='(a)' ) 'Linear molecules'

      ELSE IF ( INDEX ( molecules, 'atom') /= 0 ) THEN

         molecule_option = atom
         WRITE ( unit=output_unit, fmt='(a)' ) 'Atoms'

      ELSE IF ( INDEX ( molecules, 'chain') /= 0 ) THEN

         molecule_option = chain
         WRITE ( unit=output_unit, fmt='(a)' ) 'Atoms in a chain'

      ELSE

         WRITE ( unit=error_unit, fmt='(a,a)') 'Unrecognized molecules option: ', molecules
         STOP 'Error in initialize'

      END IF

      IF ( n <= 0 ) THEN ! Test for unspecified N

         IF ( nc <= 0 ) THEN ! Test for unphysical nc

            WRITE ( unit=error_unit, fmt='(a,2i15)') 'nc must be positive', nc
            STOP 'Error in initialize'

         ELSE ! Deduce N from nc

            WRITE ( unit=output_unit, fmt='(a,t40,i15)') 'nc = ', nc
            n = 4*nc**3
            WRITE ( unit=output_unit, fmt='(a,t40,i15)') 'n = ', n

         END IF

      ELSE ! N has been specified

         WRITE ( unit=output_unit, fmt='(a,t40,i15)') 'n = ', n

      END IF ! End test for unspecified N

      IF ( velocities ) THEN
         WRITE ( unit=output_unit, fmt='(a)' ) 'Velocities option selected'

         ! Inertia should be positive, even for atoms
         IF ( inertia < tol ) THEN
            WRITE ( unit=output_unit, fmt='(a,f15.6)' ) 'Warning, inertia = ', inertia
            WRITE ( unit=output_unit, fmt='(a)' ) 'Resetting to 1 '
            inertia = 1.0
         END IF

      ELSE
         WRITE ( unit=output_unit, fmt='(a)' ) 'No velocities option selected'
      END IF

      SELECT CASE ( molecule_option )

      CASE ( nonlinear )
         WRITE ( unit=output_unit, fmt='(a)' ) 'Periodic boundary conditions'

      CASE ( linear )

         WRITE ( unit=output_unit, fmt='(a)' ) 'Periodic boundary conditions'

         IF ( length < tol ) THEN
            WRITE ( unit=output_unit, fmt='(a,f15.6)' ) 'Warning, length = ', length
         END IF

      CASE ( atom )

         WRITE ( unit=output_unit, fmt='(a)' ) 'Periodic boundary conditions'
         IF ( ABS(length) > tol ) THEN
            WRITE ( unit=output_unit, fmt='(a,f15.6)' ) 'Warning, length = ', length
            WRITE ( unit=output_unit, fmt='(a)' ) 'Resetting to zero'
            length = 0.0
         END IF

      CASE ( chain )

         WRITE ( unit=output_unit, fmt='(a)' ) 'NO periodic boundary conditions'
         IF ( ABS(length) > tol ) THEN
            WRITE ( unit=output_unit, fmt='(a,f15.6)' ) 'Warning, length = ', length
            WRITE ( unit=output_unit, fmt='(a)' ) 'Resetting to zero'
            length = 0.0
         END IF
         IF ( velocities ) THEN
            IF ( constraints ) THEN
               WRITE ( unit=output_unit, fmt='(a)' ) 'Velocities constrained relative to bonds'
            ELSE
               WRITE ( unit=output_unit, fmt='(a)' ) 'Velocities not constrained relative to bonds'
            END IF
         END IF

      END SELECT

      ! Allocate arrays appropriately according to molecule_option
      CALL allocate_arrays ( quaternions = ( molecule_option == nonlinear ) )

      ! Initialize positions and optionally velocities

      IF ( soft ) WRITE ( unit=output_unit, fmt='(a)' ) 'Soft option selected - no overlap checking'

      IF ( molecule_option == chain ) THEN

         WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Bond length', bond
         CALL chain_positions ( bond, soft ) ! Random with chosen bond length
         IF ( velocities ) CALL chain_velocities ( temperature, constraints )

      ELSE

         ! Periodic boundaries apply
         ! Box length is deduced from density
         box = ( REAL(n) / density ) ** ( 1.0/3.0 )
         WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Density',    density
         WRITE ( unit=output_unit, fmt='(a,t40,f15.6)' ) 'Box length', box

         IF ( lattice ) THEN
            CALL fcc_positions ( box, length, soft ) ! On a lattice within box
         ELSE
            CALL ran_positions ( box, length, soft ) ! Random within box
         END IF

         IF ( velocities ) CALL ran_velocities ( temperature, inertia )

      END IF

      ! Write out configuration

      WRITE ( unit=output_unit, fmt='(a,a)' ) 'Writing configuration to filename ', filename

      SELECT CASE ( molecule_option )

      CASE ( atom )

         ! Write out coordinates in same units as box
         IF ( velocities ) THEN
            CALL write_cnf_atoms ( filename, n, box, r, v )
         ELSE
            CALL write_cnf_atoms ( filename, n, box, r )
         END IF

      CASE ( linear, nonlinear )

         ! Write out coordinates in same units as box
         IF ( velocities ) THEN
            CALL write_cnf_mols ( filename, n, box, r, e, v, w )
         ELSE
            CALL write_cnf_mols ( filename, n, box, r, e )
         END IF

      CASE ( chain )

         ! We do not use periodic boundaries for this system
         ! Instead, use "box" variable to store bond length
         IF ( velocities ) THEN
            CALL write_cnf_atoms ( filename, n, bond, r, v )
         ELSE
            CALL write_cnf_atoms ( filename, n, bond, r )
         END IF

      END SELECT

      CALL deallocate_arrays

   END SUBROUTINE initialize
END MODULE initialize_driver
