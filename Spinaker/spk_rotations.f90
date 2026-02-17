module spk_rotations
    !> Module contains rotations of 3*N_atom quantities with Rodriguez formula
    use spk_input_lattice, only: N_atom
    use spk_vector_lina, only: norm, cross
    use spk_precision_const, only: xp, C_RODRIGUEZ_TAYLOR_THRESHOLD, C_EPSILON

    implicit none
    private

    interface rotate_vector_between_tangent_spaces
        ! Interface differs between the representations of the input spin configuration
        module procedure rotate_vector_between_tangent_spaces_1darray
        module procedure rotate_vector_between_tangent_spaces_subsys
    end interface rotate_vector_between_tangent_spaces

    interface rotate_spinconfiguration
        module procedure rotate_spinconfiguration_1darray_scaler
        module procedure rotate_spinconfiguration_1darray
        module procedure rotate_spinconfiguration_subsys_disp
        module procedure rotate_spinconfiguration_subsys_disp_scaler
    end interface rotate_spinconfiguration

    public :: rotate_vector_between_tangent_spaces, rotate_spinconfiguration, calc_rotation_axis3D
    contains
        subroutine rotate_vector_between_tangent_spaces_1darray(spin, displacement, vector, vector_rotated)
            !> This routine performs the rotation of a vector from the tangent space of the current spin configuration
            !> Towards the tangent space of the configuration which is build by applying the displacement towards this
            !> current configuration.
            !> A simplification of general rodriguez formula was used which makes use of the fact that vector is in the
            !> tangent space of spin. But (different to rotation of spins) this vector is not perpendicular to the rotation
            !> axis defined by displacement and spin
            !================================Input Variable=============================================================
            real(kind=xp), intent(in) :: spin(3*N_atom)                 !< Spin configuration
            real(kind=xp), intent(in) :: displacement(3*N_atom)         !< displacement for all spins in 3N space
                                                                        !< each displacement(3*i-2:3*i) defines the disp.
                                                                        !< for a spin i towards the spin i in the next
                                                                        !< iteration. The rotation axis is defined by
                                                                        !< omega_i = (spin_i x disp_i) / |disp_i|
                                                                        !< the angle of rotation is given by theta_i = |disp_i|
            real(kind=xp), intent(in) :: vector(3*N_atom)               !< vector in the tangent space of the current
                                                                        !< spin configuration
            !================================Output Variable============================================================
            real(kind=xp), intent(out) :: vector_rotated(3*N_atom)      !< vector rotated into the tangent space frame of
                                                                        !< the displaced configuration
            !================================Local Variable=============================================================
            real(kind=xp) :: angle                                      !< angle of rotation of current spin
            real(kind=xp) :: rotation_axis_nn(3)                        !< the rotation axis but not normed
            integer :: i                                                !< iterate over the atoms
            vector_rotated = vector
            if (.not.(all(abs(vector)<=C_EPSILON))) then
                do i=1,N_atom
                    angle = norm(displacement(3*i-2:3*i))
                    rotation_axis_nn = cross(spin(3*i-2:3*i), displacement(3*i-2:3*i))
                    if (angle>=C_RODRIGUEZ_TAYLOR_THRESHOLD) then
                        vector_rotated(3*i-2:3*i) = vector(3*i-2:3*i) * cos(angle) &
                                & - spin(3*i-2:3*i) * (displacement(3*i-2) * vector(3*i-2) &
                                        & +displacement(3*i-1) * vector(3*i-1) &
                                        & +displacement(3*i) * vector(3*i)) * (sin(angle) / angle) &
                                & + rotation_axis_nn*(rotation_axis_nn(1) * vector(3*i-2) &
                                        &+rotation_axis_nn(2) * vector(3*i-1) &
                                        &+rotation_axis_nn(3) * vector(3*i)) * ((1 - cos(angle)) / (angle * angle))
                    else
                        vector_rotated(3*i-2:3*i) = vector(3*i-2:3*i)*cos(angle) &
                                & - spin(3*i-2:3*i) * (displacement(3*i-2) * vector(3*i-2) &
                                        & +displacement(3*i-1) * vector(3*i-1) &
                                        & +displacement(3*i) * vector(3*i)) * (1.0d0 - (1.0d0 / 6.0d0) &
                                        & * (angle * angle) + (1.0d0 / 120.0d0) &
                                        & * (angle * angle * angle * angle)) &
                                & + rotation_axis_nn*(rotation_axis_nn(1) * vector(3*i-2) &
                                        &+rotation_axis_nn(2) * vector(3*i-1) &
                                        &+rotation_axis_nn(3) * vector(3*i)) * (0.5d0 - (1.0d0 / 24.0d0) &
                                        & * (angle * angle) + (1.0d0 / 720.0d0) &
                                        & * (angle * angle * angle * angle))
                        if (any(abs(vector(3*i-2:3*i))>=C_EPSILON)) then
                            vector_rotated(3*i-2:3*i) = vector_rotated(3*i-2:3*i) * norm(vector(3*i-2:3*i)) / &
                                    & norm(vector_rotated(3*i-2:3*i))
                        else
                            vector_rotated(3*i-2:3*i)  = vector(3*i-2:3*i)
                        end if
                    end if
                end do
            end if


        end subroutine rotate_vector_between_tangent_spaces_1darray

        subroutine rotate_vector_between_tangent_spaces_subsys(spin, subsys_idx, displacement, vector, vector_rotated)
            !> This routine performs the rotation of a vector from the tangent space of the current spin configuration
            !> Towards the tangent space of the configuration which is build by applying the displacement towards this
            !> current configuration.
            !> A simplification of general rodriguez formula was used which makes use of the fact that vector is in the
            !> tangent space of spin. But (different to rotation of spins) this vector is not perpendicular to the rotation
            !> axis defined by displacement and spin. The displacement as well as the vectors are only defined in a sub-
            !> system of the spin configuration
            !================================Input Variable=============================================================
            integer, intent(in) :: subsys_idx(:)                    !< Index Array defining the subsystem
            real(kind=xp), intent(in) :: spin(:)                    !< Spin configuration (full 3N system
            real(kind=xp), intent(in) :: displacement(:)            !< Displacement for all spins in subsystem
                                                                    !< Each displacement(3*i-2:3*i) defines the disp.
                                                                    !< for a spin i towards the spin i in the next
                                                                    !< iteration. The rotation axis is defined by
                                                                    !< omega_i = (spin_i x disp_i) / |disp_i|
                                                                    !< the angle of rotation is given by theta_i = |disp_i|
            real(kind=xp), intent(in) :: vector(:)                  !< vector in the tangent space of the current
                                                                    !< spin configuration, defined for subsystem
            !================================Output Variable============================================================
            real(kind=xp), intent(out) :: vector_rotated(:)         !< vector rotated into the tangent space frame of
                                                                    !< the displaced configuration
            !================================Local Variable=============================================================
            integer :: N_sub                                        !< Number of atoms in the subsystem
            real(kind=xp) :: angle                                  !< angle of rotation of current spin
            real(kind=xp) :: rotation_axis_nn(3)                    !< the rotation axis but not normed
            integer :: i                                            !< iterate over the atoms
            N_sub = size(subsys_idx)
            vector_rotated = vector
            if (.not.(all(abs(vector)<=C_EPSILON))) then
                do i=1,N_sub
                    angle = norm(displacement(3*i-2:3*i))
                    rotation_axis_nn = cross(spin(3*subsys_idx(i)-2:3*subsys_idx(i)), displacement(3*i-2:3*i))
                    if (angle>=C_RODRIGUEZ_TAYLOR_THRESHOLD) then
                        vector_rotated(3*i-2:3*i) = vector(3*i-2:3*i) * cos(angle) &
                                & - spin(3*subsys_idx(i)-2:3*subsys_idx(i)) * (displacement(3*i-2) * vector(3*i-2) &
                                        & +displacement(3*i-1) * vector(3*i-1) &
                                        & +displacement(3*i) * vector(3*i)) * (sin(angle) / angle) &
                                & + rotation_axis_nn*(rotation_axis_nn(1) * vector(3*i-2) &
                                        &+rotation_axis_nn(2) * vector(3*i-1) &
                                        &+rotation_axis_nn(3) * vector(3*i)) * ((1 - cos(angle)) / (angle * angle))
                    else
                        vector_rotated(3*i-2:3*i) = vector(3*i-2:3*i)*cos(angle) &
                                & - spin(3*subsys_idx(i)-2:3*subsys_idx(i)) * (displacement(3*i-2) * vector(3*i-2) &
                                        & +displacement(3*i-1) * vector(3*i-1) &
                                        & +displacement(3*i) * vector(3*i)) * (1.0d0 - (1.0d0 / 6.0d0) &
                                        & * (angle * angle) + (1.0d0 / 120.0d0) &
                                        & * (angle * angle * angle * angle)) &
                                & + rotation_axis_nn*(rotation_axis_nn(1) * vector(3*i-2) &
                                        &+rotation_axis_nn(2) * vector(3*i-1) &
                                        &+rotation_axis_nn(3) * vector(3*i)) * (0.5d0 - (1.0d0 / 24.0d0) &
                                        & * (angle * angle) + (1.0d0 / 720.0d0) &
                                        & * (angle * angle * angle * angle))
                        if (any(abs(vector(3*i-2:3*i))>=C_EPSILON)) then
                            vector_rotated(3*i-2:3*i) = vector_rotated(3*i-2:3*i) * norm(vector(3*i-2:3*i)) / &
                                    & norm(vector_rotated(3*i-2:3*i))
                        else
                            vector_rotated(3*i-2:3*i)  = vector(3*i-2:3*i)
                        end if
                    end if
                end do
            end if


        end subroutine rotate_vector_between_tangent_spaces_subsys

        subroutine rotate_spinconfiguration_1darray_scaler(spin, displacement, scaler)
            !> Performs the rotation of the spin configuration along displacement following Rodriguez formula. In the
            !> limit of very small displacements the taylor expansion of the rodriguez formula is used. The angle of rotation
            !> of all spins can be scaled.
            !================================Input Variable=============================================================
            real(kind=xp), intent(in) :: displacement(3*N_atom)     !< displacement vector in tangent space of input spin
            real(kind=xp), intent(in) :: scaler                     !< scaling of rotation for all spins
            !================================Inout Variable=============================================================
            real(kind=xp), intent(inout) :: spin(3*N_atom)          !< input spin configuration, changed on output
            !================================Local Variable=============================================================
            integer :: i                                            !< iterates the atoms
            real(kind=xp) :: angle                                  !< rotation angle for spin i
            real(kind=xp) :: new_spin(3)                            !< new spin i
            real(kind=xp) :: displacement_dummy(3*N_atom)           !< local copy of displacement will be scaled

            displacement_dummy = scaler * displacement
            do i=1,N_atom
                angle = norm(displacement_dummy(3*i-2:3*i))
                if (angle>=C_RODRIGUEZ_TAYLOR_THRESHOLD) then
                    new_spin=spin(3*i-2:3*i) * cos(angle) + displacement_dummy(3*i-2:3*i) * sin(angle) / angle
                else
                    new_spin=spin(3*i-2:3*i) * cos(angle) + displacement_dummy(3*i-2:3*i) &
                            & * (1.0d0 - (1.0d0 / 6.0d0) * angle**2 + (1.0d0 / 120.0d0) * angle**4)
                end if
                spin(3*i-2:3*i) = new_spin / norm(new_spin)
            end do
        end subroutine rotate_spinconfiguration_1darray_scaler

        subroutine rotate_spinconfiguration_1darray(spin, displacement)
            !> Performs the rotation of the spin configuration along displacement following Rodriguez formula. In the
            !> limit of very small displacements the taylor expansion of the rodriguez formula is used.
            !================================Input Variable=============================================================
            real(kind=xp), intent(in) :: displacement(3*N_atom)     !< displacement vector in tangent space of input spin
            !================================Inout Variable=============================================================
            real(kind=xp), intent(inout) :: spin(3*N_atom)          !< input spin configuration, changed on output
            !================================Local Variable=============================================================
            integer :: i                                            !< iterates the atoms
            real(kind=xp) :: angle                                  !< rotation angle for spin i
            real(kind=xp) :: new_spin(3)                            !< new spin i

            do i=1,N_atom
                angle = norm(displacement(3*i-2:3*i))
                if (angle>=C_RODRIGUEZ_TAYLOR_THRESHOLD) then
                    new_spin=spin(3*i-2:3*i) * cos(angle) + displacement(3*i-2:3*i) * sin(angle) / angle
                else
                    new_spin=spin(3*i-2:3*i) * cos(angle) + displacement(3*i-2:3*i) &
                            & * (1.0d0 - (1.0d0 / 6.0d0) * angle**2 + (1.0d0 / 120.0d0) * angle**4)
                end if
                spin(3*i-2:3*i) = new_spin / norm(new_spin)
            end do
        end subroutine rotate_spinconfiguration_1darray

        subroutine rotate_spinconfiguration_subsys_disp(spin,subsys_idx,displacement)
            !> Performs the rotation of the spin configuration along displacement following Rodriguez formula. In the
            !> limit of very small displacements the taylor expansion of the rodriguez formula is used. The displacement
            !> is only defined in a subsystem of the spin latttice. Thus only the spins belonging to that subsystem will
            !> be rotated
            !================================Input Variable=============================================================
            integer, intent(in) :: subsys_idx(:)             !< Indices of spins belonging to the subsystem
            real(kind=xp), intent(in) :: displacement(:)     !< Displacement vector in tangent space of input spin in a
                                                             !< Subsystem of spin
            !================================Inout Variable=============================================================
            real(kind=xp), intent(inout) :: spin(:)          !< Input spin configuration, changed on output
            !================================Local Variable=============================================================
            integer :: N_sub                                 !< Number of spins in subsystem
            integer :: i                                     !< Iterates the atoms
            real(kind=xp) :: angle                           !< Rotation angle for spin i
            real(kind=xp) :: new_spin(3)                     !< New spin i

            N_sub = size(subsys_idx)
            do i=1,N_sub
                angle = norm(displacement(3*i-2:3*i))
                if (angle>=C_RODRIGUEZ_TAYLOR_THRESHOLD) then
                    new_spin=spin(3*subsys_idx(i)-2:3*subsys_idx(i)) * cos(angle) + displacement(3*i-2:3*i) * sin(angle) / angle
                else
                    new_spin=spin(3*subsys_idx(i)-2:3*subsys_idx(i)) * cos(angle) + displacement(3*i-2:3*i) &
                            & * (1.0d0 - (1.0d0 / 6.0d0) * angle**2 + (1.0d0 / 120.0d0) * angle**4)
                end if
                spin(3*subsys_idx(i)-2:3*subsys_idx(i)) = new_spin / norm(new_spin)
            end do
        end subroutine rotate_spinconfiguration_subsys_disp

        subroutine rotate_spinconfiguration_subsys_disp_scaler(spin,subsys_idx,displacement,scaler)
            !> Performs the rotation of the spin configuration along displacement following Rodriguez formula. In the
            !> limit of very small displacements the taylor expansion of the rodriguez formula is used. The displacement
            !> is only defined in a subsystem of the spin latttice. Thus only the spins belonging to that subsystem will
            !> be rotated
            !================================Input Variable=============================================================
            integer, intent(in) :: subsys_idx(:)             !< Indices of spins belonging to the subsystem
            real(kind=xp), intent(in) :: displacement(:)     !< Displacement vector in tangent space of input spin in a
                                                             !< Subsystem of spin
            real(kind=xp), intent(in) :: scaler              !< Scaler for the displacement
            !================================Inout Variable=============================================================
            real(kind=xp), intent(inout) :: spin(:)          !< Input spin configuration, changed on output
            !================================Local Variable=============================================================
            real(kind=xp), allocatable :: displacement_dummy(:)
            integer :: N_sub                                 !< Number of spins in subsystem
            integer :: i                                     !< Iterates the atoms
            real(kind=xp) :: angle                           !< Rotation angle for spin i
            real(kind=xp) :: new_spin(3)                     !< New spin i

            N_sub = size(subsys_idx)
            allocate(displacement_dummy(size(displacement)))
            displacement_dummy = displacement * scaler
            do i=1,N_sub
                angle = norm(displacement_dummy(3*i-2:3*i))
                if (angle>=C_RODRIGUEZ_TAYLOR_THRESHOLD) then
                    new_spin=spin(3*subsys_idx(i)-2:3*subsys_idx(i)) * cos(angle) + displacement_dummy(3*i-2:3*i) &
                            & * sin(angle) / angle
                else
                    new_spin=spin(3*subsys_idx(i)-2:3*subsys_idx(i)) * cos(angle) + displacement_dummy(3*i-2:3*i) &
                            & * (1.0d0 - (1.0d0 / 6.0d0) * angle**2 + (1.0d0 / 120.0d0) * angle**4)
                end if
                spin(3*subsys_idx(i)-2:3*subsys_idx(i)) = new_spin / norm(new_spin)
            end do
        end subroutine rotate_spinconfiguration_subsys_disp_scaler

        function calc_rotation_axis3D(vec_i,vec_f)
            !> Calculates a rotation axis for rotation from vec_i towards the orientation of vec_f
            !================================================Input Variable=============================================
            real(kind=xp), intent(in) :: vec_i(3),vec_f(3)          !< Initial and final vector of rotation
            !================================================Output Variable============================================
            real(kind=xp) :: calc_rotation_axis3D(3)   !< Calculated axis of rotation
            !================================================Local Variable=============================================
            real(kind=xp) :: norm_final                             !< Norm of final
            real(kind=xp) :: stretch=1.0d10                         !< stretch the final vector if it is too small
            real(kind=xp) :: vec_f_local(3)                         !< local copy of vec_f
            real(kind=xp) :: y(3)
            real(kind=xp) :: x(3)
            norm_final = norm(vec_f)
            if (norm_final<1.0d-12) then
                vec_f_local = vec_f * stretch
            else
                vec_f_local = vec_f
            end if
            ! gives back the absolut value of a times the sign of f_in(1)
            y(1) = sign(vec_f_local(1),vec_f(1))
            y(2) = sign(vec_f_local(2),vec_f(2))
            y(3) = sign(vec_f_local(3),vec_f(3))
            x = cross(vec_i,y)
            calc_rotation_axis3D = x / norm(x)

        end function calc_rotation_axis3D

end module spk_rotations
