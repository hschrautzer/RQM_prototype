module spk_geodesic
    !> This module contains utilities to calculate geodesic properties of the S² manifolds representing the true
    !> configuration space of the spin-hamiltonians.
    use spk_precision_const, only: xp, C_EPSILON
    use spk_input_lattice, only: N_atom
    use spk_vector_lina, only: norm, cross, angle_between3D
    use spk_math_const, only: pi
    use spk_mtprng, only: mtprng_rand_real1, mtprng_state
    use spk_rotations, only: calc_rotation_axis3D
    implicit none
    private
    public :: calculate_geodesic_distance, calculate_geodesic_path_single_spins, &
            & calculate_geodesic_pathlength_cumulated, geodesic_distance_matrix
    contains
        subroutine geodesic_distance_matrix(spin_configs,n,distance_matrix)
            !> Calculates the adjacency matrix for spin configurations
            !===================================Input Variable==========================================================
            real(kind=xp), intent(in) :: spin_configs(:,:)                      !< Spin Configurations
            integer, intent(in) :: n                                            !< Number of spin configs.
            !===================================Output Variable=========================================================
            real(kind=xp), intent(out) :: distance_matrix(:,:)                  !< Distance matrix (symmetric)
            !===================================Local Variable==========================================================
            integer :: i,j                                                      !< Indices for spin configs.
            real(kind=xp) :: dist                                               !< distance between pair
            distance_matrix = 0.0d0
            do i=1,n
                do j=i,n
                    call calculate_geodesic_distance(spin_configs(:,i),spin_configs(:,j),dist)
                    distance_matrix(i,j) = dist
                    distance_matrix(j,i) = dist
                end do
            end do

        end subroutine geodesic_distance_matrix

        subroutine calculate_geodesic_pathlength_cumulated(n_images, path, pathlen)
            !> Calculates the length of some given path, by calculating the geodesic distances between adjacent images
            !> and cumulutating
            !================================Input Variable=============================================================
            integer, intent(in) :: n_images                     !< Number of images of the path
            real(kind=xp), intent(in) :: path(3*N_atom,n_images)!< the path itself
            !================================Output Variable============================================================
            real(kind=xp), intent(out) :: pathlen(n_images)     !< cumulated path length
            !================================Local Variable=============================================================
            integer :: i_image                                  !< iterator for image
            real(kind=xp) :: distance                           !< distance between two images
            real(kind=xp) :: local_spin_k(3*N_atom)             !< dummy quantities for subroutine
            real(kind=xp) :: local_spin_p(3*N_atom)             !< dummy quantities for subroutine

            pathlen = 0.0d0

            do i_image=2,n_images
                local_spin_k = path(:,i_image)
                local_spin_p = path(:,i_image-1)
                call calculate_geodesic_distance(local_spin_p, local_spin_k, distance)
                pathlen(i_image) = pathlen(i_image-1) + distance
            end do
        end subroutine calculate_geodesic_pathlength_cumulated

        subroutine calculate_geodesic_distance(spinimage_1,spinimage_2,distance)
            !> Calculates the geodesic distance between two representations of the same spin lattice.
            !===========================================Input Variable==================================================
            real(kind=xp), intent(in) :: spinimage_1(3*N_atom)          !< spin lattice 1
            real(kind=xp), intent(in) :: spinimage_2(3*N_atom)          !< spin lattice 2
            !===========================================Output Variable=================================================
            real(kind=xp), intent(out) :: distance                      !< geodesic distance between the two images
            !===========================================Local Variable==================================================
            integer :: i                                                !< Iterator for i-th atom
            real(kind=xp) :: spin1(3), spin2(3)                         !< i-th spin from lattice 1 and from lattice 2
            real(kind=xp) :: norm_cp                                    !< norm of the cross product of spin1 and spin2
            real(kind=xp) :: dotprod                                    !< dot product of spin 1 and spin2
            real(kind=xp) :: dist_i                                     !< geodesic distance of spin1 and spin2
            real(kind=xp) :: dist_total                                 !< total geodesic distance
            dist_total = 0.0d0

            do i=1,N_atom
                spin1 = spinimage_1(3*i-2:3*i)
                spin2 = spinimage_2(3*i-2:3*i)
                norm_cp = norm(cross(spin1, spin2))
                dotprod = spin1(1) * spin2(1) + spin1(2) * spin2(2) + spin1(3) * spin2(3)
                dist_i = atan2(norm_cp, dotprod)
                dist_total = dist_total + dist_i ** 2
            end do
            distance = sqrt(dist_total)
        end subroutine calculate_geodesic_distance

        subroutine calculate_geodesic_path_single_spins(spin_i, spin_f, n_images, axis, path, state)
            !> Calculates the geodesic path for a single spin interpolated a n images.
            !================================Input Variable=============================================================
            real(kind=xp), intent(in) :: spin_i(3), spin_f(3)               !< Initial and Final Spin
            integer, intent(in) :: n_images                                 !< number of interpolation steps
            real(kind=xp), intent(inout) :: axis(3)                         !< axis of rotation for the spin
            type(mtprng_state), intent(inout) :: state
            !================================Output Variable============================================================
            real(kind=xp), intent(out) :: path(3,n_images)                  !< Interpolation of the initial towards the
                                                                            !< final spin
            !================================Local Variable=============================================================
            real(kind=xp) :: total_angle                                    !< total angle between initial and final spin
            integer :: i,j                                                  !< iterator
            real(kind=xp) :: dummy_vec(3)                                   !< dummy vector
            real(kind=xp) :: dp_spin_i_axis                                 !< Dot product between the initial spin and
                                                                            !< the axis.
            real(kind=xp) :: norm_axis                                      !< norm of the axis
            real(kind=xp) :: dtheta                                         !< Theta increment angle. Depends on the
                                                                            !< number of images.
            real(kind=xp) :: theta                                          !< current value of theta

            total_angle = angle_between3D(spin_i, spin_f)
            if (total_angle<C_EPSILON) then
                !In case the angle between the two spins is zero apply a linear interpolation
                do i=1,3
                    dummy_vec(i) = (spin_f(i) - spin_i(i)) / (n_images-1)
                end do
                path(:,1) = spin_i
                path(:,n_images) = spin_f
                do i=2,n_images-1
                    do j=1,3
                        path(j,i) = spin_i(j) + real((i-1),kind=xp) * dummy_vec(j)
                    end do
                    path(:,i) = path(:,i) / norm(path(:,i))
                end do
            elseif (abs(total_angle-pi)<C_EPSILON) then !Rotate around random axis perpendicular to the initial spin_i
                !In case the angle between the two spins is 180 degree the axis of rotation is undefined.
                dp_spin_i_axis = axis(1)*spin_i(1) + axis(2)*spin_i(2) + axis(3)*spin_i(3)
                !Get the orthogonal part of axis relative to spin_i
                axis(:) = axis(:) - dp_spin_i_axis*spin_i(:)
                norm_axis = norm(axis)
                do while (norm_axis<C_EPSILON)
                    do j=1,3
                        !Create each component of the axis randomly by the following procedure. The function sign(A,B)
                        !returns the value of A with the sign of B. B is a random variable from the interval [-1,1].
                        !Thus, the first part represent a randomly shuffled sign for the second part. The second part
                        !represents a random variable from the intervall [1,2]. Thus, we have random axis component in
                        ![-2,-1] u [1,2].
                        axis(j) = sign(1d0,2d0*mtprng_rand_real1(state)-1d0)*(mtprng_rand_real1(state)+1d0)
                    end do
                    dp_spin_i_axis = axis(1)*spin_i(1) + axis(2)*spin_i(2) + axis(3)*spin_i(3)
                    !Get the orthogonal part of axis relative to spin_i
                    axis(:) = axis(:) - dp_spin_i_axis*spin_i(:)
                    norm_axis = norm(axis)
                end do
                !Normalize the axis
                axis = axis / norm_axis
                !Now axis is a random vector orientated perpendicular to spin_i
                dtheta = pi/(n_images-1)
                path(:,1) = spin_i(:)
                path(:,n_images) = spin_f(:)
                do i=2,n_images-1
                    theta = (i-1)*dtheta
                    !Application of Rodriguez Formula: s_rot = s * cos(theta) + (a x s) * sin(theta) for the case that
                    !the rotation axis a and the vector to rotate (s) are perpendicular.
                    path(1,i) = spin_i(1)*cos(theta) + sin(theta)*(axis(2)*spin_i(3)-axis(3)*spin_i(2))
                    path(2,i) = spin_i(2)*cos(theta) - sin(theta)*(axis(1)*spin_i(3)-axis(3)*spin_i(1))
                    path(3,i) = spin_i(3)*cos(theta) + sin(theta)*(axis(1)*spin_i(2)-axis(2)*spin_i(1))
                    path(:,i) = path(:,i) / norm(path(:,i))
                end do
            else !Calculate axis based on the two configuration spin_i and spin_f
                axis = calc_rotation_axis3D(spin_i, spin_f)
                dtheta = total_angle/(n_images-1)
                path(:,1) = spin_i(:)
                path(:,n_images) = spin_f(:)
                do i=2,n_images-1
                    theta = (i-1)*dtheta
                    path(1,i) = spin_i(1)*cos(theta) + sin(theta)*(axis(2)*spin_i(3)-axis(3)*spin_i(2))
                    path(2,i) = spin_i(2)*cos(theta) - sin(theta)*(axis(1)*spin_i(3)-axis(3)*spin_i(1))
                    path(3,i) = spin_i(3)*cos(theta) + sin(theta)*(axis(1)*spin_i(2)-axis(2)*spin_i(1))
                    path(:,i) = path(:,i) / norm(path(:,i))
                end do
            end if
        end subroutine calculate_geodesic_path_single_spins


end module spk_geodesic
