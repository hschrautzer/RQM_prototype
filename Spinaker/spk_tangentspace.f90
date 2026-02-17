module spk_tangentspace
    !> Contains routines to project vectors from and to the tanget space of a spin configuration. Also routines
    !> for operating in the tangent space. Also creation of tangent-space basis.
    use spk_input_lattice, only: N_atom
    use spk_precision_const, only: xp, C_EPSILON
    use spk_vector_lina, only: norm, cross
    use spk_bilinear_neighbors, only: max_totalneigh_allbds, bilinear_neigh_allgrp_allbds
    use spk_coordinatesystems, only: cartesian_to_unitspherical
    use spk_products, only: AB
    use spk_sorting, only: is_in

    use omp_lib, only: omp_get_wtime

    implicit none
    private
    interface vec_to_embeddingspace
        module procedure vec_to_embeddingspace_Natom
        module procedure vec_to_embeddingspace_N
    end interface vec_to_embeddingspace

    interface to_tangentspace_grassmann_manifold
        module procedure to_tangentspace_grassmann_manifold_2D
        module procedure to_tangentspace_grassmann_manifold_1D
    end interface to_tangentspace_grassmann_manifold

    interface vec_to_tangentspace
        ! Interface differs between the representations of the input spin configuration
        module procedure vec_to_tangentspace_1darray
        module procedure vec_to_tangentspace_by_basis
        module procedure vec_to_tangentspace_subsys
    end interface vec_to_tangentspace

    interface create_tangentspace_basis
        module procedure create_tangentspace_basis_subsys
        module procedure create_tangentspace_basis_full
    end interface

    interface project_hessian3N_hessian2N_dense
        ! Interface for the computation of the projection from the 3N hessian towards the 2N hessian where all matrices
        ! are in dense format. However, this can be performed naive by multiplying all elements or only using the
        ! elements which are actually interacting according to the neighbor list.
        module procedure project_hessian3N_hessian2N_dense_all_elements
        module procedure project_hessian3N_hessian2N_dense_neighbors
        module procedure project_hessian3N_hessian2N_subsystem_dense_neighbors
    end interface project_hessian3N_hessian2N_dense

    public :: vec_to_tangentspace, vec_to_embeddingspace, create_tangentspace_basis, project_hessian3N_bilinear_dense, &
            & create_tangentspace_spherical_basis, project_3x3_to_2x2_block,to_tangentspace_stiefel_manifold, &
            & to_tangentspace_grassmann_manifold, project_hessian3N_hessian2N_dense
    contains
        subroutine vec_to_tangentspace_1darray(spin, vec, vec_tangentspace)
            !> Projects a 3*N_atom vector the tangent space of the spin configuration also in 3*N_atom representation by
            !> calculating vec_tangentspace = vec - dot(vec,spin) * spin
            !================================Input Variable=============================================================
            real(kind=xp), intent(in) :: spin(3*N_atom), vec(3*N_atom)
            !================================Output Variable============================================================
            real(kind=xp), intent(out) :: vec_tangentspace(3*N_atom)
            !================================Local Variable=============================================================
            integer :: i_atom
            real(kind=xp) :: dotprod
            vec_tangentspace = 0.0d0
            do i_atom=1,N_atom
                dotprod = vec(3*i_atom-2) * spin(3*i_atom-2) + vec(3*i_atom-1) * spin(3*i_atom-1) + vec(3*i_atom) * spin(3*i_atom)
                vec_tangentspace(3*i_atom-2:3*i_atom) = vec(3*i_atom-2:3*i_atom) - dotprod * spin(3*i_atom-2:3*i_atom)
            end do
        end subroutine vec_to_tangentspace_1darray

        subroutine vec_to_tangentspace_subsys(spin,subsys_idxlist ,vec, vec_tangentspace)
            !> Projects a 3*N_sub vector the tangent space of the spin configuration in 3*N_atom representation by
            !> calculating vec_tangentspace = vec - dot(vec,spin) * spin
            !================================Input Variable=============================================================
            real(kind=xp), intent(in) :: spin(:), vec(:)
            integer, intent(in) :: subsys_idxlist(:)
            !================================Output Variable============================================================
            real(kind=xp), intent(out) :: vec_tangentspace(:)
            !================================Local Variable=============================================================
            integer :: i, N_sub
            real(kind=xp) :: dotprod
            N_sub = size(subsys_idxlist)
            vec_tangentspace = 0.0d0
            do i=1,N_sub
                dotprod = vec(3*i-2) * spin(3*subsys_idxlist(i)-2) + vec(3*i-1) * spin(3*subsys_idxlist(i)-1) + &
                        & vec(3*i) * spin(3*subsys_idxlist(i))
                vec_tangentspace(3*i-2:3*i) = vec(3*i-2:3*i) - dotprod * spin(3*subsys_idxlist(i)-2:3*subsys_idxlist(i))
            end do
        end subroutine vec_to_tangentspace_subsys

        subroutine vec_to_tangentspace_by_basis(basis,vec,vec_tangentspace)
            !> Projects a 3*N vector to the tangent space of the spin configuration using a provided basis of this
            !> tangent space (can also be used for subsystems
            !================================Input Variable=============================================================
            real(kind=xp), intent(in) :: basis(:,:,:)                   !< tangent space basis
            real(kind=xp), intent(in) :: vec(:)                         !< vector in tangent space
            !================================Output Variable============================================================
            real(kind=xp), intent(out) :: vec_tangentspace(:)           !< vector in embedding 3N euclidean space
            !================================Local Variable=============================================================
            integer :: shape_basis(3)
            integer :: i
            shape_basis = shape(basis)
            vec_tangentspace = 0.0d0
            do i=1,shape_basis(3)
                vec_tangentspace(2*i-1) = basis(1,1,i)*vec(3*i-2) + basis(2,1,i)*vec(3*i-1) + basis(3,1,i)*vec(3*i)
                vec_tangentspace(2*i) = basis(1,2,i)*vec(3*i-2) + basis(2,2,i)*vec(3*i-1) + basis(3,2,i)*vec(3*i)
            end do
        end subroutine vec_to_tangentspace_by_basis

        subroutine vec_to_embeddingspace_Natom(basis, vec_tangentspace, vec)
            !> Projects a vector from the tangent space (described by the tangent space basis of the spin config.) to the
            !> embedding space
            !================================Input Variable=============================================================
            real(kind=xp), intent(in) :: basis(3,2,N_atom)                  !< tangent space basis
            real(kind=xp), intent(in) :: vec_tangentspace(2*N_atom)         !< vector in tangent space
            !================================Output Variable============================================================
            real(kind=xp), intent(out) :: vec(3*N_atom)                     !< vector in embedding 3N euclidean space
            !================================Local Variable=============================================================
            integer :: i
            vec = 0.0d0
            do i=1,N_atom
                vec(3*i-2)=basis(1,1,i)*vec_tangentspace(2*i-1)+basis(1,2,i)*vec_tangentspace(2*i)
                vec(3*i-1)=basis(2,1,i)*vec_tangentspace(2*i-1)+basis(2,2,i)*vec_tangentspace(2*i)
                vec(3*i)=basis(3,1,i)*vec_tangentspace(2*i-1)+basis(3,2,i)*vec_tangentspace(2*i)
            end do
        end subroutine vec_to_embeddingspace_Natom

        subroutine vec_to_embeddingspace_N(basis,N,vec_tangentspace, vec)
            !> Projects a vector from the tangent space (described by the tangent space basis of the spin config.) to the
            !> embedding space
            !================================Input Variable=============================================================
            integer, intent(in) :: N                                   !< Number of spins involved
            real(kind=xp), intent(in) :: basis(3,2,N)                  !< tangent space basis
            real(kind=xp), intent(in) :: vec_tangentspace(2*N)         !< vector in tangent space
            !================================Output Variable============================================================
            real(kind=xp), intent(out) :: vec(3*N)                     !< vector in embedding 3N euclidean space
            !================================Local Variable=============================================================
            integer :: i
            vec = 0.0d0
            do i=1,N
                vec(3*i-2)=basis(1,1,i)*vec_tangentspace(2*i-1)+basis(1,2,i)*vec_tangentspace(2*i)
                vec(3*i-1)=basis(2,1,i)*vec_tangentspace(2*i-1)+basis(2,2,i)*vec_tangentspace(2*i)
                vec(3*i)=basis(3,1,i)*vec_tangentspace(2*i-1)+basis(3,2,i)*vec_tangentspace(2*i)
            end do
        end subroutine vec_to_embeddingspace_N



        subroutine create_tangentspace_spherical_basis(spin, basis)
            !> Creates a basis of the tangent space with basis vectors equal to the spherical unit vector e_phi and e_theta
            !================================Input Variable=============================================================
            real(kind=xp), intent(in) :: spin(3*N_atom)                         !< Spin configuration
            !================================Output Variable============================================================
            real(kind=xp), intent(out) :: basis(3,2,N_atom)                     !< basis for each spin (3x2 block)
            !================================Local Variable=============================================================
            real(kind=xp) :: theta, phi
            real(kind=xp) :: e_phi(3), e_theta(3)
            real(kind=xp) :: ephiz
            integer :: i, j

            ephiz = 0.0d0
            do i=1,N_atom
                call cartesian_to_unitspherical(spin(3*i-2),spin(3*i-1),spin(3*i),phi,theta)
                e_phi = (/-sin(phi), cos(phi), ephiz /)
                e_theta = (/cos(theta)*cos(phi),cos(theta)*sin(phi),-sin(theta)/)
                e_phi = e_phi / norm(e_phi)
                e_theta = e_theta / norm(e_theta)
                do j=1,3
                    basis(j,1,i) = e_phi(j)
                    basis(j,2,i) = e_theta(j)
                end do
            end do

        end subroutine create_tangentspace_spherical_basis

        subroutine create_tangentspace_basis_full(spin, basis)
            !> Create a tangent space basis for a given spin configuration by gram-schmidt-orthogonalization. This is the
            !> updated version of the old code updated by H.S. October 2021. Previously it was build on a probilistic
            !> calculation. Also this basis stores only the blocks required. Since the tanget space is 2N dimensional we
            !> only need two basis vectors (xi, eta) to span the tangent space of each spin. Therefore the dimension of
            !> the basis is given by (3,2,N_atom)
            !================================Input Variable=============================================================
            real(kind=xp), intent(in) :: spin(:)                                !< Spin configuration
            !================================Output Variable============================================================
            real(kind=xp), intent(out) :: basis(:,:,:)                          !< basis for each spin (3x2 block)
            !================================Local Variable=============================================================
            integer :: i,j
            real(kind=xp) :: xi(3), eta(3)                                      !< basis vectors for each spin
            basis = 0.0d0
            do i=1,N_atom
                if (abs(spin(3*i))>=0.5d0) then
                    !Gram Schmidt for v=(1,0,0) explicitely written:
                    xi = (/1.0-spin(3*i-2)*spin(3*i-2), -spin(3*i-2)*spin(3*i-1), -spin(3*i-2)*spin(3*i)/)
                else
                    !Gram Schmidt for v=(0,0,1) explicitely written:
                    xi = (/-spin(3*i)*spin(3*i-2), -spin(3*i)*spin(3*i-1), 1.0-spin(3*i)*spin(3*i)/)
                end if
                xi = xi / norm(xi)
                eta = cross(xi,spin(3*i-2:3*i))
                eta = eta / norm(eta)
                do j=1,3
                    basis(j,1,i) = xi(j)
                    basis(j,2,i) = eta(j)
                end do
            end do
        end subroutine create_tangentspace_basis_full

        subroutine create_tangentspace_basis_subsys(spin, idx_list, basis)
            !> Creates tangent space basis in subsystem based on index list
            !===================================Input Variable==========================================================
            real(kind=xp), intent(in) :: spin(:)                        !< Full spin configuration
            integer, intent(in) :: idx_list(:)                          !< Indices of subsystem atoms
            !===================================Output Variable====================================== ===================
            real(kind=xp), intent(out) :: basis(:,:,:)                  !< Basis for sub-system
            !===================================Local Variable==========================================================
            real(kind=xp) :: xi(3), eta(3)
            integer :: i, j
            basis = 0.0d0
            do i=1,size(idx_list)
                if (abs(spin(3*idx_list(i)))>=0.5d0) then
                    !Gram Schmidt for v=(1,0,0) explicitely written:
                    xi = (/1.0-spin(3*idx_list(i)-2)*spin(3*idx_list(i)-2), &
                            & -spin(3*idx_list(i)-2)*spin(3*idx_list(i)-1), -spin(3*idx_list(i)-2)*spin(3*idx_list(i))/)
                else
                    !Gram Schmidt for v=(0,0,1) explicitely written:
                    xi = (/-spin(3*idx_list(i))*spin(3*idx_list(i)-2), -spin(3*idx_list(i))*spin(3*idx_list(i)-1), &
                            & 1.0-spin(3*idx_list(i))*spin(3*idx_list(i))/)
                end if
                xi = xi / norm(xi)
                eta = cross(xi,spin(3*idx_list(i)-2:3*idx_list(i)))
                eta = eta / norm(eta)
                do j=1,3
                    basis(j,1,i) = xi(j)
                    basis(j,2,i) = eta(j)
                end do
            end do
        end subroutine create_tangentspace_basis_subsys


        subroutine project_3x3_to_2x2_block(hessian3N_block,basis_i,basis_j,hessian2N_block)
            !> Calculates the projection of the 3x3 hessian element block ij to 2N space
            !================================Input Variable=============================================================
            real(kind=xp),intent(in) :: hessian3N_block(3,3)            !< Input 3N hessian
            real(kind=xp),intent(in) :: basis_i(3,2)                    !< Basis of spin i
            real(kind=xp),intent(in) :: basis_j(3,2)                    !< Basis of spin j
            !================================Inout Variable=============================================================
            real(kind=xp),intent(inout) :: hessian2N_block(2,2)         !< Inout 2N hessian
            !================================Local Variable=============================================================
            hessian2N_block(1,1) = basis_i(1,1) * (hessian3N_block(1,1)*basis_j(1,1) &
                                                & + hessian3N_block(1,2)*basis_j(2,1) &
                                                & + hessian3N_block(1,3)*basis_j(3,1)) &
                                & + basis_i(2,1) * (hessian3N_block(2,1)*basis_j(1,1) &
                                                & + hessian3N_block(2,2)*basis_j(2,1) &
                                                & + hessian3N_block(2,3)*basis_j(3,1)) &
                                & + basis_i(3,1) * (hessian3N_block(3,1)*basis_j(1,1) &
                                                & + hessian3N_block(3,2)*basis_j(2,1) &
                                                & +hessian3N_block(3,3)*basis_j(3,1))
            hessian2N_block(1,2) = basis_i(1,1) * (hessian3N_block(1,1)*basis_j(1,2) &
                                                & + hessian3N_block(1,2)*basis_j(2,2) &
                                                & + hessian3N_block(1,3)*basis_j(3,2)) &
                                & + basis_i(2,1) * (hessian3N_block(2,1)*basis_j(1,2) &
                                                & + hessian3N_block(2,2)*basis_j(2,2) &
                                                & + hessian3N_block(2,3)*basis_j(3,2)) &
                                & + basis_i(3,1) * (hessian3N_block(3,1)*basis_j(1,2) &
                                                & + hessian3N_block(3,2)*basis_j(2,2) &
                                                & + hessian3N_block(3,3)*basis_j(3,2))
            hessian2N_block(2,1) = basis_i(1,2) * (hessian3N_block(1,1)*basis_j(1,1) &
                                                & + hessian3N_block(1,2)*basis_j(2,1) &
                                                & + hessian3N_block(1,3)*basis_j(3,1)) &
                                & + basis_i(2,2) * (hessian3N_block(2,1)*basis_j(1,1) &
                                                & + hessian3N_block(2,2)*basis_j(2,1) &
                                                & + hessian3N_block(2,3)*basis_j(3,1)) &
                                & + basis_i(3,2) * (hessian3N_block(3,1)*basis_j(1,1) &
                                                & + hessian3N_block(3,2)*basis_j(2,1) &
                                                & + hessian3N_block(3,3)*basis_j(3,1))
            hessian2N_block(2,2) = basis_i(1,2) * (hessian3N_block(1,1)*basis_j(1,2) &
                                                & + hessian3N_block(1,2)*basis_j(2,2) &
                                                & + hessian3N_block(1,3)*basis_j(3,2)) &
                                & + basis_i(2,2) * (hessian3N_block(2,1)*basis_j(1,2) &
                                                & + hessian3N_block(2,2)*basis_j(2,2) &
                                                & + hessian3N_block(2,3)*basis_j(3,2)) &
                                & + basis_i(3,2) * (hessian3N_block(3,1)*basis_j(1,2) &
                                                & + hessian3N_block(3,2)*basis_j(2,2) &
                                                & + hessian3N_block(3,3)*basis_j(3,2))
        end subroutine project_3x3_to_2x2_block

        subroutine project_hessian3N_bilinear_dense(hessian3N,sparsity_pattern3N, spin, basis, gradient, hessian2N)
            !> Projects the hessian in 3N space towards the tangent space representations. The 3N-hessian is here given
            !> in dense matrix formulation. The projection is performed as follows:
            !> The matrix block H3N_ij is transformed to the block H2N_ij by:
            !> H2N_ij = U_i^T (H3N_ij - Shape_i) U_j
            !> Here Shape is the Shape Operator accounting for the curvature of the 2N configuration space. It is 3Nx3N
            !> diagonal matrix with Shape_i = dot(s_i,grad_i) * unity3x3 on each of the three diagonal entries of the block ii.
            !> Afterwards the projections are applied:
            !> U_i is a 3x2 matrix with
            !>       |xi_i^x eta_i^x|
            !> U_i = |xi_i^y eta_i^x| and U_j is the corresponding transposed matrix for basis of spin j
            !>       |xi_i^z eta_i^x|
            !> This can also be written for the whole matrix H_2N = U^T (H_3N - Shape ) U with U beeing block-diagonal
            !> and is 3N x 2N.

            !> This method now does NOT perform the whole 3N x 3N matrix multiplications. Instead it just uses the ex.
            !> blocks to save time. Therefore we again need the adjacency information.
            !================================Input Variable=============================================================
            real(kind=xp), intent(in) :: hessian3N(3*N_atom,3*N_atom)
            integer, intent(in) :: sparsity_pattern3N(N_atom,max_totalneigh_allbds+1)
            real(kind=xp), intent(in) :: basis(3,2,N_atom)
            real(kind=xp), intent(in) :: gradient(3*N_atom)
            real(kind=xp), intent(in) :: spin(3*N_atom)
            !================================Output Variable============================================================
            real(kind=xp), intent(out) :: hessian2N(2*N_atom, 2*N_atom)
            !================================Local Variable=============================================================
            integer :: i, j, k_spars_count
            real(kind=xp) :: lambda

            hessian2N = 0.0d0

            do i=1,N_atom
                ! Determine the js
                do k_spars_count = 1,max_totalneigh_allbds + 1
                    j = sparsity_pattern3N(i,k_spars_count)
                    ! If j=0 we reached the maximum of the neighbors and we can exit the loop
                    if (j==0) exit
                    if (j==i) then
                        ! Add diagonal block
                        ! Calculate the shape operator
                        lambda = -1.0d0 * (gradient(3*i-2) * spin(3*i-2) &
                                      & + gradient(3*i-1) * spin(3*i-1) &
                                      & + gradient(3*i) * spin(3*i))
                        hessian2N(2*i-1,2*j-1) = basis(1,1,i) * ((hessian3N(3*i-2,3*j-2)+lambda)*basis(1,1,j) &
                                                                & + hessian3N(3*i-2,3*j-1)*basis(2,1,j) &
                                                                & + hessian3N(3*i-2,3*j)*basis(3,1,j)) &
                                             & + basis(2,1,i) * (hessian3N(3*i-1,3*j-2)*basis(1,1,j) &
                                                                & + (hessian3N(3*i-1,3*j-1)+lambda)*basis(2,1,j) &
                                                                & + hessian3N(3*i-1,3*j)*basis(3,1,j)) &
                                             & + basis(3,1,i) * (hessian3N(3*i,3*j-2)*basis(1,1,j) &
                                                                & + hessian3N(3*i,3*j-1)*basis(2,1,j) &
                                                                & +(hessian3N(3*i,3*j)+lambda)*basis(3,1,j))
                        hessian2N(2*i-1,2*j) = basis(1,1,i) * ((hessian3N(3*i-2,3*j-2)+lambda)*basis(1,2,j) &
                                                                & + hessian3N(3*i-2,3*j-1)*basis(2,2,j) &
                                                                & + hessian3N(3*i-2,3*j)*basis(3,2,j)) &
                                           & + basis(2,1,i) * (hessian3N(3*i-1,3*j-2)*basis(1,2,j) &
                                                                & + (hessian3N(3*i-1,3*j-1)+lambda)*basis(2,2,j) &
                                                                & + hessian3N(3*i-1,3*j)*basis(3,2,j)) &
                                           & + basis(3,1,i) * (hessian3N(3*i,3*j-2)*basis(1,2,j) &
                                                                & + hessian3N(3*i,3*j-1)*basis(2,2,j) &
                                                                & + (hessian3N(3*i,3*j)+lambda)*basis(3,2,j))
                        hessian2N(2*i,2*j-1) = basis(1,2,i) * ((hessian3N(3*i-2,3*j-2)+lambda)*basis(1,1,j) &
                                                                & + hessian3N(3*i-2,3*j-1)*basis(2,1,j) &
                                                                & + hessian3N(3*i-2,3*j)*basis(3,1,j)) &
                                           & + basis(2,2,i) * (hessian3N(3*i-1,3*j-2)*basis(1,1,j) &
                                                                & + (hessian3N(3*i-1,3*j-1)+lambda)*basis(2,1,j) &
                                                                & + hessian3N(3*i-1,3*j)*basis(3,1,j)) &
                                           & + basis(3,2,i) * (hessian3N(3*i,3*j-2)*basis(1,1,j) &
                                                                & + hessian3N(3*i,3*j-1)*basis(2,1,j) &
                                                                & + (hessian3N(3*i,3*j)+lambda)*basis(3,1,j))
                        hessian2N(2*i,2*j) = basis(1,2,i) * ((hessian3N(3*i-2,3*j-2)+lambda)*basis(1,2,j) &
                                                                & + hessian3N(3*i-2,3*j-1)*basis(2,2,j) &
                                                                & + hessian3N(3*i-2,3*j)*basis(3,2,j)) &
                                         & + basis(2,2,i) * (hessian3N(3*i-1,3*j-2)*basis(1,2,j) &
                                                                & + (hessian3N(3*i-1,3*j-1)+lambda)*basis(2,2,j) &
                                                                & + hessian3N(3*i-1,3*j)*basis(3,2,j)) &
                                         & + basis(3,2,i) * (hessian3N(3*i,3*j-2)*basis(1,2,j) &
                                                                & + hessian3N(3*i,3*j-1)*basis(2,2,j) &
                                                                & + (hessian3N(3*i,3*j)+lambda)*basis(3,2,j))
                    else
                        ! Add offdiagonal block
                        hessian2N(2*i-1,2*j-1) = basis(1,1,i) * ((hessian3N(3*i-2,3*j-2))*basis(1,1,j) &
                                                                & + hessian3N(3*i-2,3*j-1)*basis(2,1,j) &
                                                                & + hessian3N(3*i-2,3*j)*basis(3,1,j)) &
                                             & + basis(2,1,i) * (hessian3N(3*i-1,3*j-2)*basis(1,1,j) &
                                                                & + (hessian3N(3*i-1,3*j-1))*basis(2,1,j) &
                                                                & + hessian3N(3*i-1,3*j)*basis(3,1,j)) &
                                             & + basis(3,1,i) * (hessian3N(3*i,3*j-2)*basis(1,1,j) &
                                                                & + hessian3N(3*i,3*j-1)*basis(2,1,j) &
                                                                & + (hessian3N(3*i,3*j))*basis(3,1,j))
                        hessian2N(2*j-1,2*i-1) = hessian2N(2*i-1,2*j-1)
                        hessian2N(2*i-1,2*j) = basis(1,1,i) * ((hessian3N(3*i-2,3*j-2))*basis(1,2,j) &
                                                                & + hessian3N(3*i-2,3*j-1)*basis(2,2,j) &
                                                                & + hessian3N(3*i-2,3*j)*basis(3,2,j)) &
                                           & + basis(2,1,i) * (hessian3N(3*i-1,3*j-2)*basis(1,2,j) &
                                                                & + (hessian3N(3*i-1,3*j-1))*basis(2,2,j) &
                                                                & + hessian3N(3*i-1,3*j)*basis(3,2,j)) &
                                           & + basis(3,1,i) * (hessian3N(3*i,3*j-2)*basis(1,2,j) &
                                                                & + hessian3N(3*i,3*j-1)*basis(2,2,j) &
                                                                & + (hessian3N(3*i,3*j))*basis(3,2,j))
                        hessian2N(2*j,2*i-1) = hessian2N(2*i-1,2*j)
                        hessian2N(2*i,2*j-1) = basis(1,2,i) * ((hessian3N(3*i-2,3*j-2))*basis(1,1,j) &
                                                                & + hessian3N(3*i-2,3*j-1)*basis(2,1,j) &
                                                                & + hessian3N(3*i-2,3*j)*basis(3,1,j)) &
                                           & + basis(2,2,i) * (hessian3N(3*i-1,3*j-2)*basis(1,1,j) &
                                                                & + (hessian3N(3*i-1,3*j-1))*basis(2,1,j) &
                                                                & + hessian3N(3*i-1,3*j)*basis(3,1,j)) &
                                           & + basis(3,2,i) * (hessian3N(3*i,3*j-2)*basis(1,1,j) &
                                                                & + hessian3N(3*i,3*j-1)*basis(2,1,j) &
                                                                & + (hessian3N(3*i,3*j))*basis(3,1,j))
                        hessian2N(2*j-1,2*i) = hessian2N(2*i,2*j-1)
                        hessian2N(2*i,2*j) = basis(1,2,i) * ((hessian3N(3*i-2,3*j-2))*basis(1,2,j) &
                                                                & + hessian3N(3*i-2,3*j-1)*basis(2,2,j) &
                                                                & + hessian3N(3*i-2,3*j)*basis(3,2,j)) &
                                           & + basis(2,2,i) * (hessian3N(3*i-1,3*j-2)*basis(1,2,j) &
                                                                & + (hessian3N(3*i-1,3*j-1))*basis(2,2,j) &
                                                                & + hessian3N(3*i-1,3*j)*basis(3,2,j)) &
                                           & + basis(3,2,i) * (hessian3N(3*i,3*j-2)*basis(1,2,j) &
                                                                & + hessian3N(3*i,3*j-1)*basis(2,2,j) &
                                                                & + (hessian3N(3*i,3*j))*basis(3,2,j))
                        hessian2N(2*j,2*i) = hessian2N(2*i,2*j)
                    end if
                end do
            end do
        end subroutine project_hessian3N_bilinear_dense

        subroutine project_hessian3N_hessian2N_dense_all_elements(hessian3N,basis,shape_op,hessian2N,interface_char)
            !> Projects the dense 3N hessian by subtracting the shape operator and applying the projection matrices U_i
            !> and U_j (tangent space basis vectors). This is the naive implementation where we run over all elements in
            !> a row and only exploit the symmetry of the matrix. Not efficient at all
            !===================================Input Variable==========================================================
            real(kind=xp), intent(in) :: hessian3N(:,:)                         !< Input 3N Hessian
            real(kind=xp), intent(in) :: basis(:,:,:)                           !< Input 3x2 Projection Matrices
            real(kind=xp), intent(in) :: shape_op(:)                            !< Shape Operator
            character, intent(in) :: interface_char                             !< Character for interface decision
            !===================================Output Variable=========================================================
            real(kind=xp), intent(out) :: hessian2N(:,:)                        !< Output 2N hessian
            !===================================Local Variable==========================================================
            integer :: i,j                                                      !< Iterators for atoms
            real(kind=xp) :: dummy3x3_block(3,3)                                !< Dummy blocks for simplyfie notations
            real(kind=xp) :: dummy2x2_block(2,2)

            hessian2N = 0.0d0

            do i=1,N_atom
                do j=i,N_atom
                    if (i==j) then
                        ! Diagonal Elements
                        dummy3x3_block(1,1) = hessian3N(3*i-2,3*i-2) - shape_op(i)
                        dummy3x3_block(1,2) = hessian3N(3*i-2,3*i-1)
                        dummy3x3_block(1,3) = hessian3N(3*i-2,3*i)
                        dummy3x3_block(2,1) = hessian3N(3*i-1,3*i-2)
                        dummy3x3_block(2,2) = hessian3N(3*i-1,3*i-1) - shape_op(i)
                        dummy3x3_block(2,3) = hessian3N(3*i-1,3*i)
                        dummy3x3_block(3,1) = hessian3N(3*i,3*i-2)
                        dummy3x3_block(3,2) = hessian3N(3*i,3*i-1)
                        dummy3x3_block(3,3) = hessian3N(3*i,3*i) - shape_op(i)
                        dummy2x2_block = 0.0d0
                        call project_3x3_to_2x2_block(dummy3x3_block,basis(:,:,i),basis(:,:,i),dummy2x2_block)
                        hessian2N(2*i-1:2*i,2*i-1:2*i) = dummy2x2_block
                    else
                        ! Off diagonal Elements
                        dummy3x3_block(1,1) = hessian3N(3*i-2,3*j-2)
                        dummy3x3_block(1,2) = hessian3N(3*i-2,3*j-1)
                        dummy3x3_block(1,3) = hessian3N(3*i-2,3*j)
                        dummy3x3_block(2,1) = hessian3N(3*i-1,3*j-2)
                        dummy3x3_block(2,2) = hessian3N(3*i-1,3*j-1)
                        dummy3x3_block(2,3) = hessian3N(3*i-1,3*j)
                        dummy3x3_block(3,1) = hessian3N(3*i,3*j-2)
                        dummy3x3_block(3,2) = hessian3N(3*i,3*j-1)
                        dummy3x3_block(3,3) = hessian3N(3*i,3*j)
                        dummy2x2_block = 0.0d0
                        call project_3x3_to_2x2_block(dummy3x3_block,basis(:,:,i),basis(:,:,j),dummy2x2_block)
                        hessian2N(2*i-1:2*i,2*j-1:2*j) = dummy2x2_block
                        hessian2N(2*j-1,2*i-1) = hessian2N(2*i-1,2*j-1)
                        hessian2N(2*j-1,2*i) = hessian2N(2*i,2*j-1)
                        hessian2N(2*j,2*i-1) = hessian2N(2*i-1,2*j)
                        hessian2N(2*j,2*i) = hessian2N(2*i,2*j)
                    end if
                end do
            end do

        end subroutine project_hessian3N_hessian2N_dense_all_elements

        subroutine project_hessian3N_hessian2N_dense_neighbors(hessian3N,basis,shape_op,hessian2N)
            !> Projects the dense 3N hessian by subtracting the shape operator and applying the projection matrices U_i
            !> and U_j (tangent space basis vectors). This is the slightly advanced method vs just calculating the
            !> matrix products for each element in the row by computing the products only for the neighbors listed in
            !> the neighbor table.
            !===================================Input Variable==========================================================
            real(kind=xp), intent(in) :: hessian3N(:,:)                         !< Input 3N Hessian
            real(kind=xp), intent(in) :: basis(:,:,:)                           !< Input 3x2 Projection Matrices
            real(kind=xp), intent(in) :: shape_op(:)                            !< Shape Operator
            !===================================Output Variable=========================================================
            real(kind=xp), intent(out) :: hessian2N(:,:)                        !< Output 2N hessian
            !===================================Local Variable==========================================================
            integer :: i,j                                                      !< Iterators for atoms
            integer :: k_neigh                                                  !< Neighbor iterator
            real(kind=xp) :: dummy3x3_block(3,3)                                !< Dummy blocks for simplyfie notations
            real(kind=xp) :: dummy2x2_block(2,2)
            integer :: shape_neigh(2)                                           !< Shape of the neighor list
            hessian2N = 0.0d0

            do i=1,N_atom
                ! Diagonal Elements
                dummy3x3_block(1,1) = hessian3N(3*i-2,3*i-2) - shape_op(i)
                dummy3x3_block(1,2) = hessian3N(3*i-2,3*i-1)
                dummy3x3_block(1,3) = hessian3N(3*i-2,3*i)
                dummy3x3_block(2,1) = hessian3N(3*i-1,3*i-2)
                dummy3x3_block(2,2) = hessian3N(3*i-1,3*i-1) - shape_op(i)
                dummy3x3_block(2,3) = hessian3N(3*i-1,3*i)
                dummy3x3_block(3,1) = hessian3N(3*i,3*i-2)
                dummy3x3_block(3,2) = hessian3N(3*i,3*i-1)
                dummy3x3_block(3,3) = hessian3N(3*i,3*i) - shape_op(i)
                dummy2x2_block = 0.0d0
                call project_3x3_to_2x2_block(dummy3x3_block,basis(:,:,i),basis(:,:,i),dummy2x2_block)
                hessian2N(2*i-1:2*i,2*i-1:2*i) = dummy2x2_block
                ! Off diagonal Elements
                shape_neigh = shape(bilinear_neigh_allgrp_allbds)
                do k_neigh=1,shape_neigh(2)
                    j = bilinear_neigh_allgrp_allbds(i, k_neigh)
                    if (j==0) then
                        cycle
                    end if
                    dummy3x3_block(1,1) = hessian3N(3*i-2,3*j-2)
                    dummy3x3_block(1,2) = hessian3N(3*i-2,3*j-1)
                    dummy3x3_block(1,3) = hessian3N(3*i-2,3*j)
                    dummy3x3_block(2,1) = hessian3N(3*i-1,3*j-2)
                    dummy3x3_block(2,2) = hessian3N(3*i-1,3*j-1)
                    dummy3x3_block(2,3) = hessian3N(3*i-1,3*j)
                    dummy3x3_block(3,1) = hessian3N(3*i,3*j-2)
                    dummy3x3_block(3,2) = hessian3N(3*i,3*j-1)
                    dummy3x3_block(3,3) = hessian3N(3*i,3*j)
                    dummy2x2_block = 0.0d0
                    call project_3x3_to_2x2_block(dummy3x3_block,basis(:,:,i),basis(:,:,j),dummy2x2_block)
                    hessian2N(2*i-1:2*i,2*j-1:2*j) = dummy2x2_block
                end do
            end do

        end subroutine project_hessian3N_hessian2N_dense_neighbors

        subroutine project_hessian3N_hessian2N_subsystem_dense_neighbors(hessian3N,idx_list,basis,shape_op,hessian2N)
            !> Projects the dense 3N hessian by subtracting the shape operator and applying the projection matrices U_i
            !> and U_j (tangent space basis vectors). This is the slightly advanced method vs just calculating the
            !> matrix products for each element in the row by computing the products only for the neighbors listed in
            !> the neighbor table. Does that only for atoms in the subsystem.
            !===================================Input Variable==========================================================
            real(kind=xp), intent(in) :: hessian3N(:,:)                         !< Input 3N Hessian of subsystem
            real(kind=xp), intent(in) :: basis(:,:,:)                           !< Input 3x2 Projection Matrices
            real(kind=xp), intent(in) :: shape_op(:)                            !< Shape Operator for subsystem atoms
            integer, intent(in) :: idx_list(:)                                  !< Subsystem List
            !===================================Output Variable=========================================================
            real(kind=xp), intent(out) :: hessian2N(:,:)                        !< Output 2N hessian for subsystem
            !===================================Local Variable==========================================================
            integer :: i,j                                                      !< Iterators for atoms
            integer :: k_neigh                                                  !< Neighbor iterator
            real(kind=xp) :: dummy3x3_block(3,3)                                !< Dummy blocks for simplify notations
            real(kind=xp) :: dummy2x2_block(2,2)
            integer :: shape_neigh(2)                                           !< Shape of the neighor list
            integer :: l_neigh_index_subsys
            hessian2N = 0.0d0

            do i=1,size(idx_list)
                ! Diagonal Elements
                dummy3x3_block(1,1) = hessian3N(3*i-2,3*i-2) - shape_op(i)
                dummy3x3_block(1,2) = hessian3N(3*i-2,3*i-1)
                dummy3x3_block(1,3) = hessian3N(3*i-2,3*i)
                dummy3x3_block(2,1) = hessian3N(3*i-1,3*i-2)
                dummy3x3_block(2,2) = hessian3N(3*i-1,3*i-1) - shape_op(i)
                dummy3x3_block(2,3) = hessian3N(3*i-1,3*i)
                dummy3x3_block(3,1) = hessian3N(3*i,3*i-2)
                dummy3x3_block(3,2) = hessian3N(3*i,3*i-1)
                dummy3x3_block(3,3) = hessian3N(3*i,3*i) - shape_op(i)
                dummy2x2_block = 0.0d0
                call project_3x3_to_2x2_block(dummy3x3_block,basis(:,:,i),basis(:,:,i),dummy2x2_block)
                hessian2N(2*i-1:2*i,2*i-1:2*i) = dummy2x2_block
                ! Off diagonal Elements
                shape_neigh = shape(bilinear_neigh_allgrp_allbds)
                do k_neigh=1,shape_neigh(2)
                    j = bilinear_neigh_allgrp_allbds(idx_list(i), k_neigh)
                    if (.not.(is_in(idx_list,j,l_neigh_index_subsys))) cycle
                    j = l_neigh_index_subsys
                    dummy3x3_block(1,1) = hessian3N(3*i-2,3*j-2)
                    dummy3x3_block(1,2) = hessian3N(3*i-2,3*j-1)
                    dummy3x3_block(1,3) = hessian3N(3*i-2,3*j)
                    dummy3x3_block(2,1) = hessian3N(3*i-1,3*j-2)
                    dummy3x3_block(2,2) = hessian3N(3*i-1,3*j-1)
                    dummy3x3_block(2,3) = hessian3N(3*i-1,3*j)
                    dummy3x3_block(3,1) = hessian3N(3*i,3*j-2)
                    dummy3x3_block(3,2) = hessian3N(3*i,3*j-1)
                    dummy3x3_block(3,3) = hessian3N(3*i,3*j)
                    dummy2x2_block = 0.0d0
                    call project_3x3_to_2x2_block(dummy3x3_block,basis(:,:,i),basis(:,:,j),dummy2x2_block)
                    hessian2N(2*i-1:2*i,2*j-1:2*j) = dummy2x2_block
                end do
            end do

        end subroutine project_hessian3N_hessian2N_subsystem_dense_neighbors

        subroutine to_tangentspace_grassmann_manifold_2D(Z,Y,Z_tspace)
            !> The orthogonal projection onto the tangent space of Gr(n,p) is given by
            !> Z_tspace = (I-YY^T)*Z
            !> Z_tspace = Z - Y (Y^T * Z) = Z- Y * C
            !===================================Input Variable==========================================================
            real(kind=xp), intent(in) :: Z(:,:)                             !< (n x p) matrix to be projected
            real(kind=xp), intent(in) :: Y(:,:)                             !< (n x p) matrix repr. equiv. class of mat.
                                                                            !< corresponding to a point on Gr(n,p)
            !===================================Output Variable=========================================================
            real(kind=xp), intent(out) :: Z_tspace(:,:)                     !< projection on Gr(n,p) at Y
            !===================================Local Variable==========================================================
            integer :: n,p,shapeY(2)                                        !< Matrix dimension
            real(kind=xp), allocatable :: C(:,:)                            !< Product Y^T * Z
            real(kind=xp), allocatable :: dummyvec(:)
            integer :: i_p, i_p2
            shapeY=shape(Y)
            n = shapeY(1)
            p = shapeY(2)
            allocate(dummyvec(n))
            allocate(C(p,p))
            Z_tspace = 0.0d0
            do i_p=1,p
                do i_p2=1,p
                    C(i_p,i_p2) = DOT_PRODUCT(Y(:,i_p),Z(:,i_p2))
                end do
            end do
            do i_p=1,p
                dummyvec = 0.0d0
                do i_p2=1,p
                    dummyvec = dummyvec +  Y(:,i_p2) * C(i_p2,i_p)
                end do
                Z_tspace(:,i_p) = Z(:,i_p) - dummyvec
            end do
        end subroutine to_tangentspace_grassmann_manifold_2D

        subroutine to_tangentspace_grassmann_manifold_1D(Z,Y)
            !> The orthogonal projection onto the tangent space of Gr(n,p) is given by
            !> Z_tspace = (I-YY^T)*Z
            !> Z_tspace = Z - Y (Y^T * Z) = Z- Y * C
            !===================================Input Variable==========================================================
            real(kind=xp), intent(in) :: Y(:,:)                             !< (n x p) matrix repr. equiv. class of mat.
                                                                            !< corresponding to a point on Gr(n,p)
            !===================================Inout Variable==========================================================
            real(kind=xp), intent(inout) :: Z(:)                            !< (n x p) matrix to be projected,
                                                                            !< represented as n*p array
            !===================================Local Variable==========================================================
            integer :: n,p,shapeY(2)                                        !< Matrix dimension
            real(kind=xp), allocatable :: C(:,:)                            !< Product Y^T * Z
            real(kind=xp), allocatable :: dummyvec(:)
            integer :: i_p, i_p2
            shapeY=shape(Y)
            n = shapeY(1)
            p = shapeY(2)
            allocate(dummyvec(n))
            allocate(C(p,p))
            do i_p=1,p
                do i_p2=1,p
                    C(i_p,i_p2) = DOT_PRODUCT(Y(:,i_p),Z(n*(i_p2-1)+1:n*i_p2))
                end do
            end do
            do i_p=1,p
                dummyvec = 0.0d0
                do i_p2=1,p
                    dummyvec = dummyvec +  Y(:,i_p2) * C(i_p2,i_p)
                end do
                Z(n*(i_p-1)+1:n*i_p) = Z(n*(i_p-1)+1:n*i_p) - dummyvec
            end do
        end subroutine to_tangentspace_grassmann_manifold_1D

        subroutine to_tangentspace_stiefel_manifold(Z,Y,Z_tspace)
            !> Projects any n x p Matrix (e.g. in Rayleigh Quotient Optimization the invariant subspace with eigenvectors
            !> of length n (probably 3*N_atom) and subspace dimension p) to the Tangent space of the Stiefel manifold at
            !> Point Y. The theory is described in
            !> "THE GEOMETRY OF ALGORITHMS WITH ORTHOGONALITY CONSTRAINTS"
            !> by ALAN EDELMAN, TOMA ́S A. ARIAS, AND STEVEN T. SMITH
            !> The formula is given in equation 2.4: Z_tspace = Y skew(Y^T * Z) + (I-YY^T) * Z
            !> Since p<<n the above formula is reformulated as Y * skew(Y^T * Z) + Z - Y * (Y^T * Z)
            !>                                                 Y * (skew(Y^T * Z) - Y^T * Z) + Z
            !>                                                 Y * K + Z
            !> with skew(A)= (A-A^T)/2
            !===================================Input Variable==========================================================
            real(kind=xp), intent(in) :: Z(:,:)                     !< n x p matrix which shall be projected to tangent
                                                                    !< space of Stiefel Manifold at a certain point Y
            real(kind=xp), intent(in) :: Y(:,:)                     !< Point on Stiefel Manifold. A point on a Stiefel
                                                                    !< manifold is defined by a nxp Matrix with ortho-
                                                                    !< normal columns
            !===================================Output Variable=========================================================
            real(kind=xp), intent(out) :: Z_tspace(:,:)             !< n x p matrix in Tangent Space of Stiefel Manifold
            !===================================Local Variable==========================================================
            integer :: p,n                                          !< Dimensions
            integer :: i_p, i2_p, i_n                               !< Iterators
            real(Kind=xp), allocatable :: K(:,:)                    !< skewY^T * Z + Y^T * Z
            real(kind=xp) :: dummyfloat                             !< Dummy Float
            real(kind=xp), allocatable :: dummymat(:,:)             !< Dummy Matrix
            !---------Comparison for debug purposes----------
            real(kind=xp), allocatable :: YTZ(:,:), skewYTZ(:,:), YskewYTZ(:,:), I(:,:),YYT(:,:),I_YYT(:,:),I_YYTZ(:,:), YT(:,:)
            real(kind=xp), allocatable :: comparison_result(:,:)
            real(kind=xp), allocatable :: ZT(:,:),ZTY(:,:)

            !real(kind=xp) :: ini_time, fin_time

            p = size(Z,DIM=2)
            n = size(Z,DIM=1)
            allocate(K(p,p))
            allocate(dummymat(n,p))
            allocate(ZT(p,n),ZTY(p,p))
            !allocate(YT(p,n))
            !allocate(YTZ(p,p),skewYTZ(p,p),YskewYTZ(n,p),I(n,n),YYT(n,n),I_YYT(n,n),I_YYTZ(n,p),comparison_result(n,p))

            !write(*,*) "Calculation of skew (YTZ) and YTZ -> calculation of K"
            !ini_time = omp_get_wtime()

            do i_p=1,p
                do i2_p=1,p
                    K(i_p,i2_p) = (DOT_PRODUCT(Y(:,i_p),Z(:,i2_p)) - DOT_PRODUCT(Y(:,i2_p),Z(:,i_p)))/2.0d0 - &
                            & DOT_PRODUCT(Y(:,i_p),Z(:,i2_p))
                end do
            end do
            !fin_time = omp_get_wtime()
            !write(*,*) "time: ", fin_time - ini_time
            !write(*,*) "==============="

            !write(*,*) "Calculation of Y * K + Z"
            !ini_time = omp_get_wtime()
            do i_n=1,n
                do i_p=1,p
                    dummyfloat = 0.0d0
                    do i2_p=1,p
                        dummyfloat = dummyfloat + Y(i_n,i2_p) * K(i2_p,i_p)
                    end do
                    Z_tspace(i_n,i_p) = dummyfloat + Z(i_n,i_p)
                end do
            end do
            !Z_tspace = dummymat

            ! Test if Z_tspace is really in tangent space
            ! and if inner product is zero 0!=trace(Z^T * Y)
            do i_n=1,n
                do i_p=1,p
                    ZT(i_p,i_n) = Z_tspace(i_n,i_p)
                end do
            end do
            call AB(ZT, p, n, Y, p, ZTY)
            write(*,*) ZTY
            dummyfloat = 0.0d0
            do i_p=1,p
                dummyfloat = dummyfloat + ZTY(i_p,i_p)
            end do
            write(*,*) "Trace: ",dummyfloat


            !fin_time = omp_get_wtime()
            !write(*,*) "time: ", fin_time - ini_time
            if (.False.) then
                ! Y skew(Y^T * Z) + (I-YY^T) * Z
                ! n x p * skew(p x n * n x p) + (n x n + n x p * p x n) * n x p
                !YTZ(p,p)
                !skewYTZ(p,p)
                !YskewYTZ(n,p)
                !I(n,n)
                !YYT(n,n)
                !I_YYT(n,n)
                !I_YYTZ(n,p)
                write(*,*) "Calc of YT"
                do i_n=1,n
                    do i_p=1,p
                        YT(i_p,i_n) = Y(i_n,i_p)
                    end do
                end do
                write(*,*) "Calc of YTZ"
                call AB(YT, p, n, Z, p, YTZ)
                write(*,*) "Calc of skew YTZ:"
                do i_p=1,p
                    do i2_p=1,p
                        skewYTZ(i_p,i2_p) = (YTZ(i_p,i2_p)-YTZ(i2_p,i_p))/2.0d0
                    end do
                end do
                write(*,*) "Calc of Y * skew YTZ:"
                call AB(Y, n, p, skewYTZ, p, YskewYTZ)
                write(*,*) "Calc of identity"
                I=0.0d0
                do i_n = 1, n
                    I(i_n,i_n) = 1.0d0
                end do
                write(*,*) "Calc of YYT"
                call AB(Y, n, p, YT, n, YYT)
                write(*,*) "Calc of I - YYT"
                I_YYT = I - YYT
                write(*,*) "Calc of (I-YY^T) * Z"
                call AB(I_YYT, n, n, Z, p, I_YYTZ)
                write(*,*) "Calc of result"
                comparison_result = YskewYTZ + I_YYTZ

                write(*,*) "Comparison of result"
                do i_p = 1, p
                    write(*,*) norm(comparison_result(:,i_p) - Z_tspace(:,i_p))
                end do

            end if

        end subroutine to_tangentspace_stiefel_manifold
end module spk_tangentspace
