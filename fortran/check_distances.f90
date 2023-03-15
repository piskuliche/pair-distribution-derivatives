program check_distance
    ! This program checks the distance between molecules
    ! and prints out the molecule ID if it is less than a cutoff
    implicit none
    integer :: i, j, k, m
    integer :: nmols, nframes, natoms, ntypes
    integer :: ichunk, nchunks, frames_per_chunk
    integer :: file_test, dim
    
    real :: cutoff, drsq, dr



    real, allocatable :: L(:,:)
    real, allocatable :: r(:,:,:)
    real, allocatable :: mass(:)
    integer, allocatable :: neighbors(:,:,:)
    integer, allocatable :: neigh_count(:,:)

    character(len=1024) :: out_fname
    character(len=1024) :: format_string
    character(len=1024) :: runname, subdir
    character(len=20) :: cmd_line_char, mol_char
    
    if (command_argument_count() .ne. 1) then
        write(*,*) 'Error: Need to include cutoff distance on command line'
        stop
    endif

    call get_command_argument(1, cmd_line_char)
    read(cmd_line_char,*) cutoff

    write(*,*) "Cutoff for this run is ", cutoff

    
    !Read input file that provides L, nmols, nframes, and cutoff
    call read_input(nmols, natoms, nframes, frames_per_chunk, dim, runname)

    ! Write to subdir variable the path
    subdir = trim(runname) // '_' // trim(cmd_line_char)

    call system("mkdir -p " // trim(subdir) // "/neighbors")

    nchunks = nframes / frames_per_chunk
    
    
    !Read mass file that provides the masses and the types
    open(12, file="masses.dat", status="old")

    read(12,*) ntypes
    write(*,*) ntypes
    allocate(mass(ntypes))



    do i = 1, ntypes
        write(*,*) i
        read(12,*) mass(i)
    end do

    close(12)
    
    !Allocate memory for the coordinates
    allocate(r(frames_per_chunk,nmols,3))
    !Allocate the neighbor list
    allocate(neighbors(frames_per_chunk,nmols,1000))
    allocate(neigh_count(frames_per_chunk,nmols))
    allocate(L(frames_per_chunk,3))
    
    
    !Open the file
    open(21, file="dump.lammpsdump", status="old")

    do j = 1, nmols
        !write(out_fname, "(A9,A21,I0,A4)") adjustl(trim(runname)), "/neighbors/neighbors_", j, ".dat"
        write(mol_char,*) j
        out_fname = trim(subdir) // '/neighbors/neighbors_' // trim(adjustl(mol_char)) // trim('.dat')
        open(22, file=trim(out_fname))
        write(22, *) "# Neighbors for molecule ", j
        close(22)
    end do
    
    
    do ichunk = 1, nchunks
        write(*,*) "Reading chunk ", ichunk, " of ", nchunks

        do i = 1, frames_per_chunk
            call read_frame(21, natoms, nmols, ntypes, mass, r(i,:,:), L(i,:))
        end do
        

        neighbors = 0
        neigh_count = 0
        !$omp parallel do private(i,j,k,dr) shared(r, L, neighbors, neigh_count)
        do i = 1, frames_per_chunk
            do j = 1, nmols
                do k = j+1, nmols
                    !Check the distance between the molecules
                    dr = pbc_dist(pbc_components(r(i,j,:), r(i,k,:), L(i,:)))

                    if (dr < cutoff) then
                        neighbors(i,j, neigh_count(i,j) + 1) = k
                        neigh_count(i,j) = neigh_count(i,j) + 1
                        neighbors(i,k, neigh_count(i, k) + 1) = j
                        neigh_count(i,k) = neigh_count(i,k) + 1
                    end if

                end do ! End k loop
                
                !Ensure that at least one neighbor was found
                if (neigh_count(i,j) == 0) then
                    print *, "Error: No neighbors found for molecule ", j
                    stop
                end if
            end do ! End j loop
        end do ! End i loop
        !$omp end parallel do
        
        do i=1,frames_per_chunk
            ! Write the neighbors to a separate file for each molecule
            call write_neighbors(nmols, neighbors(i,:,:), neigh_count(i,:))
        end do ! End frames_per_chunk loop

    end do ! End ichunk loop
    
    close(21)
    write(*,*) "for a total of ", nframes, " frames"

    contains

        function pbc_components(r1, r2, L)
            ! This function calculates the distance between two points
            ! using the minimum image convention
            real, dimension(3), intent(in) :: r1, r2, L
            real, dimension(3) :: pbc_components

            pbc_components = r1 - r2 - L*anint((r1-r2)/L)

        end function pbc_components

        function pbc_dist(pbc_comp)
            ! This function calculates the distance between two points
            ! using the minimum image convention
            real, dimension(3), intent(in) :: pbc_comp
            real :: pbc_dist
            pbc_dist = sqrt(sum(pbc_comp**2))

        end function pbc_dist

        subroutine read_input(nmols, natoms, nframes, frames_per_chunk, dim, runname)
            ! This subroutine reads in the input file
            ! and sets the values of L, nmols, and nframes
            implicit none
            integer, intent(out) :: nmols, nframes
            integer, intent(out) :: dim, natoms, frames_per_chunk

            character(len=1024), intent(out) :: runname
            
            write(*,*) "Reading input file input.dat"
            open(1, file="input.dat", status="old")

            read(1,*)
            read(1,*) runname
            read(1,*)
            read(1,*) nmols, natoms
            read(1,*)
            read(1,*) nframes, frames_per_chunk
            read(1,*)
            read(1,*) dim
            close(1)

        end subroutine read_input

        subroutine read_frame(fileobj, natoms, nmols, ntypes, mass, rcom, L)
            ! This subroutine takes an already opened file and reads the next frame
            ! From this frame, it reads the coordinates of all of the atoms, as well is the molecule ids
            ! and it calculates the center of mass of each molecule, which it stores in the array rcom
            ! The file format for the frame has the following on each line: id, molid, typeid, q, x, y, z

            ! Args:
            !   fileobj (integer): The file object that is already opened
            !   natoms (integer): The number of atoms in the system
            !   nmols (integer): The number of molecules in the system
            !   ntypes (integer): The number of types of atoms in the system
            !   mass (1D array): The mass of each type of atom

            ! Returns:
            !   rcom (2D array): The center of mass of each molecule
            !   L (1D array): The box size
            
            implicit none
            integer, intent(in) :: nmols, natoms, ntypes, fileobj
            real, dimension(ntypes), intent(in) :: mass

            real, dimension(3), intent(out) :: L
            real, dimension(nmols,3), intent(out) :: rcom

            real, dimension(3) :: rtmp, b_low, b_high
            real, dimension(nmols,3) :: refcoord
            real, dimension(nmols) :: total_mass

            integer :: i, j, k
            integer :: id, typeid, molid
            real :: q, x, y, z

            ! Zero rcom
            rcom = 0.0
            total_mass = 0.0
            do i=1,5
                read(fileobj,*)
            end do 
            do i=1,3
                read(fileobj,*) b_low(i), b_high(i)
                L(i) = b_high(i) - b_low(i)
            end do
            read(fileobj,*)
            !Loop over the atoms
            do i = 1, natoms
                !Read in the atom information
                read(fileobj,*) id, molid, typeid, (rtmp(j), j = 1, 3)
                if (rcom(molid,1) == 0.0) then
                    do j = 1, 3
                        refcoord(molid,j) = rtmp(j)
                    end do
                else 
                    ! Wrap the x, y, z coordinate to the minimum image distance
                    rtmp(:) = refcoord(molid,:) + pbc_components(rtmp(:), refcoord(molid,:), L(:))
                end if
                ! Ensure that the molecule id is in the range 1 to nmols
                if (molid < 1 .or. molid > nmols) then
                    print *, "Error: Molecule id ", molid, " is out of range"
                    stop
                end if

                ! Ensure that the type id is in the range 1 to ntypes
                if (typeid < 1 .or. typeid > ntypes) then
                    print *, "Error: Type id ", typeid, " is out of range"
                    stop
                end if

                !Add the atom to the center of mass
                do j = 1, 3
                    rcom(molid,j) = rcom(molid,j) + rtmp(j)*mass(typeid)
                end do

                total_mass(molid) = total_mass(molid) + mass(typeid)
            end do
            !Divide by the total mass to get the center of mass
            do i = 1, nmols
                do j=1, 3
                    rcom(i,j) = rcom(i,j)/total_mass(i)
                end do
            end do

        end subroutine read_frame

        subroutine write_neighbors(nmols, neighbors, neigh_count)
            ! This subroutine writes the neighbors to a file
            implicit none
            integer, intent(in) :: nmols
            integer, dimension(nmols, 1000), intent(in) :: neighbors
            integer, dimension(nmols), intent(in) :: neigh_count

            integer :: i, j

            character(len=1024) :: out_fname
            character(len=20) :: mol_char

            do i = 1, nmols
                !write(out_fname, "(A9,A21,I0,A4)") adjustl(trim(runname)), "/neighbors/neighbors_", i, ".dat"
                write(mol_char,*) i
                out_fname = trim(subdir) // '/neighbors/neighbors_' // trim(adjustl(mol_char)) // trim('.dat')
                open(22, file=trim(out_fname), status='old', position='append')
                write(22, '(1000I5)') (int(neighbors(i, j)), j=1, neigh_count(i))
                close(22)
            end do

        end subroutine write_neighbors

end program check_distance
    




