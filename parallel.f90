! ----------------------------------------------------------------------------------%%
!
! parallel.f90
!
! this is the example i used to refresh myself on the basics of message passing.
! this program sums a vector by splitting it into chunks for each slave processor,
! each slave processor gives a partial sum back to the master processor, and a grand
! total sum is printed.
!
! mpif90 -I/usr/local/include -c parallel.f90
! mpif90 -I/usr/local/include -o parallel.o -L/usr/local/lib -lmpich
! mpirun -np 4 ./a.out
! ----------------------------------------------------------------------------------%%

program parallel
	
IMPLICIT NONE

include 'mpif.h'

! DECLARE MPI STUFF
!integer :: ierr
integer, parameter :: max_rows = 10000000
integer, parameter :: send_data_tag = 2001, return_data_tag = 2002
integer :: my_id, root_process, ierr, status(MPI_STATUS_SIZE)
integer :: num_procs, an_id, num_rows_to_receive
integer :: avg_rows_per_process, num_rows, num_rows_to_send
integer :: end_row, sender, start_row, num_rows_received
real :: vector(max_rows), vector2(max_rows), partial_sum, sum
integer :: i




! PLAYING WITH MESSAGE PASSING

! process #0 is the root process
root_process = 0

! initialize a process
call MPI_INIT ( ierr )

! find out the process ID and how many processes were started so far
call MPI_COMM_RANK (MPI_COMM_WORLD, my_id, ierr)
call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)

! what to do if the process is the root process
if (my_id .eq. root_process) then
	
	! put number of rows in vector here
	num_rows = 100
	avg_rows_per_process = num_rows / num_procs
	
	! initialize the vector
	do I = 1,num_rows
		vector(i) = float(i)
	end do
	
	! distribute a portion of the vector to each child process
	do an_id = 1, num_procs -1
        start_row = ( an_id * avg_rows_per_process) + 1
        end_row = start_row + avg_rows_per_process - 1
        if (an_id .eq. (num_procs - 1)) end_row = num_rows
        num_rows_to_send = end_row - start_row + 1
		
		! send size to the an_id
		! MPI_SEND(initialAddress, numElements, dataType, rankOfDest, msgTag, comHandle, ierr)
        call MPI_SEND( num_rows_to_send, 1, MPI_INT, &
		an_id, send_data_tag, MPI_COMM_WORLD, ierr)
		
		! send vector to the an_id
		! MPI_SEND(initialAddress, numElements, dataType, rankOfDest, msgTag, comHandle, ierr)
        call MPI_SEND( vector(start_row), num_rows_to_send, MPI_REAL, &
		an_id, send_data_tag, MPI_COMM_WORLD, ierr)
     end do
	 
	 ! calculate the sum of the values in the segment assigned to the root process
	 sum = 0.0
	 do i = 1, avg_rows_per_process
		 sum = sum + vector(i)
	 end do
	 
	 write(*,*) sum
	 write(*,*) "calculated by root process"
	 
	 ! finally, collect the partial sums from the slave processes, print them, then add them
	 ! to the grand sum, and print it
	 do an_id = 1, num_procs -1
		
		call MPI_RECV( partial_sum, 1, MPI_REAL, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
		sender = status(MPI_SOURCE)
		write(*,*) "partial sum:"
		write(*,*) partial_sum
		write(*,*) "from process:"
		write(*,*) sender
		write(*,*) " "

        sum = sum + partial_sum 
     end do
	 
	 write(*,*) "the GRAND TOTAL SUM:"
	 write(*,*) sum
	
! what to do if the process is a slave process 
else
	! here is a slave process, must receive vector segment and store it
	! in a local vector, vector2
	 
	! receive size
	! MPI_RECV(output, maxNumElements, dataType, rankOfSource, msgTag, comHandle, status, ierr)
	call MPI_RECV ( num_rows_to_receive, 1 , MPI_INT, &
	root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
	
	! receive chunk, save in local vector2
	call MPI_RECV ( vector2, num_rows_to_receive, MPI_REAL, &
	root_process, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
	num_rows_received = num_rows_to_receive

	! Calculate the sum of my portion of the vector,
	partial_sum = 0.0
	do i = 1, num_rows_received
		partial_sum = partial_sum + vector2(i)
	end do

	! and, finally, send my partial sum to the root process.
	call MPI_SEND( partial_sum, 1, MPI_REAL, root_process, &
	return_data_tag, MPI_COMM_WORLD, ierr)

endif
	 
call MPI_FINALIZE ( ierr )
!stop

end program parallel
