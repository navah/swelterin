# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# User specific aliases and functions
#
# set up the intel compilers
#. /share/apps/intel/bin/iccvars.sh intel64
#. /share/apps/intel/bin/ifortvars.sh intel64

# set up the PGI compilers
. /share/apps/pgi/linux86-64/10.9/pgi.sh

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/navah/local/lib
module load intel/intel-12
module load mpi/mvapich/intel
