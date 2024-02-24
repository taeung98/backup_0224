#!/bin/bash
. /opt/sge/default/common/settings.sh

module() { eval `/usr/bin/modulecmd bash $*`; }
export -f module
module use /opt/Modules/modulefiles

#module load gsl/gcc-4.8.5/2.7.1
#module load python/gcc-4.8.5/3.12.1
#module load openmpi/intel-2021.4.0/4.1.0

module load intel/2021.4.0
module load hdf5/gcc-4.8.5/1.12.0

export CFLAGS="-std=c99"
export LD_LIBRARY_PATH="/home/taeung/.local/lib:/home/taeung/GSL/lib:$LD_LIBRARY_PATH"
export PATH="/home/taeung/.local/bin:$PATH"
export PS1='[$PWD:] \[\033[1;37m\]'

alias qstf="qstat -xml -u '*' | tr '\n' ' ' | sed 's#<job_list[^>]*>#\n#g'   | sed 's#<[^>]*>##g' | grep ' ' | column -t"
alias qsta="qstat -u '*'"
alias ls="ls --color"
alias py="python3"
