# simple shell script to be used with SGE for launching glide
# tell SGE to change into current working directory
#$-cwd

export OMP_NUM_THREADS=$NSLOTS

echo Launching $1 on $HOSTNAME using $NSLOTS slots

$*
