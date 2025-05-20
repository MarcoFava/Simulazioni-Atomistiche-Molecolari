# this is the makefile

initialize:
	echo '&nml density=0.1, nc=3 /' | ./initialize.x

run: # stamp last column of the row starting with 'Run averages'
	echo '&nml nblock=10, nstep=1000 /' | ./mc_nvt_hs.x | grep 'Run averages' | awk '{print $NF}'

all: initialize run
