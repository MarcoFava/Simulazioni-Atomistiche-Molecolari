#!/bin/bash

remote=s327286@bora.units.it:
local=./bora
folder_r=.
folder_l=project
# excludes=(--exclude env/ --exclude .f2py-jit/ --exclude __pycache__/)
excludes=(--exclude env/ --exclude .pantarei/ --exclude .f2py-jit/ --exclude *__pycache__/)

case $1 in
	up)
	rsync -uva "${excludes[@]}" "$local/$folder_l" "$remote../$folder_r"
	;;

	down)
	rsync -uva "${excludes[@]}" "$remote$folder_r" "$local/$folder_l"
	;;
esac
