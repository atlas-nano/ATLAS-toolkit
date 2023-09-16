#!/usr/bin/env bash
# Title         : GetStates.sh
# Description   : Retrieve highest contributing states to computed XAS spectrum
# Author        : Sebastiaan P. Huber
# Date          : 2015-03-18
# Version       : 0.0.1    
# Usage         : bash getstates.sh
# Notes         : 

# Set the version of the script
version="0.0.1"

# Script usage and help documentation
usage()
{
	echo -n "
 $(basename $0) [OPTIONS]... DIRECTORY

 DIRECTORY              Specify the directory in which the script will search for the \"*stick.0\" file
                        If not specified the current working directory will be used

 Options:
  -a, --atom            Filter for states by atom number e.g. B288 will only show states for excited boron atom #288
  -m, --min-energy      Set the minimum energy for which bands are printed
  -M, --max-energy      Set the maximum energy for which bands are printed
  -b, --bandnumber      Filter all states with exactly this band number. Overrides all other energy and bandnumber filters
  -n, --max-bandnumber  Set the maximum band number that is to be printed. Script tries to get it from nscf.in or defaults to 100
  -I, --intensity       Filter out all states with a intensity below this value, relative to largest oscillator strength
  -f, --format          Format output as space delimited string
  -c, --count           Limit the number of states that are printed per atom
  -s, --spin            Select only sticks from certain spin channel 1: spin up   2: spin down   0: both (default)
  -q, --quiet           Suppress all output
  -v, --verbose         Set ouput verbosity to high
      --debug           Set debug mode
      --help            Display this help and exit
      --version         Output version information and exit
      
"
}

# Default values
opt_atom=
opt_spin=0
opt_format=0
opt_intensity=0.4
opt_count=10
opt_bandnumber=
opt_max_bandnumber=0
OPT_MAX_BANDNUMBER=500
opt_min_energy=-1E-300
opt_max_energy=1E300

# Read the options and set the corresponding variables
while [[ $1 = -?* ]]; do
	case $1 in
		-a|--atom)           shift; opt_atom=$1 ;;
		-m|--min-energy)     shift; opt_min_energy=$1 ;;
		-M|--max-energy)     shift; opt_max_energy=$1 ;;
		-b|--bandnumber)     shift; opt_bandnumber=$1 ;;
		-n|--max-bandnumber) shift; opt_max_bandnumber=$1 ;;
		-I|--intensity)      shift; opt_intensity=$1 ;;
		-f|--format)         opt_format=1 ;;
		-c|--count)          shift; opt_count=$1 ;;
		-s|--spin)           shift; opt_spin=$1 ;;
		-q|--quiet)          LOG_LEVEL=0 ;;
		-v|--verbose)        LOG_LEVEL=6 ;;
		--debug)             LOG_LEVEL=7 ;;
		--help)              usage >&2; exit 0 ;;
		--version)           echo "$(basename $0) $version"; exit 0 ;;
		--endopts)           shift; break ;;
		*)                   print_die "Invalid option: $1\nSee '$(basename $0) --help' for more information" ;;
	esac
	shift
done
args="$@"

# Main body definition
main()
{	
	if [ -z $1 ]; then
		dir=$(pwd)
	else
		dir="$1"
	fi

	for directory in $(find "$dir" -name '*.stick.*' -printf '%h\n' | sort -u); do

		basename=$(basename $directory)
		if [ -n "$opt_atom" ] && [ "$basename" != "$opt_atom" ]; then
			continue
		fi

		if [ $opt_max_bandnumber -eq 0 ]; then
			nscf_file=$(find "$directory" -type f -name "*nscf.in")

			if [ ! -r "$nscf_file" ]; then
				printf "%s\n" "Error: could not find the nscf.in file setting band number cutoff to default: $OPT_MAX_BANDNUMBER"
				max_bandnumber=$OPT_MAX_BANDNUMBER
			else
				max_bandnumber=$(sed -nr 's/^.*nbnd=([0-9]*).*$/\1/p' "$nscf_file")
			fi
		else
			max_bandnumber=$opt_max_bandnumber
		fi


		nscf_out=$(find "$directory" -type f -name "*nscf.out")
		if [ ! -r "$nscf_out" ]; then
			printf "%s\n" "Error: could not find the nscf.out file for this atom so can not determine number of electrons... skipping"
			continue
		fi
		occupied_bands=$(grep "number of electrons" "$nscf_out" | awk '{printf "%d\n", $5/2}')

		atom=${directory##*/}
		file=$(find "$dir" -type f -name "*.stick.0" -path "*/$atom/*")

		[ -f "$file" ] || continue;

		# This version assumes the following column definition of stick files produced by the code
		# 1: Band number
		# 2: K-point
		# 3: Spin
		# 4: Band energy
		# 5: Oscillator strength
		if [ -z $opt_bandnumber ]; then
			states=$(awk -v min_bandnumber=$occupied_bands -v max_bandnumber=$max_bandnumber -v m=$opt_min_energy -v M=$opt_max_energy -v spin=$opt_spin '{
					if (spin == 0) {
						if ($2 == 1 && min_bandnumber < $1 && $1 <= max_bandnumber && m < $4 && $4 < M) {printf "%d %d %19.12E %19.12E\n", $1, $3, $4, $5}
					} else {
						if ($2 == 1 && min_bandnumber < $1 && $1 <= max_bandnumber && m < $4 && $4 < M && $3 == spin) {printf "%d %d %19.12E %19.12E\n", $1, $3, $4, $5}
					}
				}' "$file" | sort -gr -k4 | head -n $opt_count || true)
		else
			states=$(awk -v bandnumber=$opt_bandnumber -v spin=$opt_spin '{
					if (spin == 0) {
						if ($2 == 1 && $1 == bandnumber) {printf "%d %d %19.12E %19.12E\n", $1, $3, $4, $5}
					} else {
						if ($2 == 1 && $1 == bandnumber && $3 == spin) {printf "%d %d %19.12E %19.12E\n", $1, $3, $4, $5}
					}
				}' "$file" | sort -gr -k4 | head -n $opt_count || true)
		fi
		ostrength=$(echo "$states" | awk '(NR==1){print $4}')
		
		for state in "$states"; do
		if [ $opt_format -eq 1 ]; then
			# Produce output in format that is readable by XAS_state.sh template script
			# format: (atom#bandnumber#spin, ...)
			echo "$state" | \
			awk -v atom=$atom -v I=$opt_intensity -v os=$ostrength '{
				if ($4/os>I) {printf "%s#%d#%d ", atom, $1, $2}
			}'
			echo
		else
			printf "%-6s%-6s%-6s%-12s%-20s\n" "Atom" "Spin" "Band" "Energy" "Relative strength"
			echo "$state" | \
			awk -v atom=$atom -v I=$opt_intensity -v os=$ostrength '{
				if ($4/os>I) {printf "%-4s  %4d  %4d  %10f  %10f\n", atom, $2, $1, $3, $4/os}
			}'
			echo
		fi
		done
	done

	exit 0;
}

# Execute the main body of the script
main $args

exit 0;
