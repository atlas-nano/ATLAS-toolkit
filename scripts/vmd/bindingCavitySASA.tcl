# Prints the SASA of a cavity in a molecule each timestep
proc print_cavity_sasa {mol sys cav} {
	# use frame 0 for the reference
	#set protein [atomselect $mol $sys frame 0]
	set protein [atomselect top "protein" frame 0]
	# the cavity
	#set cavity [atomselect $mol $cav]
	set cavity [atomselect top "protein and resid 60 to 63 67 74 139 215"]

	set num_steps [molinfo $mol get numframes]
	set fw [open "cav_sasa.dat" w]
	for {set frame 0} {$frame < $num_steps} {incr frame} {
		# get the correct frame
		$protein frame $frame
		$protein update
		$cavity frame $frame
		$cavity update
		set cav_sasa [ measure sasa 1.4 $protein -restrict $cavity ]
		# print the sasa
		puts "SASA of $frame is $cav_sasa"
		puts $fw "$frame $cav_sasa"
	}
	close $fw
}
