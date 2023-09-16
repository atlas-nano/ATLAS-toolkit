proc diffusion2 { sel f1 } {
    set fout [open $f1 w]
    set nf [molinfo [$sel molid] get numframes]

    for { set i 0 } { $i < $nf } { incr i } {
        $sel frame $i
        lappend cxy [$sel get {x y}]
    }

    set ds 0.0
    set lds 0.0
    puts $fout "0 0.0"

    for { set i 1 } { $i < [expr $nf ] } { incr i } {
        for { set j 0 } { $j < [expr $nf - $i] } { incr j } {
            set k [expr $j + $i]
            set ds [expr $ds + [veclength [vecsub [lindex [lindex $cxy $k] 0] \
                [lindex [lindex $cxy $j] 0]]]]
        }
        puts $fout "$i [expr $ds/($nf-$i)]"
        set ds 0.0
    }
    close $fout
} 
