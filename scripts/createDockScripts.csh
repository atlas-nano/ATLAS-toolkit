#!/bin/tcsh
#!/bin/csh

if ($#argv < 3) then
echo "usage: $0 [dealanized protein] [ligand] [spheres] [save prefix (optional)]"
exit(1)
endif

set ligand = `/ul/tpascal/scripts/getFileAbsPath.pl $argv[2]`
set ligFullname = `basename $ligand .bgf`
set spheres = `/ul/tpascal/scripts/getFileAbsPath.pl $argv[3]`
set fullpath = `/ul/tpascal/scripts/getFileAbsPath.pl $argv[1]`
set molname = `basename $fullpath .bgf`
set protein_dir = `dirname $fullpath`
set protein_dealanized = $fullpath
set protein_alanized = ${protein_dir}/${molname}.ala.sph.bgf
set ligdir = `dirname $ligand`
set ligname = `basename $ligFullname _dreidii`
set molname = `basename $molname _dreidii`

if ((-e $1) == 0) then
echo "Error: Cannot locate protein file $1!\n"
exit(1)
else if ((-e $protein_alanized) == 0) then
echo "Error: Cannot locate alanized protein file $protein_alanized"
exit(1)
else if ((-e $ligand) == 0) then
echo "Error: Cannot locate ligand file $ligand\n"
exit(1)
else if ((-e $spheres) == 0) then
echo "Error: Cannot locate spheres file $spheres"
exit(1)
else if ($#argv > 3) then
set molname = $argv[4]
endif

#Create input files
cat >> ${molname}_${ligname}_GenDock.in <<DATA;
@DarwinDock
	protein_full_path 	$protein_alanized
	ligand_full_path	${ligdir}/${ligFullname}.bgf
	spheres_full_path 	$spheres
	diversity		1.2
	step_size		5000
	fraction_families	0.10
	final_poses		100
	complete_threshold	0.05
	bump_cutoff		6
@ScreamUnifiedBindingSite
	dealanized_protein_path $protein_dealanized
	scoring_energy 		unified_cavity
	reference_energies 	interaction
	reference_energies 	total
	reference_energies 	local_cavity
	reference_energies 	full_delphi
	reference_energies 	partial_delphi
@BindSiteMinimize
	bindsite_minimization_steps	50
        bindsite_minimization_cmm_level	0
	unified_bind_site_radius	4.0
	keep_percent_bindsiteanneal	50
	scoring_energy 		total
	reference_energies 	interaction
	reference_energies 	full_delphi
	reference_energies 	unified_cavity
	reference_energies 	local_cavity
	reference_energies 	partial_delphi
@ComplexMinimize
	scoring_energy 		total
	reference_energies 	interaction
	reference_energies 	full_delphi
	reference_energies 	unified_cavity
	reference_energies 	local_cavity
	reference_energies 	partial_delphi

DATA

cat >> ${molname}_${ligname}_neutral_GenDock.in <<DATA;
@DarwinDock
        protein_full_path       $protein_alanized
        ligand_full_path        ${ligdir}/${ligFullname}.bgf
        spheres_full_path       $spheres
        diversity               1.2
        step_size               5000
        fraction_families       0.10
        final_poses             100
        complete_threshold      0.05
        bump_cutoff             6
@ScreamUnifiedBindingSite
        dealanized_protein_path $protein_dealanized
        scoring_energy          unified_cavity
        reference_energies      interaction
        reference_energies      total
        reference_energies      local_cavity
        reference_energies      full_delphi
        reference_energies      partial_delphi
@Neutralize
        neutral_definitions     /ul/victor/SCREAM/SCREAM_v2.0/SCREAM/src/Neutralize/Standard.def
        neutral_ligand          $ligand
        scoring_energy          full_delphi
        reference_energies      interaction
        reference_energies      total
        reference_energies      unified_cavity
        reference_energies      local_cavity
        reference_energies      partial_delphi
@BindSiteMinimize
        bindsite_minimization_steps     50
        bindsite_minimization_cmm_level 0
        unified_bind_site_radius        4.0
        keep_percent_bindsiteanneal     50
        scoring_energy          total
        reference_energies      interaction
        reference_energies      full_delphi
        reference_energies      unified_cavity
        reference_energies      local_cavity
        reference_energies      partial_delphi
@ComplexMinimize
        scoring_energy          total
        reference_energies      interaction
        reference_energies      full_delphi
        reference_energies      unified_cavity
        reference_energies      local_cavity
        reference_energies      partial_delphi

DATA

#create script files
cat >> ${molname}_${ligname}_GenDock.script <<DATA;
#PBS -l nodes=1:ppn=1
#PBS -l walltime=960:00:00
#PBS -q workq
#PBS -j oe
#PBS -N ${molname}_${ligname}_GenDock
#!/bin/csh

cd $PWD
/project/Biogroup/Software/Docking/GenDock.pl -i ${molname}_${ligname}_GenDock.in -n ${molname}_${ligname}_GenDock_run

DATA

cat >> ${molname}_${ligname}_neutral_GenDock.script <<DATA;
#PBS -l nodes=1,walltime=960:00:00
#PBS -q workq
#PBS -j oe
#PBS -N ${molname}_${ligname}_neutral_GenDock
#!/bin/csh

cd $PWD
/project/Biogroup/Software/Docking/GenDock.pl -i ${molname}_${ligname}_neutral_GenDock.in -n ${molname}_${ligname}_neutral_GenDock_run

DATA


