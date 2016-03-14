#!/usr/bin/perl -w
# Extract velocity dispersion profiles from star1.list file

###DOESN'T WORK !!!  BECAUSE STRUCTURE OF list FILE DEPENDS ON NUMBER OF COMPONENTS!

#use strict;

use Getopt::Long;

my $time;

$i_step = 0;
$base_name = "sigma_profile";

while (<>) {
  chomp;
  if (/ +GAMMA= *([^ ]+) *XBIN2B=/) {
    $gamma=$1;
  }
  elsif (/^ *XNTOTA= *([^ ]+) *XIMF=/) {
    $Nstar_tot=$1;
  }
  elsif (/^ *T= *(.+) DT=/) {
    $time=$1;
    #print "T = ",$time,"\n";
  }
  elsif (/^ *RADIUS /) { # header line announcing data block

    # find column numbers for various components
    # column zero is the number of the shell, column 1 is the radius
    s/([A-Z]) +([0-9])/$1_$2/g;
    @headers = split;

    @columns_num_rho  = (); # densites
    @columns_num_sigr = (); # radial velocity dispersions
    @columns_num_sigt = (); # tangential velocity dispersions

    for  ($i_col=0; $i_col<=$#headers; ++$i_col) {
      if ( $headers[$i_col] =~ /^RHO_([1-9][0-9]*)$/ ) {
	$columns_num_rho[$1-1] = $i_col+1;
	#print STDOUT $headers[$i_col]," found in column ",$i_col+1,"\n";
      }
      elsif ( $headers[$i_col] =~ /^SIGR_([1-9][0-9]*)$/ ) {
	$columns_num_sigr[$1-1] = $i_col+1;
	#print STDOUT $headers[$i_col]," found in column ",$i_col+1,"\n";
      }
      elsif ( $headers[$i_col] =~ /^SIGT_([1-9][0-9]*)$/ ) {
	$columns_num_sigt[$1-1] = $i_col+1;
	#print STDOUT $headers[$i_col]," found in column ",$i_col+1,"\n";
      }
    }

    # read block of data
    # ==================

    @rho_data=()  if (@columns_num_rho );
    @sigt_data=() if (@columns_num_sigt);
    @sigr_data=() if (@columns_num_sigr);
    @radius=();

    while (<>,  /^ *$/) {} # first, skip possible blank lines
    $i_shell = 0;
  DATA_READ: while (<>) {
      s/\*/ /g;
      chomp;
      @data_line = split;
      last DATA_READ if ($#data_line != $#headers+1);
      $radius[$i_shell]    = $data_line[1];
      $rho_data[$i_shell]  = [ @data_line[@columns_num_rho ] ] if (@columns_num_rho ); # densites
      $sigt_data[$i_shell] = [ @data_line[@columns_num_sigt] ] if (@columns_num_sigt); # radial velocity dispersions
      $sigr_data[$i_shell] = [ @data_line[@columns_num_sigr] ] if (@columns_num_sigr); # tangential velocity dispersions
      $i_shell++;
    }
    $i_shell--;

    # Do we have complete density and velocity data?
    if (@rho_data && @sigt_data && @sigr_data) {
      $nb_comp = ($#columns_num_rho>$#columns_num_sigt  ) ? $#columns_num_rho+1 : $#columns_num_sigt+1;
      $nb_comp = ($nb_comp         >$#columns_num_sigr+1) ? $nb_comp            : $#columns_num_sigr+1;
      ($nb_comp>1) and print STDOUT "Found complete information for ",$nb_comp," components\n"
	           or  print STDOUT "Found complete information for 1 component\n";
      
      # Compute average velocity dispersion
      @sigr_avrg =();
      @sigt_avrg =();
      @tot_dens  =();
      for ($j_shell = 0; $j_shell<=$i_shell; $j_shell++) {
	$sigr_avrg[$j_shell] = 0;
	$sigt_avrg[$j_shell] = 0;
	$tot_dens[$j_shell]  = 0;
	for  ($i_comp=0; $i_comp<=$nb_comp-1; ++$i_comp) {
	  $tot_dens[$j_shell]=$tot_dens[$j_shell]+$rho_data[$j_shell][$i_comp];
	  $sigr_avrg[$j_shell] = $sigr_avrg[$j_shell] + $rho_data[$j_shell][$i_comp]*$sigr_data[$j_shell][$i_comp]**2;
	  $sigt_avrg[$j_shell] = $sigt_avrg[$j_shell] + $rho_data[$j_shell][$i_comp]*$sigt_data[$j_shell][$i_comp]**2;
	}
	$sigr_avrg[$j_shell] = sqrt($sigr_avrg[$j_shell]/$tot_dens[$j_shell]);
	$sigt_avrg[$j_shell] = sqrt($sigt_avrg[$j_shell]/$tot_dens[$j_shell]);
      }
      
      # write data to a file
      # ====================
      
      $file_name = $base_name . $i_step . ".asc";
      print STDOUT "T = ",$time," file : ",$file_name,"\n";
      open(OUTFILE,">$file_name");
      print OUTFILE "# Velocity dispersion profile from SPEDI simulation\n";
      print OUTFILE "# N-body units are used\n";
      print OUTFILE "# Nstar = ",$Nstar_tot,"\n";
      print OUTFILE "# gamma_Coulomb = ",$gamma,"\n";
      print OUTFILE "# Time = $time\n";
      print OUTFILE "# Number of mass components = ",$nb_comp,"\n";
      print OUTFILE "# 1: Radius 2: Rho_tot 3: SigAvrgR 4: SigAvrgT";
      if ($nb_comp>1) {
	for  ($i_comp=1; $i_comp<=$nb_comp; ++$i_comp) {
	  print OUTFILE $i_comp+4,": SigR$i_comp ";
	}
	for  ($i_comp=1; $i_comp<=$nb_comp; ++$i_comp) {
	  print OUTFILE $i_comp+$nb_comp+4,": SigT$i_comp ";
	}
      }
      print OUTFILE "\n";
      for ($j_shell = 0; $j_shell<=$i_shell; $j_shell++) {
	print OUTFILE $radius[$j_shell]," ",sprintf('%12.5e 'x3,$tot_dens[$j_shell],$sigr_avrg[$j_shell],$sigt_avrg[$j_shell]);
	if ($nb_comp>1) {
	  $prov=$sigr_data[$j_shell];
	  print OUTFILE " ",join(" ",@$prov);
	  $prov=$sigt_data[$j_shell];
	  print OUTFILE " ",join(" ",@$prov);
	}
	print OUTFILE "\n";
      }
      close(OUTFILE);
      print "\n";
      ++$i_step;
      
      @rho_data=();
      @sigt_data=();
      @sigr_data=();
      @radius=();
    }
  }
}
