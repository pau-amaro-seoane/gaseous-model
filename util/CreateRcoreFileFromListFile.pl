#!/usr/bin/perl -w
# Extract velocity dispersion profiles from star1.list file

###DOESN'T WORK !!!  BECAUSE STRUCTURE OF list FILE DEPENDS ON NUMBER OF COMPONENTS!

#use strict;

use Getopt::Long;

my $time;
my $pi=3.141592653589793116;

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

    while (<>,  /^ *$/) {} # first, skip possible blank lines
    $_ = <>;
    s/\*/ /g;
    chomp;
    @data_line = split;
    @rho_data  = @data_line[@columns_num_rho ] if (@columns_num_rho ); # central densites
    @sigt_data = @data_line[@columns_num_sigt] if (@columns_num_sigt); # central radial velocity dispersions
    @sigr_data = @data_line[@columns_num_sigr] if (@columns_num_sigr); # central tangential velocity dispersions
    
    # Do we have complete density and velocity data?
    if (@rho_data && @sigt_data && @sigr_data) {
      $nb_comp = $#rho_data+1;
      ($nb_comp>1) and print STDOUT "Found complete information for ",$nb_comp," components\n"
	           or  print STDOUT "Found complete information for 1 component\n";
      # Compute average velocity dispersion
      $sigv2_avrg = 0;
      $tot_dens   = 0;
      for ($i_comp=0; $i_comp<=$nb_comp-1; ++$i_comp) {
	$tot_dens   = $tot_dens   + $rho_data[$i_comp];
	$sigv2_avrg = $sigv2_avrg + $rho_data[$i_comp]*($sigr_data[$i_comp]**2+2*$sigt_data[$i_comp]**2);
      }
      $sigv2_avrg = $sigv2_avrg/$tot_dens;
      $Rcore_avrg = sqrt(3*$sigv2_avrg/(4*$pi*$tot_dens));
      @Rcore = 0*@rho_data;
      for ($i_comp=0; $i_comp<=$nb_comp-1; ++$i_comp) {
	$Rcore[$i_comp] = sqrt(3*($sigr_data[$i_comp]**2+2*$sigt_data[$i_comp]**2)/(4*$pi*$tot_dens));
      }

      # Write Lagr. radii to a file
      # ===========================
      $file_name = 'Rcore.asc';
      if ($header_written) { # Append to file
	print STDOUT "> Adding data to $file_name (T=$time)...\n";
      } else { # Create file and write header
	print STDOUT "> Creating $file_name...\n";
	open(OUTFILE,">$file_name");
	print OUTFILE "# Core radii from SPEDI simulation\n";
	print OUTFILE "# N-body units are used\n";
	print OUTFILE "# Nstar = ",$Nstar,"\n";
	print OUTFILE "# gamma_Coulomb = ",$gamma,"\n";
	print OUTFILE "# Number of mass components = ",$nb_comp,"\n";
       	print OUTFILE "# 1: Time 2: Rcore_avrg ";
	if ($nb_comp>1) {
	  for  ($i_comp=1; $i_comp<=$nb_comp; ++$i_comp) {
	    print OUTFILE $i_comp+2,": Rcore$i_comp ";
	  }
	}
	print OUTFILE "\n";
	$header_written = 1;
      }
      print OUTFILE sprintf("%16.9e",$time)," ",sprintf("%12.5e",$Rcore_avrg);
      if ($nb_comp>1) {
	print OUTFILE sprintf(" %12.5e"x$nb_comp,@Rcore);
      }
      print OUTFILE "\n";
      @rho_data=();
      @sigt_data=();
      @sigr_data=();
    }
  }
}
close OUTFILE;

