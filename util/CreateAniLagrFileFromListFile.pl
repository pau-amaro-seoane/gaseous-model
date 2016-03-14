#!/usr/bin/perl -w

# Create a file depicting the evolution of the anisotropy at selected Lagrang radii from star1.list file


#use strict;

use Getopt::Long;

@frac_lag = @ARGV;

(@frac_lag) or @frac_lag = ( 0.05, 0.1, 0.2, 0.5, 0.75, 0.95 );

my $time;

$i_step = 0;

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
    $totmass_found=0;

    for  ($i_col=0; $i_col<=$#headers; ++$i_col) {
      if ( $headers[$i_col] =~ /^TOTMASS$/ ) {
	$column_totmass = $i_col+1;
	$totmass_found=1;
	print STDOUT $headers[$i_col]," found in column ",$i_col+1,"\n";
      }
      if ( $headers[$i_col] =~ /^RHO_([1-9][0-9]*)$/ ) {
	$columns_num_rho[$1-1] = $i_col+1;
	print STDOUT $headers[$i_col]," found in column ",$i_col+1,"\n";
      }
      elsif ( $headers[$i_col] =~ /^SIGR_([1-9][0-9]*)$/ ) {
	$columns_num_sigr[$1-1] = $i_col+1;
	print STDOUT $headers[$i_col]," found in column ",$i_col+1,"\n";
      }
      elsif ( $headers[$i_col] =~ /^SIGT_([1-9][0-9]*)$/ ) {
	$columns_num_sigt[$1-1] = $i_col+1;
	print STDOUT $headers[$i_col]," found in column ",$i_col+1,"\n";
      }
    }

    # read block of data
    # ==================

    @rho_data=()  if (@columns_num_rho );
    @sigt_data=() if (@columns_num_sigt);
    @sigr_data=() if (@columns_num_sigr);
    @mass_data=() if ($totmass_found);
    @radius=();

    while (<>,  /^ *$/) {} # first, skip possible blank lines
    $i_shell = 0;
  DATA_READ: while (<>) {
      s/\*/ /g;
      chomp;
      @data_line = split;
      last DATA_READ if ($#data_line != $#headers+1);
      $radius[$i_shell]    = $data_line[1];
      $mass_data[$i_shell] = $data_line[$column_totmass]       if ($totmass_found);
      $rho_data[$i_shell]  = [ @data_line[@columns_num_rho ] ] if (@columns_num_rho ); # densites
      $sigt_data[$i_shell] = [ @data_line[@columns_num_sigt] ] if (@columns_num_sigt); # radial velocity dispersions
      $sigr_data[$i_shell] = [ @data_line[@columns_num_sigr] ] if (@columns_num_sigr); # tangential velocity dispersions
      $i_shell++;
    }
    $i_shell--;

    # Do we have complete density and velocity data?
    if ($totmass_found && @rho_data && @sigt_data && @sigr_data) {
      @prov = @{$rho_data[0]};
      $nb_comp = $#prov+1;
      ($nb_comp>1) and print STDOUT "Found complete information for ",$nb_comp," components\n"
	           or  print STDOUT "Found complete information for 1 component\n";

      # Compute average velocity dispersion
      @sigr_avrg =();
      @sigt_avrg =();
      @anis_avrg =();
      for ($j_shell = 0; $j_shell<=$i_shell; $j_shell++) {
	$sigr_avrg[$j_shell] = 0;
	$sigt_avrg[$j_shell] = 0;
	$tot_dens=0;
	for  ($i_comp=0; $i_comp<=$nb_comp-1; ++$i_comp) {
	  $tot_dens=$tot_dens+$rho_data[$j_shell][$i_comp];
	  $sigr_avrg[$j_shell] = $sigr_avrg[$j_shell] + $rho_data[$j_shell][$i_comp]*$sigr_data[$j_shell][$i_comp]**2;
	  $sigt_avrg[$j_shell] = $sigt_avrg[$j_shell] + $rho_data[$j_shell][$i_comp]*$sigt_data[$j_shell][$i_comp]**2;
	}
	$sigr_avrg[$j_shell] = sqrt($sigr_avrg[$j_shell]/$tot_dens);
	$sigt_avrg[$j_shell] = sqrt($sigt_avrg[$j_shell]/$tot_dens);
	$anis_avrg[$j_shell] = 2 - 2*($sigt_avrg[$j_shell]/$sigr_avrg[$j_shell])**2;
      }
      
      # Find specified Lagrangian radii for current mass profile
      # and compute average stellar mass
      # ========================================================

      # Get total mass
      $tot_mass = $mass_data[$#mass_data];
      $limit_mass = 0.0;
      $i_shell = 0;
      @ani_lag=0*@frac_lag;
      for ($i_lag=0; $i_lag<=$#frac_lag; $i_lag++) {
	$limit_mass = $frac_lag[$i_lag]*$tot_mass;
	while ( ($mass_data[$i_shell] < $limit_mass) && ($mass_data[$i_shell] < $tot_mass) ) {
	  $i_shell++;
	}
	# get masses in each component
	if ($i_shell==0) {
	  $ani_lag[$i_lag]=$anis_avrg[0];
	} elsif ($i_shell==$#mass_data) {
	  $ani_lag[$i_lag]=$anis_avrg[$#mass_data];
	} else { # linear interpolation
	  $M_inf=$mass_data[$i_shell-1];
	  $M_sup=$mass_data[$i_shell  ];
	  $a_inf=$anis_avrg[$i_shell-1];
	  $a_sup=$anis_avrg[$i_shell  ];
	  if ($M_sup>$M_inf) { # Linear interpolation
	    print STDERR "a_inf= ",$a_inf," a_sup=",$a_sup," M_inf=",$M_inf," M_sup=",$M_sup," limit_mass=",$limit_mass,"\n";
	    $ani_lag[$i_lag] = $a_inf + ($a_sup-$a_inf)/(log($M_sup)-log($M_inf))*(log($limit_mass)-log($M_inf));
	  } else {
	    $ani_lag[$i_lag] = $a_inf;
	  }
	}
      }
      
      # Write Lagr. radii to a file
      # ===========================
      if ($file_started) { # Append to file
	print STDERR "> Adding data to $file_name (T=$time)...\n";
      } else { # Create file and write header
	$file_name = 'lagr_anis.asc';
	print STDERR "> Creating $file_name...\n";
	open(OUTFILE,">$file_name");
	print OUTFILE "# Lagrangian anisotropies from SPEDI simulation\n";
	print OUTFILE "# Anisotropy parameter A=2-<Vtg^2>/<Vrad^2>\n";
	print OUTFILE "# Nstar = ",$Nstar_tot,"\n";
	print OUTFILE "# gamma_Coulomb = ",$gamma,"\n";
	print OUTFILE "# 1: Time ";
	for  ($i_lag=0; $i_lag<=$#frac_lag; $i_lag++) {
	  print OUTFILE $i_lag+2,": Ani",sprintf("%6.6d ",1e5*$frac_lag[$i_lag]);
	}
	print OUTFILE "\n";
	$file_started = 1;
      }
      print OUTFILE (sprintf("%12.5e " x ($#ani_lag + 2),$time,@ani_lag),"\n");
      ++$i_step;
      
      @rho_data=();
      @sigt_data=();
      @sigr_data=();
      @radius=();
      $totmass_found=0;
    }
  }
}
