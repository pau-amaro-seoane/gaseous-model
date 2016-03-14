#!/usr/bin/perl -w

#use strict;

use Getopt::Long;

@frac_lag = @ARGV;
@ARGV=();

my $time;

$i_step = 0;
$base_name = "lagr_radii_";

while (<>) {
  chomp;
  # transform all "double precision" numbers with "D" notation (from fortran code output) to standard "E" notation
  #s/([+-]?(?:[0-9]+|\.[0-9]+|[0-9]+\.[0-9]*))(?:D|d)([+-]?[0-9]+)/$1E$2/g;

  if (/ +GAMMA= *([^ ]+) *XBIN2B=/) {
    $gamma=$1;
  }
  elsif (/^ *XNTOTA= *([^ ]+) *XIMF=/) {
    $Nstar=$1;
  }
  elsif (/^ *T= *(\S+) DT= *\S+ *RHO= *(\S+)/) {
    $time=$1;
    $rhoc=$2;
    #print "T = ",$time," Rho_c = ",$rhoc,"\n";
  }
  elsif (/^ *RADIUS +TOTMASS .* MASS/) { # header line announcing data about Mr(R) for all components

    # find column numbers for various components
    # column zero is the number of the shell, column 1 is the radius
    s/([A-Z]) +([0-9])/$1_$2/g;
    @headers = split;
    @columns_num = ();
    for  ($i_col=0; $i_col<=$#headers; ++$i_col) {
      if ( $headers[$i_col] =~ /^TOTMASS$/ ) {
	$columns_num[0] = $i_col+1;
      } elsif ( $headers[$i_col] =~ /^MASS_([1-9][0-9]*)$/ ) {
	$columns_num[$1] = $i_col+1;
      }
    }
    $nb_comp = $#columns_num;
    ($nb_comp == 1) and @columns_num = $columns_num[0];

    # read block of data
    # ==================

    @mass_data=();
    @radius=();

    while (<>,  /^ *$/) {} # first, skip possible blank lines
    $i_shell = 0;
  DATA_READ: while (<>) {
      s/\*/ /g;
      chomp;
      @data_line = split;
      last DATA_READ if ($#data_line < $#headers);
      $radius[$i_shell] = $data_line[1];
      $mass_data[$i_shell] = [ @data_line[@columns_num] ]; # array of references to arrays
      $i_shell++;
    }
    $i_shell--;
    $n_shells=$i_shell;

    # Find Lagrangian radii for current mass profiles
    # ===============================================
    $nmax = $nb_comp;
    $nmax = 0 if ($nb_comp==1);
    for ($i_comp=0; $i_comp<=$nmax; ++$i_comp) {

      # Get total mass for the component
      $tot_comp_mass = $mass_data[$n_shells][$i_comp];
      $limit_mass = 0.0;
      $i_shell = 0;
      for ($i_lag=0; $i_lag<=$#frac_lag; $i_lag++) {
	$limit_mass = $frac_lag[$i_lag]*$tot_comp_mass;
	while ( ($mass_data[$i_shell][$i_comp] < $limit_mass) && ($mass_data[$i_shell][$i_comp] < $tot_comp_mass) ) {
	  $i_shell++;
	}
	if ($i_shell==0) {
	  $r_lag[$i_lag] = $radius[$i_shell];
	} elsif ($i_shell==$n_shells) {
	  $r_lag[$i_lag] = $radius[$n_shells];
	} else { # linear interpolation
	  $r_inf=$radius[$i_shell-1];
	  $r_sup=$radius[$i_shell];
	  $m_inf=$mass_data[$i_shell-1][$i_comp];
	  $m_sup=$mass_data[$i_shell][$i_comp];
	  if ($m_sup>$m_inf) { # Linear interpolation
	    # Interpolation in the log gives smoother Lagr. radii evolution
	    $r_lag[$i_lag] = exp(log($r_inf) + (log($r_sup)-log($r_inf))/(log($m_sup)-log($m_inf))*(log($limit_mass)-log($m_inf)));
	  } else {
	    $r_lag[$i_lag] = $r_inf;
	  }
	}
      }
      # Write Lagr. radii to a file
      # ===========================
      $file_name = $base_name . $i_comp . '.asc';
      if ($header_written[$i_comp]) { # Append to file
	print STDOUT "> Adding data to $file_name (T=$time)...\n";
	open(OUTFILE,">>$file_name");
      } else { # Create file and write header
	print STDOUT "> Creating $file_name...\n";
	open(OUTFILE,">$file_name");
	print OUTFILE "# Lagrangian radii from SPEDI simulation\n";
	print OUTFILE "# N-body units are used\n";
	print OUTFILE "# Nstar = ",$Nstar,"\n";
	print OUTFILE "# gamma_Coulomb = ",$gamma,"\n";
       	print OUTFILE "# Mass component : $i_comp (0 means all components)\n";
       	print OUTFILE "# 1: Time 2: Rho_c 3: TotalComponentMass ";
	for  ($i_lag=0; $i_lag<=$#frac_lag; $i_lag++) {
	  print OUTFILE $i_lag+4,": R",sprintf("%7.7d ",1e6*$frac_lag[$i_lag]);
	}
	print OUTFILE "\n";
	$header_written[$i_comp] = 1;
      }
      print OUTFILE (sprintf("%16.9e" . ("%12.5e " x ($#r_lag + 3)),$time,$rhoc,$tot_comp_mass,@r_lag),"\n");
      close(OUTFILE);
    }

    ++$i_step;
  }
}
