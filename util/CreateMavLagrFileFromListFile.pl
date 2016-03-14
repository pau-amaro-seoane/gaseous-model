#!/usr/bin/perl -w

#use strict;

use Getopt::Long;

@frac_lag = @ARGV;
@ARGV=();

my $time;

$i_step = 0;
@N_comp=();
@M_comp=();
@m_comp=();

while (<>) {
  chomp;

  # Get global properties of components
  if (/^  *I= *(\S+) +NI= *(\S*) +Mi= *(\S+) +mi= *(\S+) *$/) {
    $icomp=$1;
    $N_comp[$icomp]=$2;
    $M_comp[$icomp]=$3;
    $m_comp[$icomp]=$4;
    print STDERR "> component $icomp: N=$N_comp[$icomp], M=$M_comp[$icomp] Msun, m=$m_comp[$icomp] Msun\n";
  }
  elsif (/ +GAMMA= *([^ ]+) *XBIN2B=/) {
    $gamma=$1;
  }
  elsif (/^ *XNTOTA= *([^ ]+) *XIMF=/) {
    $Nstar_tot=$1;
  }
  elsif (/^ *T= *(.+) DT=/) {
    $time=$1;
  }
  elsif (/^ *RADIUS +TOTMASS .* MASS/) { # header line announcing data about Mr(R) for all components

    # Find column numbers for various components
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
    print STDERR "> Data for T=$time read\n";
    $i_shell--;
    $n_shells=$i_shell;

    # Find specified Lagrangian radii for current mass profile
    # and compute average stellar mass
    # ========================================================
    $nmax = $nb_comp;
    $nmax = 0 if ($nb_comp==1);

    # Get total mass
    $tot_mass = $mass_data[$n_shells][0];
    $limit_mass = 0.0;
    $i_shell = 0;
    for ($i_lag=0; $i_lag<=$#frac_lag; $i_lag++) {
      $limit_mass = $frac_lag[$i_lag]*$tot_mass;
      while ( ($mass_data[$i_shell][0] < $limit_mass) && ($mass_data[$i_shell][0] < $tot_mass) ) {
	$i_shell++;
      }
      # get masses in each component
      for ($i_comp=1; $i_comp<=$nmax; $i_comp++) {
	if ($i_shell==0) {
	  $mr_comp[$i_comp]=$mass_data[0][$i_comp];
	} elsif ($i_shell==$n_shells) {
	  $mr_comp[$i_comp]=$mass_data[$n_shells][$i_comp];
	} else { # linear interpolation
	  $M_inf=$mass_data[$i_shell-1][0];
	  $M_sup=$mass_data[$i_shell  ][0];
	  $m_inf=$mass_data[$i_shell-1][$i_comp];
	  $m_sup=$mass_data[$i_shell  ][$i_comp];
	  if ($M_sup>$M_inf) { # Linear interpolation
	    # Interpolation in the log gives smoother Lagr. radii evolution
	    $mr_comp[$i_comp] = exp(log($m_inf) + (log($m_sup)-log($m_inf))/(log($M_sup)-log($M_inf))*(log($limit_mass)-log($M_inf)));
	  } else {
	    $mr_comp[$i_comp] = $m_inf;
	  }
	}
      }
      # Determine average mass
      $nstar=0;
      $mav_lag[$i_lag] = 0;
      for ($i_comp=1; $i_comp<=$nmax; $i_comp++) {
	$mav_lag[$i_lag] = $mav_lag[$i_lag] + $mr_comp[$i_comp];
	$nstar = $nstar + $mr_comp[$i_comp]/$m_comp[$i_comp];
      }
      $mav_lag[$i_lag] = $mav_lag[$i_lag]/$nstar;
    }

    # Write Lagr. radii to a file
    # ===========================
    if ($file_started) { # Append to file
      print STDERR "> Adding data to $file_name (T=$time)...\n";
    } else { # Create file and write header
      $file_name = 'lagr_avrg_masses.asc';
      print STDERR "> Creating $file_name...\n";
      open(OUTFILE,">$file_name");
      print OUTFILE "# Lagrangian average stellar masses from SPEDI simulation\n";
      print OUTFILE "# N-body units are used\n";
      print OUTFILE "# Nstar = ",$Nstar_tot,"\n";
      print OUTFILE "# gamma_Coulomb = ",$gamma,"\n";
      print OUTFILE "# Number of mass component : $#m_comp\n";
      print OUTFILE "# 1: Time ";
      for  ($i_lag=0; $i_lag<=$#frac_lag; $i_lag++) {
	print OUTFILE $i_lag+2,": Mav",sprintf("%6.6d ",1e5*$frac_lag[$i_lag]);
      }
      print OUTFILE "\n";
      $file_started = 1;
    }
    print OUTFILE (sprintf("%12.5e " x ($#frac_lag + 2),$time,@mav_lag),"\n");

    ++$i_step;
  }
}

close(OUTFILE);
