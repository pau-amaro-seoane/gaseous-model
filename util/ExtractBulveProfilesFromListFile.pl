#!/usr/bin/perl -w

#use strict;

use Getopt::Long;

my $time;

$i_step = 0;
$base_name = "bulve_profile";

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
  elsif (/^ *RADIUS +V-\*/) { # header line announcing data about bulk velocities for all components

    # find column numbers for various components
    # column zero is the number of the shell, column 1 is the radius
    s/V\-\*/BULVE/g;
    s/\-//g;
    s/\.//g;
    s/([A-Z]) +([0-9])/$1_$2/g;
    @headers = split;
    @columns_num = ();
    for  ($i_col=0; $i_col<=$#headers; ++$i_col) {
      if ( $headers[$i_col] =~ /^BULVE_([1-9][0-9]*)$/ ) {
	$columns_num[$1-1] = $i_col+1;
      }
    }
    $nb_comp = $#columns_num+1;
    #($nb_comp == 1) and @columns_num = $columns_num[0];

    # read block of data
    # ==================

    @bulve_data=();
    @radius=();

  SKIP_BLANK: while (<>) {
      last SKIP_BLANK if /\S/;
    }
    @data_line = split;
    $col_number = $#data_line;
    $i_shell = 0;
  DATA_READ: while (<>) {
      s/\*/ /g;
      chomp;
      @data_line = split;
      last DATA_READ if ($#data_line != $col_number);
      $radius[$i_shell] = $data_line[1  ];
      $bulve_data[$i_shell] = [ @data_line[@columns_num] ];
      $i_shell++;
    }
    $i_shell--;

    # write data to a file
    # ====================

    $file_name = $base_name . $i_step . ".asc";
    print STDOUT "T = ",$time," file : ",$file_name,"\n";
    open(OUTFILE,">$file_name");
    print OUTFILE "# Bulk velocity profile from SPEDI simulation\n";
    print OUTFILE "# N-body units are used\n";
    print OUTFILE "# Nstar = ",$Nstar_tot,"\n";
    print OUTFILE "# gamma_Coulomb = ",$gamma,"\n";
    print OUTFILE "# Time = $time\n";
    print OUTFILE "# Number of mass components = ",$nb_comp,"\n";
    print OUTFILE "# 1: Radius ";
    for  ($i_comp=1; $i_comp<=$nb_comp; ++$i_comp) {
      print OUTFILE $i_comp+1,": Bulve$i_comp ";
    }
    print OUTFILE "\n";
    for ($j_shell = 0; $j_shell<=$i_shell; $j_shell++) {
      $prov=$bulve_data[$j_shell];
      print OUTFILE $radius[$j_shell]," ",join(" ",@$prov),"\n";
    }
    close(OUTFILE);
    print "\n";
    ++$i_step;
  }
}
