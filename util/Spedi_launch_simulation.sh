#!/bin/sh
DATE=`date +%d_%m_%Y_%Hh%Mm_%Nns`
SPEDI=`which spedi`
WHICHSPEDI=`ls -l $SPEDI`

echo "Simulation started on: `date '+%d/%m/%Y %Hh%Mm %Nns'`"
echo ""
echo "Name of the output: star1.list.`date +%d_%m_%Y_%Hh%Mm_%Nns`"
echo ""
echo "Executable used:"
echo ""
echo "`echo $WHICHSPEDI`"


cat >Info_Sim_`date +%d_%m_%Y_%Hh%Mm_%Nns`.txt <<-__EOF__
Simulation started on: `date '+%d/%m/%Y %Hh%Mm %Nns'`

Name of the output: star1.list.`date +%d_%m_%Y_%Hh%Mm_%Nns`

Executable used:

`echo $WHICHSPEDI`
__EOF__

nohup $SPEDI > star1.list.`date +%d_%m_%Y_%Hh%Mm_%Nns`&
