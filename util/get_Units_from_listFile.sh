#!/bin/sh

echo "# 0: unit_r_cgs 1: unit_m_cgs 2: unit_rho_cgs 3: unit_u_cgs 4: unit_t_cgs 5: unit_lum_cgs"
gunzip --force --to-stdout | head -n 200 | grep -A3 'CONVERSION OF CALCULATION UNITS TO CGS UNITS' \
                           | tail -n 3 | sed 's/^[^0-9]*= *//; s/ .*= */ /; s/ *[^ ]*$//' | tr '\012' ' '
echo
