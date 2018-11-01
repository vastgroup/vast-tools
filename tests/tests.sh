#!/bin/bash

set -ueo pipefail

DIRBASE=$1
DIROUT=$2
SAMPLE=NewFQ_Test_Hsa-51

# Compare with expression outcome file
cmp --silent $DIRBASE/$SAMPLE/expr_out/$SAMPLE.cRPKM $DIROUT/expr_out/$SAMPLE.cRPKM || exit 1
cmp --silent $DIRBASE/$SAMPLE/to_combine/$SAMPLE.MULTI3X $DIROUT/to_combine/$SAMPLE.MULTI3X || exit 1

exit 0
