#!/bin/bash

set -ueo pipefail

DIRBASE=$1
DIROUT=$2

# Compare with expression outcome file
cmp --silent $DIRBASE/NewFQ_Test_Hsa-51/expr_out/NewFQ_Test_Hsa-51.cRPKM $DIROUT/expr_out/NewFQ_Test_Hsa-51.cRPKM || exit 1
cmp --silent $DIRBASE/NewFQ_Test_Hsa-51/to_combine/NewFQ_Test_Hsa-51.MULTI3X $DIROUT/to_combine/NewFQ_Test_Hsa-51.MULTI3X || exit 1

exit 0
