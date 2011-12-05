#! /bin/bash

# Run all Mannotator tests

echo "*** TEST 1 ***"
pushd ./test1/
./test_mannotator.sh 1> /dev/null
diff mannotatored_nr.gff mannotatored_nr.gff.expected
diff mannotatored_uniref.gff mannotatored_uniref.gff.expected
rm -rf mannotator_* mannotatored_nr.gff mannotatored_uniref.gff
popd
echo ""

echo "*** TEST 2 ***"
pushd ./test2/
./test_mannotator.sh 1> /dev/null
diff mannotatored_nr.gff mannotatored_nr.gff.expected
diff mannotatored_uniref.gff mannotatored_uniref.gff.expected
rm -rf mannotator_* mannotatored_nr.gff mannotatored_uniref.gff
popd
echo ""

echo "*** TEST 3 ***"
pushd ./test3/
./test_mannotator.sh 1> /dev/null
diff mannotatored.gff3 mannotatored.gff3.expected
rm -rf mannotator_* mannotatored.gff3
popd
echo ""
