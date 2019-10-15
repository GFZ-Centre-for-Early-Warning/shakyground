#!/bin/bash

rm -vf "test.xml"

set -x
python3 service.py quakeml_test.xml > "test.xml"
set +x

if ! [ -f "test.xml" ]; then
    echo "Output file was not created"
    exit 1
fi

# we have to shorten the grid data so xmllint has no parser problems
sed -i -e '19,500000d' "test.xml"

echo "Validating against schema ..."

if [ -x "$(command -v xmllint)" ]; then
    xmllint --noout --schema "shakemap.xsd" "test.xml"
fi

rm -vf "test.xml"

set -x
python3 service.py quakeml_test_ns.xml > "test.xml"
set +x

if ! [ -f "test.xml" ]; then
    echo "Output file was not created"
    exit 1
fi

# we have to shorten the grid data so xmllint has no parser problems
sed -i -e '19,40000d' "test.xml"

echo "Validating against schema ..."

if [ -x "$(command -v xmllint)" ]; then
    xmllint --noout --schema "shakemap.xsd" "test.xml"
fi
