#!/bin/bash

sed -i 's/ggiraphOutput/girafeOutput/g' ui.R
sed -i 's/renderggiraph/renderGirafe/g' server.R

exit
