#!/bin/bash

rsync -razvhP --exclude="*exe" --exclude="flats_testing/" ./* comet:~/dist_flats
