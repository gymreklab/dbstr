#!/bin/bash

status=`curl -s -o /dev/null -w "%{http_code}" http://webstr.ucsd.edu`
echo `date` $status >> crash_checks.log

if [ "$status" -ne "200" ]
then
       # Take any appropriate recovery action here.
	echo "webserver seems down, initiating reboot." >> check.log
	sudo reboot
fi
