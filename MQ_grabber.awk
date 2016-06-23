#!/bin/awk -f
BEGIN { print "MQ\tGrabber" }
{ print $5}
END { print " - DONE -" }