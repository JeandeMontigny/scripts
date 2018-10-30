#!/bin/bash

if [ "$#" = "3" ]; then
	path=$1
	delay=$2
	name=$3
	echo "gif file construction.."
	convert -delay $delay -loop 0 $path*.png $path/$name.gif
	if [ -f $path/$name.gif ]; then
		echo "done"
	else
		echo "problem during gif creation"
	fi
else
	echo "arg error. need 3 arg: [images_path] [delay (around 15)] [gif_name]"
fi

