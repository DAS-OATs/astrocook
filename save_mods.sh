#!/bin/sh
source ~/.bashrc
source activate astrocook 2> /dev/null

file=$1
fmod=${file/.acs/_mods.acs}

cp save_mods.json save_mods_temp.json
sed -ie "s|FILE|$file|" save_mods_temp.json
sed -ie "s|FMOD|$fmod|" save_mods_temp.json
pythonw ac_gui.py save_mods_temp.json
rm save_mods_temp.json*