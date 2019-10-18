#!/bin/bash

# Exit on any error. More complex thing could be done in future
# (see https://stackoverflow.com/questions/4381618/exit-a-script-on-error)
set -e

echo ""
echo "[INFO] Executing entrypoint..."

#---------------------
#   Save env
#---------------------
echo "[INFO] Dumping env"

# Save env vars for later usage (e.g. ssh)

env | \
while read env_var; do
  if [[ $env_var == HOME\=* ]]; then
      : # Skip HOME var
  elif [[ $env_var == PWD\=* ]]; then
      : # Skip PWD var
  else
      echo "export $env_var" >> /env.sh
  fi
done

#---------------------
#  Entrypoint command
#---------------------

if [ "$@x" == "x" ]; then
    DEFAULT_COMMAND="supervisord -c /etc/supervisor/supervisord.conf"
    echo -n "[INFO] Executing default entrypoint command: "
    echo $DEFAULT_COMMAND
    exec $DEFAULT_COMMAND
else
    echo -n "[INFO] Executing entrypoint command: "
    echo $@
    exec $@
fi 


