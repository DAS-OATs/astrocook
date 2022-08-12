if [ -z "$1" ]
  then
    pythonw -m pytest --disable-warnings -rPx tests/*.py
else
    pythonw -m pytest --disable-warnings -rPx tests/test_$1.py
fi
