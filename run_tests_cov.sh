if [ -z "$1" ]
  then
    pythonw -m pytest --cov=astrocook --cov-branch --cov-report=term-missing --disable-warnings -rPx tests/*.py
else
    pythonw -m pytest --cov=astrocook --cov-branch --cov-report=term-missing --disable-warnings -rPx tests/test_$1.py
fi
