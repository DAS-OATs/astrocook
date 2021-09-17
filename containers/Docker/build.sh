
#!/bin/bash

if [ ! -f "Dockerfile" ]; then
    # TODO: This check is weak: improve me!
    echo "Please run this script from the containers/Docker folder"
    exit 1
fi

# Move to the root directory to allow accessing the code in Docker build context
OR_DIR=$PWD
cd ../../

# Are we on a Git repo?
echo ""
git status &> /dev/null
if [[ "x$?" == "x0" ]] ; then
    CHANGES=$(git status | grep "Changes not staged for commit") 
    if [[ "x$CHANGES" == "x" ]] ; then
        TAG=$(git rev-parse HEAD | cut -c1-7)
        echo "I will tag this container with the Git short hash \"$TAG\" "
    else
        TAG="latest"
        echo "You have uncomitted changes, I will not tag this container with the Git short hash. "
    fi
else
    TAG="latest"
fi
echo ""

# Use --no-cache in case of build problems (i.e. 404 not found)
docker build  -f containers/Docker/Dockerfile . -t astrocook -t astrocook:$TAG

cd $OR_DIR
