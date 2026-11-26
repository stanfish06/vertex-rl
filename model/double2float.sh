#!/bin/bash

# Specify the directory where your .h files are located
directory="."

# Use the find command to locate all .h files and then use sed to replace "double" with "float"
find "$directory" -type f -name "*.h" -exec sed -i 's/double/float/g' {} \;

