#!/bin/bash

grep -l '\->required()' src/*.cc | xargs sed -i 's/\->required()//g'
