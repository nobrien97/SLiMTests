#!/bin/bash

while read -r id; do
    rm $2/$id
done < $1