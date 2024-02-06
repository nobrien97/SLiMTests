#!/bin/bash -l

# Purpose: Sets up a sql database for mutation data

# Path to updated version of sqlite
SQLITE3=~/Tools/sqlite/sqlite3

cat load_data.sql | $SQLITE3 standingVarMuts.db 2>/dev/null
