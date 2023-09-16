#!/usr/bin/env bash
# Title         : Util.sh
# Description   : Defines various utility functions

print ()
{
    local IFS=" ";
    printf '%s\n' "$*"
}

print_time ()
{
    local IFS=" ";
    TIMESTAMP=$(date +%H:%M:%S)
    printf '%s | %s\n' $TIMESTAMP "$*"
}

print_labeled ()
{
    local IFS=" ";
    LABEL=$1; shift
    TIMESTAMP=$(date +%H:%M:%S)
    printf '%s | %s | %s\n' $LABEL $TIMESTAMP "$*"
}

print_line ()
{   
    [ -z "$1" ] && character="=" || character="$1"
    [ -z "$2" ] && width="100"   || width="$2"

    printf -- "$character%.s" $(seq 1 $width); echo
}

print_header ()
{
    print_line
    print "$1"
    print_line
}