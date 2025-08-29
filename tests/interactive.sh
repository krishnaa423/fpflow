#!/bin/bash
salloc  --account=m3571 --qos=interactive --constraint=cpu --job-name=struct_job --nodes=4 --time=03:20:00
