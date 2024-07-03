#!/bin/bash

base_url="https://www.ncei.noaa.gov/data/oceans/crw/5km/v3.1/nc/v1.0/daily/sst/2021/"

# Loop through the date range and download the files
start_date="20210623"
end_date="20211201"

current_date="$start_date"
while [ "$current_date" -le "$end_date" ]; do
  # Format the date to match the filenames on the website
  file_url="${base_url}CRW_SST_${current_date}.nc"

  # Download the file
  wget "$file_url"

  # Increment the date by one day (handle month/year rollovers)
  current_date=$(date -d "${current_date:0:4}-${current_date:4:2}-${current_date:6:2} + 1 day" +"%Y%m%d")
done

