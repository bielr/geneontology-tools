#!/bin/bash

die() {
    echo "$@"
    exit 1
}

docker_dir="$(dirname "$0")/.."

cd "$docker_dir"

wget -t0 'http://archive.geneontology.org/latest-lite/go_weekly-assocdb-data.gz' -O- |
    gzip -d - | \
    docker-compose exec -T geneontology mysql geneontology -ugeneontology -pgeneontology

date > import_date_assocdb.txt
