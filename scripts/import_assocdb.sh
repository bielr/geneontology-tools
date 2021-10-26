#!/bin/bash

die() {
    echo "$@"
    exit 1
}

docker_dir="$(dirname "$0")/.."

cd "$docker_dir"

wget -t0 'http://archive.geneontology.org/latest-full/go_monthly-assocdb-data.gz' -O- |
    gzip -d - | \
    docker-compose exec -T geneontology mysql geneontology -ugeneontology -pgeneontology

date > "$docker_dir/import_date_assocdb.txt"
