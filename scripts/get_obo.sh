#!/bin/bash

die() {
    echo "$@"
    exit 1
}

docker_dir="$(dirname "$0")/.."

wget -t0 'http://purl.obolibrary.org/obo/go.obo' -O "$docker_dir/extra/go.obo"
date > "$docker_dir/import_date_obo.txt"
