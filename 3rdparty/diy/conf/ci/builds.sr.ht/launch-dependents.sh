#!/bin/sh

TOKEN=`cat ~/.sr.ht.token`
AUTH=Authorization:"token ${TOKEN}"

REV=`git rev-parse HEAD`
REV_SHORT=`git rev-parse --short HEAD`

curl https://raw.githubusercontent.com/diatomic/tess2/master/.build.yml |
(cat && echo "environment: { DIY_REV: $REV }") |
jq -sR '{"manifest": .,"note":"Test tess2 with DIY rev ['$REV_SHORT'](https://github.com/diatomic/diy/commit/'$REV')"}' |
curl --silent \
    -H "$AUTH" \
    -H "Content-Type: application/json" \
    -X POST \
    -d @- \
    https://builds.sr.ht/api/jobs

GH_TOKEN=`cat ~/.github.token`
GH_AUTH=Authorization:"token ${GH_TOKEN}"
curl -H "$GH_AUTH" https://raw.githubusercontent.com/mrzv/Reeber2/master/.build.yml |
(cat && echo "environment: { DIY_REV: $REV }") |
jq -sR '{"manifest": .,"note":"Test Reeber2 with DIY rev ['$REV_SHORT'](https://github.com/diatomic/diy/commit/'$REV')"}' |
curl --silent \
    -H "$AUTH" \
    -H "Content-Type: application/json" \
    -X POST \
    -d @- \
    https://builds.sr.ht/api/jobs
