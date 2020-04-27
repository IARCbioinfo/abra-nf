#!/bin/bash
cd ~/project/

git config --global user.email "follm@iarc.fr"
git config --global user.name "Circle CI_$CIRCLE_PROJECT_REPONAME_$CIRCLE_BRANCH"
git add .
git status
git commit -m "Generated DAG [skip ci]"
git push origin $CIRCLE_BRANCH

curl -H "Content-Type: application/json" --data "{\"source_type\": \"Branch\", \"source_name\": \"$CIRCLE_BRANCH\"}" -X POST https://hub.docker.com/api/build/v1/source/1ae50eb0-017b-47b8-9411-c8af0c23758e/trigger/58ef2847-a2ad-4568-9611-4456167240a3/call/
