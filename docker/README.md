### GitLab
```
# download flups locally
git clone git@git.immc.ucl.ac.be:examples/flups.git 
# build the container
docker build -t immc/murphy-ci:v1.2 .
# push it
docker login
docker push immc/murphy-ci:v1.2
rm -rf flups
```

# local test:
```
docker run -it --rm immc/murphy-ci:v1.2
```
